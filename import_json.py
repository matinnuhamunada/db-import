#!/usr/bin/env python3
"""Import an antiSMASH JSON results file into the antiSMASH database."""
import hashlib
import json
import logging
import os
import sys
import time
import traceback
import urllib
from argparse import ArgumentParser, FileType
from collections import defaultdict

import antismash
import duckdb
from antismash.modules.nrps_pks.data_structures import Prediction
from antismash.modules.nrps_pks.name_mappings import get_substrate_by_name
from Bio import Entrez

from dbimporter.common import preparation
from dbimporter.common.data import read_json
from dbimporter.common.record_data import RecordData
from dbimporter.modules import cluster_compare, clusterblast, pfams, tfbs, tigrfams

# pylint: disable=line-too-long,missing-docstring


log_format = "%(levelname)-8s %(asctime)s   %(message)s"
date_format = "%d/%m %H:%M:%S"
logging.basicConfig(format=log_format, datefmt=date_format, level=logging.DEBUG)

DB_CONNECTION = "antismash_db.duckdb"
Entrez.email = "matinnu@biosustain.dtu.dk"
REPORTED_TYPES = set()

options = antismash.config.build_config([])
RIPP_PRODUCTS = set(
    map(
        lambda r: r.name,
        filter(
            lambda r: r.category == "RiPP",
            antismash.detection.hmm_detection.get_ruleset(options).rules,
        ),
    )
)
assert RIPP_PRODUCTS

CANDIDATE_KINDS: dict[str, int] = {}


class ExistingRecordError(ValueError):
    pass


class MissingAssemblyIdError(ValueError):
    pass


def main(
    filename, db_connection, allowed_assembly_prefix="GCA,GCF", ignore_check=False
):
    """Run the import."""
    logging.info(f"Connecting to {db_connection}...")
    connection = duckdb.connect(db_connection)

    raw_data, results = read_json(filename)

    with connection.cursor() as cursor:
        try:
            short_name, _ = os.path.splitext(os.path.basename(filename))
            id_parts = short_name.split("_")
            allowed_assembly_prefix = [
                prefix.strip() for prefix in allowed_assembly_prefix.split(",")
            ]
            if id_parts[0] not in allowed_assembly_prefix:
                raise MissingAssemblyIdError("assembly ID does begin with 'GCF'/'GCA'")
            assembly_id = "_".join(id_parts[:2])

            logging.debug(f"Processing assembly_id: {assembly_id}")
            if assembly_id:
                input_basename = os.path.basename(filename)
                cursor.execute(
                    "SELECT (assembly_id) FROM antismash.filenames WHERE base_filename = ? AND assembly_id = ?",
                    (input_basename, assembly_id),
                )
                if cursor.fetchone() is not None:
                    raise ExistingRecordError()
                cursor.execute(
                    "INSERT INTO antismash.filenames (assembly_id, base_filename) VALUES (?, ?)",
                    (assembly_id, input_basename),
                )
            record_no = 0
            for rec, module_results in zip(results.records, results.results):
                raw_record = raw_data["records"][record_no]
                record_no += 1
                preparation.prepare_record(rec, raw_record["areas"], module_results)
                load_record(rec, module_results, cursor, assembly_id, record_no)
            connection.commit()
            logging.info(f"Changes committed for: {assembly_id}")
        except ExistingRecordError:
            connection.rollback()
        except Exception:
            connection.rollback()
            raise
    connection.close()


def load_record(rec, module_results, cur, assembly_id, record_no):
    """Load a record into the database using the cursor."""
    logging.debug(f"RUNNING: def load_record {record_no}...")
    if not rec.get_regions():
        return
    genome_id = get_or_create_genome(rec, cur, assembly_id)
    logging.info("genome_id: %s", genome_id)
    try:
        seq_id = get_or_create_dna_sequence(rec, cur, genome_id, record_no)
        logging.info("seq_id: %s", seq_id)
    except ExistingRecordError:
        print("skipping existing record:", rec.id)
        raise

    data = RecordData(cur, rec, seq_id, assembly_id, module_results, record_no)

    for region in sorted(rec.get_regions()):
        handle_region(data, seq_id, region)
        handle_ripps(data)

    add_tta_codons(data)

    for module in [cluster_compare, pfams, tfbs, tigrfams]:
        module.import_results(data)

    for gene in rec.get_genes():
        handle_gene(data, gene)


def get_or_create_dna_sequence(rec, cur, genome_id, record_no):
    """Fetch existing dna_sequence entry or create a new one."""
    logging.debug("RUNNING: def get_or_create_dna_sequence...")
    # record level entry
    params = {}
    params["seq"] = str(rec.seq)
    params["md5sum"] = hashlib.md5(params["seq"].encode("utf-8")).hexdigest()
    params["accession"] = rec.annotations["accessions"][0]
    params["version"] = rec.annotations.get("sequence_version", "0")
    params["definition"] = rec.description
    params["genome_id"] = genome_id
    params["record_number"] = record_no

    parameters = (
        params["seq"],
        params["md5sum"],
        params["accession"],
        params["version"],
        params["definition"],
        params["genome_id"],
        params["record_number"],
    )

    cur.execute(
        "INSERT INTO antismash.dna_sequences (dna, md5, accession, version, definition, genome_id, record_number)"
        "VALUES (?, ?, ?, ?, ?, ?, ?)"
        "RETURNING accession;",
        parameters,
    )
    return cur.fetchone()[0]


def get_or_create_genome(rec, cur, assembly_id):
    """Fetch existing genome entry or create a new one."""
    logging.debug("RUNNING: def get_or_create_genome...")
    # try:
    taxid = get_or_create_tax_id(
        cur, get_organism(rec), get_taxid(rec), get_strain(rec)
    )
    # except psycopg2.ProgrammingError:
    #    print(rec)
    #    raise
    cur.execute(
        "SELECT genome_id FROM antismash.genomes WHERE tax_id = ? AND assembly_id = ?",
        (taxid[0], assembly_id),
    )
    ret = cur.fetchone()
    if ret is None:
        cur.execute(
            "INSERT INTO antismash.genomes (tax_id, assembly_id) VALUES (?, ?) RETURNING genome_id, tax_id, assembly_id;",
            (taxid[0], assembly_id),
        )
        ret = cur.fetchone()

    return ret[0]


def get_taxid(rec):
    """Extract the taxid from a record."""
    logging.debug("RUNNING: def get_taxid...")
    for feature in rec.get_misc_feature_by_type("source"):
        if feature.type != "source":
            continue
        refs = feature.get_qualifier("db_xref")
        if refs is None:
            return 0
        for entry in refs:
            if entry.startswith("taxon:"):
                return int(entry[6:])
    return None


def get_organism(rec):
    """Extract the organism from a record."""
    for feature in rec.get_misc_feature_by_type("source"):
        organism = feature.get_qualifier("organism")
        if organism:
            return organism[0]

    return None


def get_strain(rec):
    """Extract the strain from a record."""
    for feature in rec.get_misc_feature_by_type("source"):
        strain = feature.get_qualifier("strain")
        if strain:
            return strain[0]
        serovar = feature.get_qualifier("serovar")
        if serovar:
            return serovar[0]

    return None


def handle_gene(data, gene):
    """Handle gene features."""
    region_id = None
    for region in data.record.get_regions():
        if gene.is_contained_by(region):
            region_id = data.feature_mapping[region]
            break
    if not region_id:
        return

    params = {
        "location": str(gene.location),
        "locus_tag": gene.locus_tag,
        "region_id": region_id,
    }
    parameters = (
        params["locus_tag"],
        params["location"],
        params["region_id"],
    )
    data.insert(
        """
INSERT INTO antismash.genes (locus_tag, location, region_id)
VALUES (?, ?, ?)""",
        parameters,
    )


def handle_cds(data, region_id, cds):
    """Handle CDS features."""
    params = {}
    params["location"] = str(cds.location)

    params["locus_tag"] = cds.locus_tag
    params["name"] = cds.gene
    params["product"] = cds.product
    params["protein_id"] = cds.protein_id
    params["func_class"] = str(cds.gene_function)
    params["translation"] = cds.translation
    params["region_id"] = region_id

    parameters = (
        params["func_class"],
        params["locus_tag"],
        params["name"],
        params["product"],
        params["protein_id"],
        params["location"],
        params["translation"],
        params["region_id"],
    )
    cds_id = data.insert(
        """
INSERT INTO antismash.cdss (
    functional_class_id,
    locus_tag,
    name,
    product,
    protein_id,
    location,
    translation,
    region_id
) VALUES (
    (SELECT functional_class_id FROM antismash.functional_classes WHERE name = ?),
    ?, ?, ?, ?, ?, ?, ?
)
RETURNING cds_id
""",
        parameters,
    )
    assert cds_id
    data.feature_mapping[cds] = cds_id

    genefunctions = data.module_results[
        antismash.detection.genefunctions.__name__
    ]._tools

    all_smcog_results = genefunctions.get("smcogs")
    hit = all_smcog_results.best_hits.get(cds.get_name())
    if hit:
        create_smcog_hit(data.cursor, hit, cds_id)

    all_resfam_results = genefunctions.get("resist")
    if all_resfam_results:
        hit = all_resfam_results.best_hits.get(cds.get_name())
        if hit:
            create_resfam_hit(data.cursor, hit, cds_id)

    create_profile_hits(data, cds)


def add_tta_codons(data):
    tta_results = data.module_results.get(antismash.modules.tta.__name__)
    if not tta_results:
        return
    for feature in tta_results.features:
        data.insert(
            "INSERT INTO antismash.tta_codons (location, seq_id) VALUES (?, ?)",
            (str(feature.location), data.record_id),
        )


def create_resfam_hit(cursor, hit, cds_id):
    cursor.execute(
        "SELECT resfam_id from antismash.resfams WHERE name = ?", (hit.hit_id,)
    )
    ret = cursor.fetchone()
    if not ret:
        raise ValueError("unknown resfam ID: %s?" % hit.hit_id)
    resfam_id = ret[0]
    cursor.execute(
        "INSERT INTO antismash.resfam_domains (score, evalue, resfam_id, cds_id) VALUES (?, ?, ?, ?)",
        (hit.bitscore, hit.evalue, resfam_id, cds_id),
    )


def create_smcog_hit(cur, hit, cds_id):
    """Create an smCOG hit entry."""
    smcog_name = hit.hit_id.split(":", 1)[0]
    smcog_score = hit.bitscore
    smcog_evalue = hit.evalue
    smcog_id = get_smcog_id(cur, smcog_name)
    cur.execute(
        "INSERT INTO antismash.smcog_hits (score, evalue, smcog_id, cds_id)"
        "VALUES (?, ?, ?, ?)",
        (smcog_score, smcog_evalue, smcog_id, cds_id),
    )


def create_profile_hits(data, cds):
    """Create profile hit entries for a feature."""
    detected_domains = parse_domains_detected(cds)
    for domain in detected_domains:
        domain["cds_id"] = data.feature_mapping[cds]
        # try:
        params = (
            domain["cds_id"],
            domain["name"],
            domain["evalue"],
            domain["bitscore"],
            domain["seeds"],
        )
        data.insert(
            """
INSERT INTO antismash.profile_hits (cds_id, name, evalue, bitscore, seeds)
VALUES (?, ?, ?, ?, ?)""",
            params,
        )
        # except psycopg2.IntegrityError:
        #    print(cds)
        #    print(domain)
        #    raise


def parse_domains_detected(feature):
    """Parse detected domains."""
    domains = []
    if not feature.sec_met:
        return domains

    for domain in feature.sec_met.domains:
        dom = {}
        dom["name"] = domain.name
        dom["evalue"] = domain.evalue
        dom["bitscore"] = domain.bitscore
        dom["seeds"] = domain.nseeds
        domains.append(dom)

    return domains


def get_smcog_id(cur, smcog):
    """Get the smcog_id given the smCOG identifier in a feature."""
    cur.execute("SELECT smcog_id FROM antismash.smcogs WHERE name = ?", (smcog,))
    ret = cur.fetchone()
    if ret is None:
        raise ValueError("could not find matching smCOG entry for: %s" % smcog)
    return ret[0]


def handle_ripps(data):
    """Handle RiPP predictions."""
    for protocluster in data.current_region.get_unique_protoclusters():
        if protocluster.product not in RIPP_PRODUCTS:
            continue
        for motif in data.record.get_cds_motifs():
            if motif.overlaps_with(protocluster) and isinstance(
                motif, antismash.common.secmet.Prepeptide
            ):
                handle_ripp(data, protocluster, motif)


def handle_ripp(data, protocluster, motif):
    """Handle a single RiPP prediction."""
    locus_tag = motif.locus_tag
    # since the motifs append which module created them so as to have unique locus tags...
    assert locus_tag.endswith("peptide"), locus_tag
    locus_tag, ripp_type = locus_tag.rsplit("_", 1)
    # multiple protoclusters can like the same motif, so make sure this is the right
    # combination
    if protocluster.product != ripp_type:
        return

    params = defaultdict(lambda: None)
    params["protocluster_id"] = data.feature_mapping[protocluster]
    params["locus_tag"] = locus_tag
    params["cds_id"] = data.feature_mapping[data.record.get_cds_by_name(locus_tag)]
    parse_ripp_core(motif, params)
    if params["peptide_sequence"] is None:
        return
    assert params["bridges"] is None or isinstance(params["bridges"], int), params
    parameters = (
        params["protocluster_id"],
        params["peptide_sequence"],
        params["molecular_weight"],
        params["monoisotopic_mass"],
        params["alternative_weights"],
        params["bridges"],
        params["class"],
        params["subclass"],
        params["score"],
        params["locus_tag"],
        params["cds_id"],
    )
    data.insert(
        """
INSERT INTO antismash.ripps (
    protocluster_id,
    peptide_sequence,
    molecular_weight,
    monoisotopic_mass,
    alternative_weights,
    bridges,
    class,
    subclass,
    score,
    locus_tag,
    cds_id
) VALUES (
    ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
)""",
        parameters,
    )


def parse_ripp_core(feature, params):
    """Parse RIPP core features."""
    params["monoisotopic_mass"] = feature.monoisotopic_mass
    params["molecular_weight"] = feature.molecular_weight
    params["alternative_weights"] = str(feature.alternative_weights)
    params["class"] = feature.peptide_class
    params["subclass"] = feature.peptide_subclass
    params["score"] = feature.score
    params["peptide_sequence"] = feature.core
    if feature.detailed_information:
        # go for a generic approach, since the detail qualifiers vary
        bridges = feature.detailed_information.to_biopython_qualifiers().get(
            "number_of_bridges", [None]
        )[0]
        if bridges is not None:
            params["bridges"] = int(bridges)


def handle_as_domain_subtype(data, as_domain_id: int, subtype: str) -> None:
    logging.debug(f"RUNNING: def handle_as_domain_subtype {as_domain_id} {subtype}...")
    """ Link domains to their subtypes """
    assert as_domain_id
    assert subtype

    # ensure the subtype is present in the subtype table
    data.cursor.execute(
        "SELECT subtype FROM antismash.as_domain_subtypes WHERE subtype = ?", (subtype,)
    )
    if not data.cursor.fetchone():
        logging.debug(f"Inserting new aSDomain subtype: {subtype}")
        data.insert(
            "INSERT INTO antismash.as_domain_subtypes (subtype, description) VALUES "
            "(?, ?)",
            (subtype, ""),
        )  # TODO add meaningful descriptions or just build this up as a static thing
    data.insert(
        "INSERT INTO antismash.rel_as_domain_to_subtype (as_domain_id, subtype) VALUES (?, ?)",
        (as_domain_id, subtype),
    )


def handle_asdomain(
    data, domain, module_id, function_id, predictions: dict[str, Prediction], follows
):
    """Handle aSDomain features."""
    params = {}
    params["pks_signature"] = None
    params["minowa"] = None
    params["nrps_predictor"] = None
    params["stachelhaus"] = None
    params["consensus"] = None
    params["kr_activity"] = None
    params["kr_stereochemistry"] = None
    params["location"] = str(domain.location)

    params["score"] = domain.score
    params["evalue"] = domain.evalue
    params["translation"] = domain.translation
    params["locus_tag"] = domain.locus_tag
    params["cds_id"] = data.feature_mapping[
        data.record.get_cds_by_name(domain.locus_tag)
    ]
    params["detection"] = domain.detection
    params["as_domain_profile_id"] = get_as_domain_profile_id(
        data.cursor, domain.domain
    )
    params["function_id"] = function_id
    params["module_id"] = module_id
    params["follows"] = follows

    parse_specificity(domain, params)

    # try:
    parameters = (
        params["detection"],
        params["score"],
        params["evalue"],
        params["translation"],
        params["pks_signature"],
        params["minowa"],
        params["nrps_predictor"],
        params["stachelhaus"],
        params["consensus"],
        params["kr_activity"],
        params["kr_stereochemistry"],
        params["as_domain_profile_id"],
        params["location"],
        params["cds_id"],
        params["module_id"],
        params["function_id"],
        params["follows"],
    )
    as_domain_id = data.insert(
        """
INSERT INTO antismash.as_domains (
    detection,
    score,
    evalue,
    translation,
    pks_signature,
    minowa,
    nrps_predictor,
    stachelhaus,
    consensus,
    kr_activity,
    kr_stereochemistry,
    as_domain_profile_id,
    location,
    cds_id,
    module_id,
    function_id,
    follows
) VALUES (
    ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
) RETURNING as_domain_id""",
        parameters,
    )
    # except psycopg2.ProgrammingError:
    #    print("error fetching cds_id for locus_tag", params["locus_tag"])
    #    raise
    data.feature_mapping[domain] = as_domain_id

    for subtype in domain.subtypes:
        handle_as_domain_subtype(data, as_domain_id, subtype)

    if params["consensus"] is None:
        return as_domain_id

    substrates = params["consensus"].split("|")
    for substrate in substrates:
        substrate_id = get_substrate(data.cursor, substrate)
        assert substrate_id is not None
        data.insert(
            "INSERT INTO antismash.rel_as_domains_substrates (as_domain_id, substrate_id) VALUES "
            "(?, ?)",
            (as_domain_id, substrate_id),
        )
    return as_domain_id


def get_as_domain_profile_id(cur, name):
    """Get the as_domain_profile_id for a domain by name."""
    if name is None:
        return None

    cur.execute(
        "SELECT as_domain_profile_id FROM antismash.as_domain_profiles WHERE name = ?",
        (name,),
    )
    ret = cur.fetchone()
    if ret is None:
        raise ValueError("Invalid asDomain name {!r}".format(name))

    return ret[0]


# TODO: Remove this once antiSMASH reports the right substrate names for PKSes
SUBSTRATE_TRANSLATIONS = {
    "ccmmal": "me-ccmal",
    "mmal": "me-mal",
    "ohmmal": "me-ohmal",
    "redmmal": "me-redmal",
}


def get_substrate(cur, name):
    """Get the substrate_id for a substrate by name."""
    name = get_substrate_by_name(SUBSTRATE_TRANSLATIONS.get(name, name)).short
    cur.execute("SELECT substrate_id FROM antismash.substrates WHERE name = ?", (name,))
    ret = cur.fetchone()
    if ret is None:
        raise ValueError("missing substrate in substrates table:", name)
    return ret[0]


def get_or_create_monomer(cur, name, substrate_id, modified):
    """Get the monomer_id for a monomer by name, create monomer entry if needed."""
    cur.execute(
        "SELECT monomer_id FROM antismash.monomers WHERE name = ?", (name.lower(),)
    )
    ret = cur.fetchone()
    if ret:
        return ret[0]
    cur.execute(
        "SELECT description FROM antismash.substrates WHERE substrate_id = ?",
        (substrate_id,),
    )
    ret = cur.fetchone()
    desc = ("modified " if modified else "") + ret[0]
    print("inserting new monomer:", name, " -> ", desc)
    cur.execute(
        "INSERT INTO antismash.monomers (substrate_id, name, description) VALUES (?, ?, ?) RETURNING monomer_id",
        (substrate_id, name.lower(), desc),
    )
    return cur.fetchone()[0]


def parse_specificity(feature, params):
    """Parse the feature's specificity entries."""
    for prediction in feature.specificity:
        method, pred = prediction.split(": ", 1)
        # drop any prediction which has no actual result
        if pred == antismash.modules.nrps_pks.results.UNKNOWN:
            continue
        if method == "KR activity":
            params["kr_activity"] = not pred.endswith("inactive")
        elif method == "KR stereochemistry":
            params["kr_stereochemistry"] = pred
        elif method == "NRPSpredictor":
            params["nrps_predictor"] = pred
        elif method == "Stachelhaus":
            params["stachelhaus"] = pred
        elif method == "Minowa":
            params["minowa"] = pred
        elif method == "PKS signature":
            params["pks_signature"] = pred
        elif method == "substrate consensus":
            params["consensus"] = pred
        elif method == "transATor":
            pass  # covered by generic domain subtype handling, but prevent the fallback error
        else:
            raise ValueError(f"unknown method: {method}")


def handle_region(data, sequence_id, region):
    """Handle cluster features."""
    assert region
    params = defaultdict(lambda: None)
    params["contig_edge"] = region.contig_edge
    params["location"] = str(region.location)
    params["start_pos"] = int(region.location.start)
    params["end_pos"] = int(region.location.end)
    params["sequence_id"] = sequence_id
    params["region_number"] = region.get_region_number()

    parameters = (
        params["region_number"],
        params["sequence_id"],
        params["location"],
        params["start_pos"],
        params["end_pos"],
        params["contig_edge"],
    )
    region_id = data.insert(
        """
INSERT INTO antismash.regions (region_number, accession, location, start_pos, end_pos, contig_edge)
VALUES (?, ?, ?, ?, ?, ?)
RETURNING region_id""",
        parameters,
    )
    params["region_id"] = region_id
    data.feature_mapping[region] = region_id
    data.current_region = region

    for product in region.products:
        product = product.lower()
        # try:
        nx_create_rel_regions_types(data.cursor, params, product)
        # except psycopg2.IntegrityError:
        #    raise RuntimeError("no definition in schema for product type: %" % product)

    for cds in region.cds_children:
        handle_cds(data, region_id, cds)

    for candidate in region.candidate_clusters:
        handle_candidate(data, candidate)
    for protocluster in region.get_unique_protoclusters():
        handle_protocluster(data, protocluster)
        if protocluster.product == "T2PKS":
            handle_t2pks(data, protocluster)
    link_proto_to_candidates(data)

    clusterblast.import_region_results(data, region, deferred=False)

    handle_region_nrpspks(data)

    return


def handle_region_nrpspks(data):
    if not hasattr(handle_region_nrpspks, "_domain_function_mapping"):
        data.cursor.execute("SELECT * FROM antismash.module_domain_functions")
        handle_region_nrpspks._domain_function_mapping = {
            func: func_id for func_id, func in data.cursor.fetchall()
        }

    function_ids = handle_region_nrpspks._domain_function_mapping

    modules = []
    domain_results = data.module_results[antismash.detection.nrps_pks_domains.__name__]
    all_nrps_pks_results = data.module_results[antismash.modules.nrps_pks.__name__]
    domain_predictions = all_nrps_pks_results.domain_predictions

    # handle all domains first, as modules will refer to them
    domains_to_id = {}
    for cds in data.current_region.cds_children:
        if cds not in domain_results.cds_results:
            continue
        previous = None
        module_id = None
        function = None
        domains = sorted(
            domain_results.cds_results[cds].domain_features.values(),
            reverse=cds.location.strand == -1,
        )
        for domain in domains:
            predictions = domain_predictions[domain.domain_id]
            previous = handle_asdomain(
                data, domain, module_id, function, predictions, follows=previous
            )
            domains_to_id[domain] = previous

    # then do a second pass, now that cross-CDS modules won't attempt to refer
    # to domains that are not yet processed
    for cds in data.current_region.cds_children:
        if not cds.modules:
            continue
        cds_domain_results = domain_results.cds_results[cds]
        cds_modules = cds.modules
        # if the first or last module is a cross-CDS module, don't reprocess it
        if (
            cds_modules[0].is_multigene_module()
            and cds.get_name() != cds_modules[0].parent_cds_names[0]
        ):
            cds_modules = cds_modules[1:]
        if (
            cds_modules
            and cds_modules[-1].is_multigene_module()
            and cds.get_name() != cds_modules[-1].parent_cds_names[0]
        ):
            cds_modules = cds_modules[:-1]
        assert len(cds_modules) == len(cds_domain_results.modules)
        for module, raw_module in zip(cds_modules, cds_domain_results.modules):
            modules.append(module)
            handle_module(
                data,
                raw_module,
                cds_domain_results,
                module,
                domains_to_id,
                function_ids,
                cds.get_name(),
            )

    if not modules:
        return

    nrps_results = data.module_results.get(antismash.modules.nrps_pks.__name__)
    if not nrps_results:
        return

    for candidate in data.current_region.candidate_clusters:
        for module in modules:
            if module.is_contained_by(candidate):
                data.insert(
                    "INSERT INTO antismash.rel_candidates_modules (candidate_id, module_id) VALUES (?, ?)",
                    (data.feature_mapping[candidate], data.feature_mapping[module]),
                )


def handle_module(
    data,
    raw_module,
    domain_results,
    secmet_module,
    domains_to_id,
    function_ids,
    cds_name,
):
    logging.debug(f"RUNNING: def handle_module {cds_name}...")
    statement = """
INSERT INTO antismash.modules (location, type, complete, iterative, region_id, trans_at, multi_gene)
VALUES (?, ?, ?, ?, ?, ?, ?)
RETURNING module_id"""
    assert raw_module.is_trans_at() in [True, False]
    values = {
        "location": str(secmet_module.location),
        "type": secmet_module.type,
        "complete": secmet_module.is_complete(),
        "iterative": secmet_module.is_iterative(),
        "region_id": data.current_region_id,
        "trans_at": raw_module.is_trans_at(),
        "multi_gene": secmet_module.is_multigene_module(),
    }
    parameters = (
        values["location"],
        values["type"],
        values["complete"],
        values["iterative"],
        values["region_id"],
        values["trans_at"],
        values["multi_gene"],
    )
    logging.debug(f"INSERTING: {parameters}")
    module_id = data.insert(statement, parameters)
    data.feature_mapping[secmet_module] = module_id

    singles = {
        "starter": raw_module._starter,
        "loader": raw_module._loader,
        "carrier_protein": raw_module._carrier_protein,
        "finalisation": raw_module._end,
    }
    # modification_components = {
    #    component.label
    #    for component in raw_module.components
    #    if component.is_modification()
    # }
    assert len(secmet_module.domains) == len(raw_module.components)
    for domain, component in zip(secmet_module.domains, raw_module.components):
        domain_id = domains_to_id[domain]
        function = "other"
        for label, single in singles.items():
            if single and domain is domain_results.domain_features.get(single.domain):
                function = label
                break
        else:
            if component.is_modification():
                function = "modification"

        update_statement = "UPDATE antismash.as_domains SET function_id = ?, module_id = ? WHERE as_domain_id = ?"
        try:
            data.cursor.execute(
                update_statement, (function_ids[function], module_id, domain_id)
            )
        except duckdb.duckdb.ConstraintException as e:
            # Handle Over-Eager Unique Constraint Checking
            logging.warning(
                f"Error updating as_domain due to Over-Eager Unique Constraint Checking in DuckDB:\n{e}"
            )
            logging.debug("Storing updated row as json for later insertion...")
            # Step 1: Query the existing row
            data.cursor.execute(
                "SELECT * FROM antismash.as_domains WHERE as_domain_id = ?",
                (domain_id,),
            )
            row = data.cursor.fetchone()
            if row:
                # Step 2: Convert the row to a dictionary
                row_dict = dict(
                    zip([column[0] for column in data.cursor.description], row)
                )

                # Step 3: Update the dictionary with new values
                row_dict["function_id"] = function_ids[function]
                row_dict["module_id"] = module_id

                # Step 4: Append the updated dictionary to a JSON file
                try:
                    with open("as_domain_update.json", "r+") as f:
                        try:
                            updated_as_domain_row = json.load(f)
                        except json.JSONDecodeError:
                            updated_as_domain_row = []

                        # Check if row_dict already exists with the same values
                        exists = any(item == row_dict for item in updated_as_domain_row)

                        if not exists:
                            updated_as_domain_row.append(row_dict)
                            f.seek(0)
                            f.truncate()  # Clear the file before writing
                            json.dump(updated_as_domain_row, f, indent=4)
                        else:
                            logging.debug(
                                "The row already exists in the file with the same values."
                            )
                except FileNotFoundError:
                    with open("as_domain_update.json", "w") as f:
                        json.dump([row_dict], f, indent=4)

    # don't insert the module if the current CDS is not the first CDS of a multi-CDS module,
    # otherwise it'll duplicate
    logging.debug(
        "Check if module is complete and CDS is the first CDS of a multi-CDS module"
    )
    if secmet_module.is_complete() and cds_name == secmet_module.parent_cds_names[0]:
        for substrate, monomer in set(secmet_module.get_substrate_monomer_pairs()):
            substrate_id = get_substrate(data.cursor, substrate)
            modified = substrate != monomer
            monomer_id = get_or_create_monomer(
                data.cursor, monomer, substrate_id, modified
            )
            statement = """
INSERT INTO antismash.rel_modules_monomers (module_id, substrate, monomer)
VALUES (?, ?, ?)"""
            data.insert(statement, (module_id, substrate_id, monomer_id))


def get_product_id(cur, product):
    cur.execute(
        "SELECT bgc_type_id FROM antismash.bgc_types WHERE term = ?", (product.lower(),)
    )
    ret = cur.fetchone()
    if ret is None:
        raise ValueError("missing product type from products table:", product)
    return ret[0]


def handle_candidate(data, candidate):
    polymer = None
    nrps_results = data.module_results.get(antismash.modules.nrps_pks.__name__)
    if nrps_results:
        for cand_result in nrps_results.region_predictions[
            data.current_region.get_region_number()
        ]:
            if (
                cand_result.candidate_cluster_number
                == candidate.get_candidate_cluster_number()
            ):
                polymer = cand_result.polymer
                break

    if candidate.kind not in CANDIDATE_KINDS:
        data.cursor.execute(
            "SELECT candidate_type_id FROM antismash.candidate_types WHERE description = ?",
            (str(candidate.kind).replace("_", " "),),
        )
        kind_id = data.cursor.fetchone()
        assert kind_id is not None
        CANDIDATE_KINDS[candidate.kind] = kind_id
    else:
        kind_id = CANDIDATE_KINDS[candidate.kind]
    candidate_id = data.insert(
        """
INSERT INTO antismash.candidates (region_id, location, candidate_type_id, smiles, polymer)
VALUES (?, ?, ?, ?, ?)
RETURNING candidate_id""",
        (
            data.current_region_id,
            str(candidate.location),
            kind_id[0],
            candidate.smiles_structure,
            polymer,
        ),
    )
    data.feature_mapping[candidate] = candidate_id

    for product in candidate.products:
        data.cursor.execute(
            """
INSERT INTO antismash.rel_candidates_types (candidate_id, bgc_type_id)
VALUES (?, ?)""",
            (candidate_id, get_product_id(data.cursor, product)),
        )


def handle_protocluster(data, protocluster):
    product_id = get_product_id(data.cursor, protocluster.product)
    protocluster_id = data.insert(
        """
INSERT INTO antismash.protoclusters (region_id, location, bgc_type_id)
VALUES (?, ?, ?)
RETURNING protocluster_id""",
        (
            data.feature_mapping[data.current_region],
            str(protocluster.location),
            product_id,
        ),
    )
    data.feature_mapping[protocluster] = protocluster_id


def link_proto_to_candidates(data):
    for candidate in data.current_region.candidate_clusters:
        cand_id = data.feature_mapping[candidate]
        for proto in candidate.protoclusters:
            proto_id = data.feature_mapping[proto]
            data.insert(
                """
INSERT INTO antismash.rel_candidates_protoclusters (candidate_id, protocluster_id)
VALUES (?, ?)""",
                (cand_id, proto_id),
            )


def nx_create_rel_regions_types(cur, params, product):
    """Create relation table to bgc_types."""
    assert params.get("region_id")
    assert product
    product_id = get_product_id(cur, product)
    cur.execute(
        """
SELECT * FROM antismash.rel_regions_types WHERE region_id = ? AND
    bgc_type_id = ?""",
        (params["region_id"], product_id),
    )
    ret = cur.fetchone()
    if ret is None:
        cur.execute(
            """
INSERT INTO antismash.rel_regions_types (region_id, bgc_type_id)
VALUES (?, ?)""",
            (params["region_id"], product_id),
        )


def handle_t2pks(data, protocluster):
    protocluster_id = data.feature_mapping[protocluster]
    t2pks_id = data.insert(
        "INSERT INTO antismash.t2pks (protocluster_id) VALUES (?) RETURNING t2pks_id",
        (protocluster_id,),
    )

    t2pks_results = data.module_results[
        antismash.modules.t2pks.__name__
    ].cluster_predictions[protocluster.get_protocluster_number()]

    for starter in t2pks_results.starter_units:
        statement = """
INSERT INTO antismash.t2pks_starters (t2pks_id, name, evalue, bitscore)
VALUES (?, ?, ?, ?)
RETURNING domain_id"""
        values = {
            "t2pks_id": t2pks_id,
            "name": starter.name,
            "evalue": starter.evalue,
            "bitscore": starter.score,
        }
        parameters = (
            values["t2pks_id"],
            values["name"],
            values["evalue"],
            values["bitscore"],
        )
        starter_id = data.insert(statement, parameters)
        for elongation_joined in t2pks_results.malonyl_elongations:
            for elongation in elongation_joined.name.split("|"):
                try:
                    weight = t2pks_results.molecular_weights[
                        "%s_%s" % (starter.name, elongation)
                    ]
                except KeyError:
                    print(t2pks_results.molecular_weights)
                    raise
                elongation_count = int(elongation)
                statement = """
    INSERT INTO antismash.t2pks_starter_elongation (domain_id, elongation, weight)
    VALUES (?, ?, ?)
    RETURNING combo_id"""
                data.insert(statement, (starter_id, elongation_count, weight))

    for cds_name, cds_results in t2pks_results.cds_predictions.items():
        cds_id = data.feature_mapping[data.record.get_cds_by_name(cds_name)]
        assert cds_id, cds_id
        for cds_result in cds_results:
            data.cursor.execute(
                "SELECT profile_id FROM antismash.t2pks_profiles WHERE name = ?",
                (cds_result.ptype,),
            )
            profile_id = data.cursor.fetchone()[0]

            statement = """
INSERT INTO antismash.t2pks_cds_domain (t2pks_id, cds_id, profile_id, protein_type, protein_function, evalue, bitscore)
VALUES (?, ?, ?, ?, ?, ?, ?)
RETURNING domain_id"""
            values = {
                "t2pks_id": t2pks_id,
                "cds_id": cds_id,
                "profile_id": profile_id,
                "ptype": cds_result.ptype,
                "pfunc": cds_result.pfunc,
                "evalue": cds_result.evalue,
                "score": cds_result.bitscore,
            }
            parameters = (
                values["t2pks_id"],
                values["cds_id"],
                values["profile_id"],
                values["ptype"],
                values["pfunc"],
                values["evalue"],
                values["score"],
            )
            data.insert(statement, parameters)

    for product in sorted(t2pks_results.product_classes):
        statement = """
INSERT INTO antismash.t2pks_product_classes (t2pks_id, product_class)
VALUES (?, ?)"""
        data.insert(statement, (t2pks_id, product))


def get_or_create_tax_id(cur, name, ncbi_taxid, strain):
    """Get the tax_id or create a new one."""
    logging.debug("RUNNING: def get_or_create_tax_id...")
    combined_name = f"{name} {strain}"
    cur.execute(
        "SELECT tax_id FROM antismash.taxa WHERE name = ? AND ncbi_taxid = ?",
        (combined_name, ncbi_taxid),
    )
    ret = cur.fetchone()
    if ret is None:
        lineage = get_lineage(ncbi_taxid)
        lineage["ncbi_taxid"] = ncbi_taxid
        lineage["name"] = combined_name
        lineage["strain"] = strain
        try:
            parameters = (
                lineage["ncbi_taxid"],
                lineage["superkingdom"],
                lineage["kingdom"],
                lineage["phylum"],
                lineage["class"],
                lineage["order"],
                lineage["family"],
                lineage["genus"],
                lineage["species"],
                lineage["strain"],
                lineage["name"],
            )
            cur.execute(
                """
INSERT INTO antismash.taxa (ncbi_taxid, superkingdom, kingdom, phylum, class, taxonomic_order, family, genus, species, strain, name) VALUES
    (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?) RETURNING tax_id""",
                parameters,
            )
        except KeyError:
            print("Error inserting {!r}".format(lineage))
            raise
        ret = cur.fetchone()
    return ret


def get_lineage(taxid):
    logging.debug("RUNNING: def get_lineage...")
    """Get the full lineage for a taxid from Entrez."""
    info = TAX_DUMP["mappings"].get(
        str(TAX_DUMP["deprecated_ids"].get(str(taxid), taxid))
    )
    if info:
        result = dict(info)
        result.pop("tax_id")
        return result

    api_key = os.environ.get("ASDBI_ENTREZ_API_KEY", "")
    extra_params = {}
    if api_key:
        extra_params["api_key"] = api_key
    retries = 5
    while retries > 0:
        try:
            handle = Entrez.efetch(
                db="taxonomy", id=taxid, retmode="xml", **extra_params
            )
            break
        except urllib.error.HTTPError as err:
            retries -= 1
            print("Got http error:", err, retries, "left")

    records = Entrez.read(handle)
    lineage = defaultdict(lambda: "Unclassified")
    for entry in records[0]["LineageEx"]:
        if entry["Rank"] == "no rank":
            continue
        lineage[entry["Rank"]] = entry["ScientificName"].split(" ")[-1]

    if (
        "species" not in lineage
        and "ScientificName" in records[0]
        and len(records[0]["ScientificName"].split()) == 2
    ):
        lineage["species"] = records[0]["ScientificName"].split()[-1]
    return lineage


def test_delete():
    connection = duckdb.connect(DB_CONNECTION)
    with connection.cursor() as cursor:
        cursor.execute("SELECT * FROM antismash.genomes")
        existing = cursor.fetchall()
        if existing:
            print("existing genomes")
            for row in existing:
                print(row)
            print("deleting genomes")
            cursor.execute("DELETE FROM antismash.genomes RETURNING *")
            rows = cursor.fetchall()
            if rows:
                print("rows deleted")
                print(rows)
            connection.commit()
        else:
            print("no existing genomes")

        print("tables with entries remaining:")
        tables = sorted(
            [
                "filenames",
                "t2pks_starter_elongation",
                "genomes",
                "t2pks_starters",
                "t2pks",
                "t2pks_product_classes",
                "dna_sequences",
                "ripps",
                "t2pks_cds_domain",
                "monomers",
                "tta_codons",
                "rel_modules_monomers",
                "regions",
                "module_modifications",
                "rel_as_domains_substrates",
                "protoclusters",
                "rel_regions_types",
                "candidates",
                "rel_candidates_protoclusters",
                "rel_candidates_types",
                "modules",
                "rel_candidates_modules",
                "clusterblast_hits",
                "profile_hits",
                "as_domains",
                "smcog_hits",
                "cdss",
                "genes",
                "pfam_domains",
                "pfam_go_entries",
            ]
        )
        for table in tables:
            cursor.execute("SELECT COUNT(*) FROM antismash.?" % table)
            count = cursor.fetchone()[0]
            if count > 0:
                print("", table, count)

        connection.commit()
    connection.close()


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument(
        "--db",
        default=DB_CONNECTION,
        help="DB connection string to use (default: %(default)s)",
    )
    parser.add_argument("--taxonomy", type=FileType("r"), help="Taxonomy dump as JSON")
    parser.add_argument(
        "--from-filelist",
        action="store_true",
        default=False,
        help="Passed filename is a list of filenames",
    )
    parser.add_argument(
        "--success-log",
        type=FileType("a", encoding="utf-8"),
        default=None,
        help="File to store successfully imported file names in",
    )
    parser.add_argument(
        "--error-log",
        type=FileType("a", encoding="utf-8"),
        default=None,
        help="File to store file names that failed to import in",
    )
    parser.add_argument(
        "--allowed_assembly_prefix",
        default="GCF,GCA",
        help="Allowed assembly prefixes (default: %(default)s)",
    )
    parser.add_argument(
        "--ignore_check",
        action="store_true",
        default=False,
        help="Ignore check (default: %(default)s)",
    )
    parser.add_argument("filenames", nargs="*")
    args = parser.parse_args()

    TAX_DUMP = (
        json.load(args.taxonomy)
        if args.taxonomy
        else {"mappings": {}, "deprecated_ids": {}}
    )

    total_duration = 0
    total_imports = 0
    failed = False
    successful_imports = 0

    filenames: list[str] = args.filenames
    if args.from_filelist:
        filenames = []
        for filename in args.filenames:
            with open(filename, "r", encoding="utf-8") as handle:
                content = handle.read()
                filenames.extend(content.splitlines())

    for filename in filenames:
        start_time = time.time()
        try:
            main(
                filename,
                args.db,
                allowed_assembly_prefix=args.allowed_assembly_prefix,
                ignore_check=args.ignore_check,
            )
            successful_imports += 1
            if args.success_log:
                print(filename, file=args.success_log)
        except MissingAssemblyIdError as err:
            print("failed to import", filename, ":", err)
            failed = True
            if args.error_log:
                print(filename, file=args.error_log)
        except Exception as err:
            print("failed to import", filename, ":", err)
            traceback.print_exc()
            failed = True
            if args.error_log:
                print(filename, file=args.error_log)
        finally:
            end_time = time.time()
            import_duration = end_time - start_time
            print("took", round(import_duration, 2), "seconds", end="\t")
            total_imports += 1
            total_duration += import_duration
            print(
                "average:",
                round(total_duration / total_imports, 2),
                f"for {total_imports} total imports ({successful_imports} successful)",
            )
    sys.exit(1 if failed else 0)
