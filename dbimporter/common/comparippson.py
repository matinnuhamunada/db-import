# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Database importing logic for CompaRiPPson results from various RiPP modules """

import antismash
from antismash.common import comparippson

from .record_data import RecordData

# cache MIBiG precursors, since there's very few of them comparatively
_MIBIG_IDS: dict[tuple[str, str], int] = {}
_RIPP_MODULES = (
    antismash.modules.lanthipeptides,
    antismash.modules.lassopeptides,
    antismash.modules.sactipeptides,
    antismash.modules.thiopeptides,
)
ASDB_SIMILARITY_THRESHOLD = 1.2
MIBIG_SIMILARITY_THRESHOLD = 0.2


def get_mibig_id(data, hit: comparippson.data_structures.Hit) -> int:
    """Get the TFBS regulator_id given the name of the regulator."""
    key = (hit.reference_fields["accession"], hit.reference_fields["locus"])
    if key not in _MIBIG_IDS:
        data.cursor.execute(
            """
SELECT comparippson_mibig_id
FROM antismash.comparippson_mibig_references
WHERE accession = ? AND name = ?""",
            key,
        )
        ret = data.cursor.fetchone()
        if ret is None:
            parameters = {
                hit.reference_fields["accession"],
                hit.reference_fields["locus"],
                hit.reference_fields["type"],
                hit.reference_fields["compounds"],
            }
            reference_id = data.insert(
                """
INSERT INTO antismash.comparippson_mibig_references (accession, name, product, compound)
VALUES (?, ?, ?, ?)
RETURNING comparippson_mibig_id""",
                parameters,
            )
        else:
            reference_id = ret[0]
        _MIBIG_IDS[key] = reference_id
    return _MIBIG_IDS[key]


def get_asdb_id(data, hit: comparippson.data_structures.Hit) -> int:
    """Get the TFBS regulator_id given the name of the regulator."""
    reference_name = hit.reference_fields["locus"]
    data.cursor.execute(
        "SELECT cds_id FROM antismash.cdss WHERE locus_tag = ?", (reference_name,)
    )
    ret = data.cursor.fetchone()
    if ret is None:
        raise ValueError(
            f"missing reference CDS from CompaRiPPson hit: {reference_name}"
        )
    cds_id = ret[0]

    data.cursor.execute(
        "SELECT comparippson_asdb_id FROM antismash.comparippson_asdb_references WHERE cds_id = ?",
        (cds_id,),
    )
    ret = data.cursor.fetchone()
    if ret is None:
        parameters = {hit.reference_fields["locus"], hit.reference_fields["type"]}
        reference_id = data.insert(
            """
INSERT INTO antismash.comparippson_asdb_references (name, product)
VALUES (?, ?)
RETURNING comparippson_asdb_id""",
            parameters,
        )
    else:
        reference_id = ret[0]
    return reference_id


def import_db(data: RecordData, results: comparippson.analysis.DBResults) -> None:
    for name, hits in results.hits.items():
        for hit in hits:
            if hit.similarity < min(
                MIBIG_SIMILARITY_THRESHOLD, ASDB_SIMILARITY_THRESHOLD
            ):
                continue
            data.cursor.execute(
                "SELECT cds_id, region_id FROM antismash.cdss WHERE locus_tag = ?",
                (name,),
            )
            cds_id, region_id = data.cursor.fetchone()
            params = {
                "cds_id": cds_id,
                "similarity": hit.similarity,
                "region_id": region_id,
                "mibig_id": None,
                "asdb_id": None,
            }
            if results.database.name.lower() == "mibig":
                if hit.similarity < MIBIG_SIMILARITY_THRESHOLD:
                    continue
                params["mibig_id"] = get_mibig_id(data, hit)
            else:
                assert results.database.name == "antiSMASH-DB", results.database.name
                if hit.similarity < ASDB_SIMILARITY_THRESHOLD:
                    continue
                params["asdb_id"] = get_asdb_id(data, hit)
            parameters = (
                params["cds_id"],
                params["similarity"],
                params["region_id"],
                params["mibig_id"],
                params["asdb_id"],
            )
            data.insert(
                """
INSERT INTO antismash.comparippson_hits (cds_id, similarity, region_id, comparippson_mibig_id, comparippson_asdb_id)
VALUES (?, ?, ?, ?, ?)""",
                parameters,
            )


def import_results(data: RecordData) -> None:
    for module in _RIPP_MODULES:
        results = data.module_results.get(module.__name__)
        if not results:
            continue
        assert isinstance(results.comparippson_results, comparippson.MultiDBResults)
        for db in results.comparippson_results.db_results:
            import_db(data, db)
