#!/bin/bash

echo "Checking for existing database"
rm antismash_db.duckdb
rm antismash_db.duckdb.wal

echo "Copying database"
cp /home/matinnu/datadrive_cemist/duckdb/antismash_db_duckdb/duckdb-schema/antismash_db.duckdb antismash_db_test.duckdb
source /home/matinnu/datadrive_cemist/duckdb/bgc_atlas_misc/.env
wget https://antismash-db.secondarymetabolites.org/output/GCF_008931305.1/GCF_008931305.1.json -nc

set -o errexit
set -o pipefail
set -o nounset

touch imported.txt
touch deferred_imported.txt

ERROR_FILE=import_errors.txt

rm -f $ERROR_FILE

IMPORTDIR=`dirname "$0"`
readonly SEQDIR=$1
readonly TAXONOMY=$2

echo "Importing base results"
for infile in $(find ${SEQDIR} -name "*.json"); do
    grep -q ${infile} imported.txt && echo "Skipping ${infile}" && continue
    echo "importing ${infile}"
    if $IMPORTDIR/import_json.py --taxonomy ${TAXONOMY} ${infile}; then
        echo ${infile} >> imported.txt
    else
        echo ${infile} >> $ERROR_FILE
        false  # marks the loop as a failure, but doesn't _exit_ the loop
    fi
done || { echo "Skipping deferred imports due to at least one base result import failure (see '$ERROR_FILE')"; exit 1;}
