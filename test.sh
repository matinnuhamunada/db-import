echo "STARTING TESTS"
(cd ../asdb-schema && bash init_database.sh)
echo ""
echo "TESTING import_json.py"
python import_json.py --taxonomy .test_data/asdb_cache.json .test_data/GCF_000056065.1.json --db "host='localhost' port=5432 user='postgres' password='$PGPASSWORD' dbname='antismash'"
