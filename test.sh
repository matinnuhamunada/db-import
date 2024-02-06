echo "STARTING TESTS"
(cd ../asdb-schema && bash init_database.sh)
echo ""
echo "TESTING import_json.py"
python import_json.py --taxonomy .test_data/asdb_cache.json --from-filelist .test_data/file_list.txt --db "host='localhost' port=5432 user='postgres' password='$PGPASSWORD' dbname='antismash'"
