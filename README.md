antiSMASH DB GenBank to SQL importer
====================================

Code to import antiSMASH results in GenBank format into the SQL database

## Usage
### Requirements
- antiSMASH 7.1.0.1
- install extra dependencies in the [`requirements.txt`](/home/matinnu/datadrive/asdb-import/requirements.txt)
- Generate Entrez API Key. see [here](https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/)
- Export the key as Environment variables
```bash
export ASDBI_ENTREZ_API_KEY=<your api key here>
```
- Create a database named `antismash` in postgreSQL
- Export the password as environment if necessary
```bash
export PGPASSWORD="<your postgres password>"
```
- Build the ncbi taxdump, follow instructions in https://github.com/matinnuhamunada/bgc_atlas_misc

License
-------

Under the same GNU AGPL v3.0 or later license as antiSMASH, see [`LICENSE`](LICENSE) file for details.
