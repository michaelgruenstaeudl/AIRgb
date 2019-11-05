*PlastomeIntegrityChecks*
=========================

Scripts for evaluating plastome integrity across all available plastid genomes on NCBI

## EXAMPLE USAGE
#### On Linux / MacOS
```
mkdir -p testing
DATE=$(date '+%Y_%m_%d')
AVAILTABLE=plastome_availability_table_${DATE}.csv

# SCRIPT 01: Generating plastome availability table
python 01_generate_plastome_availability_table.py -o testing/$AVAILTABLE 1>testing/Script01_${DATE}.runlog 2>&1

# SCRIPT 02: Downloading records and extracting IRs
mkdir -p testing/records_${DATE}
mkdir -p testing/data_${DATE}
python 02_download_records_and_extract_IRs.py -i testing/$AVAILTABLE -r testing/records_${DATE}/ -d testing/data_${DATE}/ 1>testing/Script02_${DATE}.runlog 2>&1



#python3 ../PlastomeIntegrityChecks/02_download_records_and_extract_aspects.py -l $(awk 'NR > 1 {print $2}' summaries.csv) -o .

# Compare IRs
#python3 ../PlastomeIntegrityChecks/03_compare_existing_IRs.py -o align_info.csv -i NC_*.tar.gz

```

<--
## FULL USAGE
```
#python3 ../PlastomeIntegrityChecks/01_generate_plastome_availability_table.py -o summaries.csv -q "Magnoliophyta[ORGN] AND 00000180000[SLEN] : 00000200000[SLEN] AND complete genome[TITLE] AND (chloroplast[TITLE] OR plastid[TITLE])"
```
-->

## VISUALIZATIONS & GRAPHS
#### TO DO
* Let us plot the distribution by how many nucleotides the individual records have misidentified the IR length. (It would be interesting to see if the mis-estimation is normally distributed or if there is a skew in the distribution.)
