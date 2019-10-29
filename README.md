*PlastomeIntegrityChecks*
=========================

Scripts for evaluating plastome integrity across all available plastid genomes on NCBI

## EXAMPLE USAGE
#### On Linux / MacOS
```
MASTERLIST=01_masterlist_2019_10_29.csv

mkdir -p testing

# Generating masterlist
python 01_generate_plastome_availability_table.py -o testing/$MASTERLIST

#python3 ../PlastomeIntegrityChecks/01_generate_plastome_availability_table.py -o summaries.csv -q "Magnoliophyta[ORGN] AND 00000180000[SLEN] : 00000200000[SLEN] AND complete genome[TITLE] AND (chloroplast[TITLE] OR plastid[TITLE])"

# Downloading records
python 02_download_records_and_extract_aspects.py -i testing/01_masterlist_2019_10_29.csv -o testing/records/

#python3 ../PlastomeIntegrityChecks/02_download_records_and_extract_aspects.py -l $(awk 'NR > 1 {print $2}' summaries.csv) -o .

# Compare IRs
#python3 ../PlastomeIntegrityChecks/03_compare_existing_IRs.py -o align_info.csv -i NC_*.tar.gz

```

## VISUALIZATIONS & GRAPHS
#### TO DO
* Let us plot the distribution by how many nucleotides the individual records have misidentified the IR length. (It would be interesting to see if the mis-estimation is normally distributed or if there is a skew in the distribution.)
