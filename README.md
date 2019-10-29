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

# Downloading records
python 02_download_records_and_extract_aspects.py -i testing/01_masterlist_2019_10_29.csv -o testing/records/

```

## VISUALIZATIONS & GRAPHS
#### TO DO
* Let us plot the distribution by how many nucleotides the individual records have misidentified the IR length. (It would be interesting to see if the mis-estimation is normally distributed or if there is a skew in the distribution.)
