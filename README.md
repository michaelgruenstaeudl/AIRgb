*AIRgb*: Accessing the Inverted Repeats of plastid genomes on GenBank
=====================================================================
A Python package for automatically accessing the inverted repeats of thousands of plastid genomes stored on NCBI GenBank

## EXAMPLE USAGE
#### SCRIPT 01: Generating plastome availability table
```
# Angiosperms
TESTFOLDER=./03_testing/angiosperms_Start2000toEnd2019
DATE=$(date '+%Y_%m_%d')
MYQUERY='complete genome[TITLE] AND (chloroplast[TITLE] OR plastid[TITLE]) AND 2000/01/01:2019/12/31[PDAT] AND 0000050000:00000250000[SLEN] NOT unverified[TITLE] NOT partial[TITLE] AND (Embryophyta[ORGN] AND Magnoliophyta[ORGN])'
AVAILTABLE=plastome_availability_table_${DATE}.tsv
mkdir -p $TESTFOLDER
```
```
# Non-angiosperm landplants
TESTFOLDER=./03_testing/nonangiosperm_landplants_Start2000toEnd2019
DATE=$(date '+%Y_%m_%d')
MYQUERY='complete genome[TITLE] AND (chloroplast[TITLE] OR plastid[TITLE]) AND 2000/01/01:2019/12/31[PDAT] AND 0000050000:00000250000[SLEN] NOT unverified[TITLE] NOT partial[TITLE] AND (Embryophyta[ORGN] NOT Magnoliophyta[ORGN])'
AVAILTABLE=plastome_availability_table_${DATE}.tsv
mkdir -p $TESTFOLDER
```
```
# Defining blacklist
if [ ! -f ./02_blacklists/BLACKLIST__master_${DATE} ]; then
    cat $(ls ./02_blacklists/BLACKLIST__* | grep -v "master") > ./02_blacklists/BLACKLIST__master_${DATE}
fi
```
```
python ./01_package/01_generate_plastome_availability_table.py -q "$MYQUERY" -o $TESTFOLDER/$AVAILTABLE --blacklist ./02_blacklists/BLACKLIST__master_${DATE} 1>>$TESTFOLDER/Script01_${DATE}.runlog 2>&1
```

#### SCRIPT 02: Downloading records and extracting IR information
```
REPRTDSTAT=reported_IR_stats_table_${DATE}.tsv
mkdir -p $TESTFOLDER/records_${DATE}
mkdir -p $TESTFOLDER/data_${DATE}
```
```
python ./01_package/02_download_records_and_extract_IRs.py -i $TESTFOLDER/$AVAILTABLE -r $TESTFOLDER/records_${DATE}/ -d $TESTFOLDER/data_${DATE}/ -o $TESTFOLDER/$REPRTDSTAT 1>>$TESTFOLDER/Script02_${DATE}.runlog 2>&1
```

#### FUTURE: SCRIPTS 03 and 04
```
# SCRIPT 03: Comparing the reported IR information
#python3 ../PlastomeIntegrityChecks/02_download_records_and_extract_aspects.py -l $(awk 'NR > 1 {print $2}' summaries.tsv) -o .

# SCRIPT 04: Comparing the IRs to generate accurate IR information
#python3 ./05_future/03_compare_existing_IRs.py -o align_info.csv -i NC_*.tar.gz

```

#### OTHER COMMANDS

##### Updating blacklist on IR loss
Updating blacklist on IR loss using accompanying script
```
python generate_blacklist.py -f BLACKLIST__IR_loss__genera.txt
```

##### Identifying IR length
Self-blasting of the plastid genome sequence in order to infer the IR length
```
INF=chloroplastGenome.gb
python gb2fas.py $INF

blastn -subject ${INF%.gb*}.fas -query ${INF%.gb*}.fas -outfmt 7 -strand 'both' |\
 awk '{ if ($4 > 10000 && $4 < 50000) print $4, $7, $8, $9, $10}' >  ${INF%.gb*}_IRs.txt

for i in *.txt; do 
 sed -n '2p' $i |\
 awk '{print "IRb_start: "$2"; IRb_end: "$3"; IRa_start: "$5"; IRa_end: "$4}' > ${i%.txt*}_summary.txt; 
done
```

<!--
## FOO BAR BAZ
```
Foo bar baz
```
-->

## TO DO
#### Visualizations
* Please generate visualization code that plots the growth in plastid genomes (as done via visulations script 01), but not with the number of plastid genomes on the y-axis but with the number of unique genera plotted. (This will illustrate if the current growth in plastid genome number is reflective of a better taxonomic survey or of higher sampling in genera and families.)
