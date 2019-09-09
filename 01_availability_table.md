**Availability table**
======================

* Author: Michael Gruenstaeudl
* Contact: m.gruenstaeudl@fu-berlin.de
* Version: 2019.09.09

##### Objective
The following code shall rapidly and automatically generate a table of plastid genome records that are currently available on NCBI and simultaneously ensure that no plastid genome is counted twice (issue about regular vs. RefSeq NC_ records)

The table shall contain in different columns:
1. the unique identifier,
2. the accession number (ACCESSION line),
3. the sequence version number (VERSION line),
4. the organism name (ORGANISM line, not DEFINITION line),
5. the sequence length (see LOCUS line),
6. the date the record went online (see LOCUS line),
7. the authors (see upper-most AUTHORS line),
8. the name of the publication (see upper-most TITLE line), and
9. the full citation of the publication (see upper-ost JOURNAL line)

##### Notes
* The commands `efetch -h` and `esummary -h` provide excellent background info regarding the available options.
* The COMMENT line contains the information which other record the reference sequence is identical to.
* For testing purposes (i.e., to work only on a handful of records), the start sequence length can be increased to 190000 (`00000190000[SLEN]`).
* Searches in the sequence databases of NCBI (nucleotide, protein, EST, GSS) allow the usage of [these fields](https://www.ncbi.nlm.nih.gov/books/NBK49540/).
* An example of a plastid genome record on NCBI can be found [here](https://www.ncbi.nlm.nih.gov/nucleotide/NC_031505.1).
* The list is currently restrict to angiosperms (Magnoliophyta).

##### Development / To Do

Currently, the code merely obtains the unique identifiers of a record (`-format uid`). It needs to be expanded to obtain the other eight items of information for each record. For example, efetch could fetch the full xml record (`-format xml` or maybe `-format gb -mode xml`) and code is added that parses the eight other items of information from the record; this involves xml parsing via the commandline. Alternatively, efetch could fetch the GenBank flatfile (`-format gb`) and code is added that parses the eight other items of information from the record; this involves xml parsing via the commandline.

Be aware that the list may contain thousands of records. The parsing of the records should, thus, be conducted one by one, not all simultaneously. It may, thus, be advisable to download a list of unique identifiers first (e.g., `uidlist.txt` in example) and then loop over those identifiers.

The final outcome shall be a table in which each row contains the parsed information of a single record. There mist be as many rows in that table as uids in the master list. There must be nine columns in the table (i.e., the nine items of information listed above).


##### Existing code

###### Generating uidlist
```
esearch -db nucleotide -query \
"Magnoliophyta[ORGN] AND \
00000100000[SLEN] : 00000200000[SLEN] AND \
complete genome[TITLE] AND\
(chloroplast[TITLE] OR plastid[TITLE]) \
" | efetch -db nucleotide -format uid > uidlist.txt
```

###### Parsing accession number for each UID
```
for i in $(cat uidlist.txt); do
  ACCN=$(esummary -db nucleotide -id $i | xmllint --xpath 'string(//Caption)' -);
  echo "$ACCN"
done;
```
