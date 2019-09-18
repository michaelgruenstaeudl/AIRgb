## E-UTILITIES CODE

esearch -db nuccore -query "biomol genomic [PROP] AND gbdiv_pln[PROP]" |
efilter -query "plastid [ALL]" |
efilter -query "100:10000000 [SLEN]" |
efilter -query "phylogenetic study [PROP]" |
#efilter -query "phylogenetic [TITL]" |                                 # Does not work
efilter -query "2010/01/01:2010/12/31[PDAT]" |
#efilter -query "srcdb_embl[PROP]"                                      # Does not work (why?)
efilter -query "srcdb_genbank[PROP]"


esearch -db nuccore -query "biomol genomic [PROP] AND gbdiv_pln[PROP]" |\
efilter -query "plastid [ALL]" |\
efilter -query "100:10000000 [SLEN]" |\
efilter -query "phylogenetic study [PROP]" |\
#efilter -query "2010/01/01:2010/12/31[PDAT]" |
efilter -query "srcdb_genbank[PROP]" |\
efetch -format gb > test1.gb


## SCRIPT

result=""
citations=`
esearch -db nuccore -query "biomol genomic [PROP] AND gbdiv_pln[PROP]" |\
efilter -query "plastid [ALL]" |\
efilter -query "100:10000000 [SLEN]" |\
efilter -query "phylogenetic study [PROP]" |\
efilter -query "srcdb_genbank[PROP]"
`
current=`for (( yr = 2015; yr >= 1995; yr -= 1 ))
  do
    echo "$citations" |
    efilter -mindate "$yr" -maxdate "$((yr+1))" -datetype PDAT |
    xtract -pattern ENTREZ_DIRECT -lbl "${yr}" -element Count
  done`
current=`echo -e "year,count\n$current"`
if [ -n "$result" ]
then
  result=`join -t $',' <(echo "$result") <(echo "$current")`
else
  result=$current
fi
done
echo "$result" > GenBank_records.csv
