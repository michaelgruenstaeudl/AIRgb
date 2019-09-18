# Evaluate Published Plastid Genomes
====================================

Author: Michael Gruenstaeudl
Contact: m.gruenstaeudl@fu-berlin.de
Version: 2019.08.30

## EVALUATION 1: Does the input plastome blast with high similarity against archived plastomes of the same plant family?

#### Blasting the genome as a whole against the NCBI nt database
```
blastn -db nt -query chloroplastGenome.fasta -remote -outfmt '7 length pident sscinames' -max_target_seqs 10 -out chloroplastGenome.out

cat chloroplastGenome.out | grep -v "#" | sort -n -r -k1,1 | head -n15 | sort -n -r -k2,2 | head -n5 >> chloroplastGenome.blastn
# OPTIONAL # cat chloroplastGenome.out | grep -v "#" | sort -grk1,1 | head -n15 | sort -grk2,2 | head -n5 >> chloroplastGenome.blastn
```

## EVALUATION 2: Do the individual sections (LSC, SSC, IR) of the input plastome blast with high similarity against archived plastomes of the same plant family?

#### Split DNA sequence into 10 equally-sized files
```
INF=plastGenome.fasta
split -d -b $(bc <<< $(tail -n1 $INF | wc -c)/10) $INF prt
```
#### Blast each region against the NCBI nt database
```
for i in $(ls prt*); do
   echo $i >> results.txt;
   blastn -db nt -query $i -remote -outfmt '7 length pident sscinames' -max_target_seqs 10 -out $i.out;
   cat $i.out | grep -v "#" | sort -n -r -k1,1 | head -n15 | sort -n -r -k2,2 | head -n5 >> results.txt;
   # OPTIONAL: # cat $i.out | grep -v "#" | sort -grk1,1 | head -n15 | sort -grk2,2 | head -n5 >> results.txt;
   rm $i;
done
```

## EVALUATION 3: Is the gene order of the input genome identical (excluding one of the inverted repeats) to the gene order of archived plastomes of the same plant family?

```
# IF ONLINE # # Extract gene order of same plant family from online record
# IF ONLINE # # awk 'a !~ $1; {a=$1}' removes duplicate, adjacent lines
# IF ONLINE # esearch -db nuccore -query "SamePlantFamily [ORGN] AND 100000:200000[SLEN]" | efetch -format ft | awk -F '\t' '$4 == "gene" {print $5}' | awk 'a !~ $1; {a=$1}' > geneOrder_SamePlantFamily_online.txt

# IF LOCAL # # Extract gene order of Same Planmt Family from local GenBank file
# IF LOCAL # # awk 'a !~ $1; {a=$1}' removes duplicate, adjacent lines
# IF LOCAL # cat SamePlantFamily.gb | grep -A1 [[:space:]]gene | grep "/gene=" | awk -F '=' '{print $2}'| sed 's/\"//g' | awk 'a !~ $1; {a=$1}' > geneOrder_SamePlantFamily_local.txt
```



## EVALUATION 5: Are there occurrences where matches smaller than x bp on the plus strand and some on the minus strand are commonplace also in archived plastomes of the same plant family?

#### Downloading the chloroplast genomes of Same Plant Family from Genbank via Entrez Direct<br />Fetching chloroplast genomes
```
esearch -db nucleotide -query "Same Plant Family [ORGN] AND 00000100000[SLEN]:00000200000[SLEN] AND chloroplast[TITLE] AND complete genome[TITLE]" | efetch -format uid > masterlist.txt

for i in $(cat masterlist.txt); do
  ACCN=$(esummary -db nucleotide -id $i | xmllint --xpath 'string(//Caption)' -)
  efetch -db nucleotide -id $i -format fasta > $ACCN.fas;
  efetch -db nucleotide -id $i -format gb > $ACCN.gb;
  # OPTIONAL # efetch -db nucleotide -id $i -format xml > $ACCN.xml;
done

# Blasting Nuphar advena and extracting info on strand match
blastn -db nt -query SamePlantFamily.fas -remote -perc_identity 100 -num_descriptions 100 -ungapped | grep -A10 "^>" | grep "Strand" | sort -u

# Blasting on the minus strand of the query only
blastn -db nt -query SamePlantFamily.fas -remote -outfmt '7 length pident sscinames' -max_target_seqs 10 -strand 'minus' -out SamePlantFamily.MINUS.blastn
```

## EVALUATION 6: Generate a gene synteny table and visually compare gene syteny
Note: .gff File kommt aus Annotations Script, da dieses erst einmal ignoriert werden soll alles mit .gff File als Lagacy markieren. Da wir aber trotzdem eine gene_Order Datei haben Michael fragen ob diese dementsprechend bearbeitet werden muss.
```
# LEGACY # ## Extract gene order from .gff file
# LEGACY # # LEGACY ## LEGACY ## LEGACY ## LEGACY ## LEGACY ## LEGACY ## LEGACY ## LEGACY ## LEGACY #for INF in $(ls *.gff); do
# LEGACY #   # Extract gene order of Nuphar advena from local GenBank file
# LEGACY #   # Important: GFF file must be ordered by start position (i.e., column 4).
# LEGACY #   grep -v "#" $INF | sort -k4 -n | awk '$3=="gene" {print $9}' | awk -F ';' '{print $1}' | sed 's/Name\=//g' > ${INF%.gff*}_GeneOrder.txt
# LEGACY #   echo ${INF%.gff*} > tmp
# LEGACY #   cat ${INF%.gff*}_GeneOrder.txt >> tmp
# LEGACY #   mv tmp ${INF%.gff*}_GeneOrder.txt
# LEGACY # done

# LEGACY # ## Combine ordered gene lists and align
# LEGACY # # Combine different columns into single file; columns are sorted by file size (ls -S)
# LEGACY # paste $(ls -S *_GeneOrder.txt) > tmp
# LEGACY # # Align columns against one another by row; order is not maintained
# LEGACY # #http://stackoverflow.com/questions/34219629/align-columns-against-one-another-by-row/34220449#34220449
# LEGACY # awk 'max<NF{max=NF} {for(i=1;i<=NF;i++){b[$i]++; a[i$i]++; }} END {for(i in b) if(c<length(i)) c=length(i); for(i in b) {for(j=1;j<=max;j++){ if(a[j i]) d = i; else d=""; printf "%-"(c+5)"s", d; } print ""; } }' tmp > tmp2

# LEGACY # # Get unique words
# LEGACY # tail -n +2 tmp | sed 's/\t/\n/g' | sort | uniq | sort > tmp_f1

# LEGACY # # Get words of first column
# LEGACY # awk '{print $1}' tmp | sort > tmp_f2

# LEGACY # # Get words from unique that are not in first column
# LEGACY # comm -23 tmp_f1 tmp_f2 > tmp_f3

# LEGACY # # Generate empty file 'tmp3'
# LEGACY # touch tmp3

# LEGACY # # Reestablish original order
# LEGACY # for i in $(awk '{print $1}' tmp); do grep ^$i tmp2 >> tmp3; done

# LEGACY # # Add separator
# LEGACY # echo "##########" >> tmp3

# LEGACY # # Add lines that are not in first column;
# LEGACY # for i in $(cat tmp_f3); do grep $i tmp_f3 >> tmp3; done

# LEGACY # # Convert multiple whitespaces into single tab
# LEGACY # cat tmp3 | sed 's/ \+ /\t/g' > Overview_firstLargest.txt

# LEGACY # # File hygiene
# LEGACY # rm tmp*

# LEGACY # # Same code as above but with reference taxon selected by name, not by file size
# LEGACY # paste $(ls *_GeneOrder.txt) > tmp
# LEGACY # awk 'max<NF{max=NF} {for(i=1;i<=NF;i++){b[$i]++; a[i$i]++; }} END {for(i in b) if(c<length(i)) c=length(i); for(i in b) {for(j=1;j<=max;j++){ if(a[j i]) d = i; else d=""; printf "%-"(c+5)"s", d; } print ""; } }' tmp > tmp2
# LEGACY # tail -n +2 tmp | sed 's/\t/\n/g' | sort | uniq | sort > tmp_f1
# LEGACY # awk '{print $1}' tmp | sort > tmp_f2
# LEGACY # comm -23 tmp_f1 tmp_f2 > tmp_f3
# LEGACY # touch tmp3
# LEGACY # for i in $(awk '{print $1}' tmp); do grep ^$i tmp2 >> tmp3; done
# LEGACY # echo "##########" >> tmp3
# LEGACY # for i in $(cat tmp_f3); do grep $i tmp_f3 >> tmp3; done
# LEGACY # cat tmp3 | sed 's/ \+ /\t/g' > Overview_firstName.txt
# LEGACY # rm tmp*
```
