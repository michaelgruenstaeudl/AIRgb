**Check IR identity**
=====================

Author: Michael Gruenstaeudl
Contact: m.gruenstaeudl@fu-berlin.de
Version: 2019.08.30

##### Objective
Checking if the IRs of a plastid genome in available online via NCBI are truly identical.
##### Notes
In file names with variables, don't use underscores. (Why?)

#### STEP 1. Obtain UID list for a search query
```
esearch -db nucleotide -query \
"Magnoliophyta[ORGN] AND \
00000100000[SLEN] : 00000200000[SLEN] AND \
chloroplast[TITLE] AND \
complete genome[TITLE]" | efetch -format uid > masterlist.txt
```
#### STEP 2. Obtain accession number for each UID, fetch all sequences by their UID, and check if IRs are labelled
The result is a list of plastid genomes that have either the term "inverted repeat A" or the term "inverted repeat B".

*Development: Why are the plastid ngenomes save in xml-format and not in regular gb-format?*

```
for i in $(cat masterlist.txt); do

  # Obtaining the accession number for each UID
  ACCN=$(esummary -db nucleotide -id $i | xmllint --xpath 'string(//Caption)' -);

  # Fetch all sequences by their UID; check if IRs are labelled
  # IMPORTANT: Not all cp genomes have IR labelled
  efetch -db nucleotide -id $i -format xml > $ACCN.xml;
  if grep -q "inverted repeat A\|inverted repeat B" $ACCN.xml;
    then echo $i >> filtered_list.txt
    else rm $ACCN.xml
  fi

done;
```

#### STEP 2. For each accession number, download the record, de-interleave the FASTA-file, extract the IRs, reverse-complement one of the IRs, compare the IRs via application 'mummer' (function 'nucmer'), and generate side-by-side comparison.
*Development:*
- Check if replacing the comparison via function ‘cmp’ through the application ‘mummer’ makes sense.
- Break into pieces.

```
for j in $(cat filtered_list.txt); do

  # Obtaining accession number
  ACCN=$(esummary -db nucleotide -id $j | xmllint --xpath 'string(//Caption)' -);

  # Download record and save as in FASTA and in GB format
  efetch -db nucleotide -id $j -format fasta > $ACCN.fas;
  efetch -db nucleotide -id $j -format gb > $ACCN.gb;

  # Deinterleaving FASTA file
  perl -MBio::SeqIO -e 'my $seqin = Bio::SeqIO->new(-fh => \*STDIN, -format => 'fasta'); while (my $seq = $seqin->next_seq) { print ">",$seq->id,"\n",$seq->seq,"\n"; }' < $ACCN.fas > $ACCN.deint.fas
  # OPTIONAL # perl -MBio::SeqIO -e 'my $seqin = Bio::SeqIO->new(-fh => \*STDIN, -format => 'fasta'); while (my $seq = $seqin->next_seq) { print ">",$seq->id,"\n",$seq->seq,"\n"; }' < $ACCN.fas > "$ACCN"_deint.fas

  # Extracting the start and stop information of the IRs
  IRASTART=$(grep -A10 "inverted repeat A" $ACCN.xml | grep -m 1 Seq-interval_from | xmllint --xpath 'string(//Seq-interval_from)' - );
  IRASTOP=$(grep -A10 "inverted repeat A" $ACCN.xml | grep -m 1 Seq-interval_to | xmllint --xpath 'string(//Seq-interval_to)' - );
  IRBSTART=$(grep -A10 "inverted repeat B" $ACCN.xml | grep -m 1 Seq-interval_from | xmllint --xpath 'string(//Seq-interval_from)' -);
  IRBSTOP=$(grep -A10 "inverted repeat B" $ACCN.xml | grep -m 1 Seq-interval_to | xmllint --xpath 'string(//Seq-interval_to)' -);

  # Extract DNA sequence
  FAS=$(sed -n '2p' $ACCN.deint.fas);
  # OPTINAL # FAS=$(sed -n '2p' $ACCN_deint.fas);

  # Generate FASTA only with IRa
  # OPTIONAL # echo >"$ACCN"_IRa > "$ACCN"_IRa.fas;
  # OPTIONAL # echo -n $(sed -n '1p' $ACCN.fas) > "$ACCN"_IRa.fas;
  echo ">"$ACCN"_IRa" > $ACCN.IRa.fas;
  IRA=${FAS:IRASTART:IRASTOP};
  # OPTIONAL # echo $IRA >> "$ACCN"_IRa.fas;
  # OPTIONAL # echo -n $(sed -n '1p' $ACCN.fas) > "$ACCN"_IRb.fas;
  # OPTIONAL # echo >"$ACCN"_IRb > "$ACCN"_IRb.fas;
  echo $IRA >> $ACCN.IRa.fas;

  # Generate reverse-complement FASTA only with IRb
  echo ">"$ACCN"_IRb" > $ACCN.IRb.fas;
  IRB=${FAS:IRBSTART:IRBSTOP};
  # OPTIONAL # echo $IRB >> "$ACCN"_IRb.fas;
  IRB=$(echo $IRB | rev | tr ATGC TACG);
  echo $IRB >> $ACCN.IRb.fas;
  # OPTIONAL # echo $IRB >> "$ACCN"_IRb.2.fas;

  # Comparison of IR FASTAs via application 'mummer', function 'nucmer'
  nucmer -maxmatch -c 100 -p $ACCN $ACCN.IRa.fas $ACCN.IRb.fas
  show-coords -r -c -l $ACCN.delta > $ACCN.coords
  show-snps $ACCN.delta > $ACCN.snps
  show-tiling $ACCN.delta > $ACCN.tiling

  # Generate side-by-side comparison
  show-aligns $ACCN.delta "$ACCN"_IRa "$ACCN"_IRb > $ACCN.alignviz
  # OPTIONAL # cat "$ACCN"_IR*.fas > tmp;
  # OPTIONAL # clustalo -i tmp > "$ACCN"_clustalo.fas;
  # OPTIONAL # # Deinterleaving alignment
  # OPTIONAL # perl -MBio::SeqIO -e 'my $seqin = Bio::SeqIO->new(-fh => \*STDIN, -format => 'fasta'); while (my $seq = $seqin->next_seq) { print ">",$seq->id,"\n",$seq->seq,"\n"; }' < "$ACCN"_clustalo.fas > "$ACCN"_clustalo_deint.fas
  # OPTIONAL # # Comparison via command 'cmp'
  # OPTIONAL # IRA_ALN=$(sed -n '2p' "$ACCN"_clustalo_deint.fas);
  # OPTIONAL # IRB_ALN=$(sed -n '4p' "$ACCN"_clustalo_deint.fas);
  # OPTIONAL # cmp -bl <(echo $IRA_ALN) <(echo $IRB_ALN) | awk '{print $1,$3,$5}';

  # File hygiene
  rm $ACCN.deint.fas;
  # OPTIONAL # rm $ACCN_deint.fas;
  # OPTIONAL # rm tmp;

done
```

### ALTERNATIVE: Check if the IRs in a set of input plastomes (in GFF-format) are identical.
```
for i in $(ls *.gff); do
  INF=$i

  # Extract sequence name
  ACCN=$(grep -oP '(?<=##DNA )(.*)' $INF)

  # Extract DNA sequence
  # Note: Shells don't expand variables in single quotes
  SEQ=$(cat $INF | tr '\n' ' ' | grep -oP '(?<=##DNA )(.*)(?=##end-DN)' | sed 's/ ##//g' | sed "s/$ACCN//g")

  # Extracting the start and stop information of the IRs
  # Note: The option '-m 1' is important.
  IRASTART=$(grep -m 1 "inverted repeat A" $INF | awk -F'\t' '{print $4}')
  IRASTOP=$(grep -m 1 "inverted repeat A" $INF | awk -F'\t' '{print $5}')
  IRBSTART=$(grep -m 1 "inverted repeat B" $INF | awk -F'\t' '{print $4}')
  IRBSTOP=$(grep -m 1 "inverted repeat B" $INF | awk -F'\t' '{print $5}')

  # Generate FASTA only with IRa
  echo ">"$ACCN"_IRa" > $ACCN.IRa.fas
  IRA=${SEQ:IRASTART:IRASTOP}
  echo $IRA >> $ACCN.IRa.fas

  # Generate reverse-complement FASTA only with IRb
  echo ">"$ACCN"_IRb" > $ACCN.IRb.fas
  IRB=${SEQ:IRBSTART:IRBSTOP}
  IRB=$(echo $IRB | rev | tr ATGC TACG)
  echo $IRB >> $ACCN.IRb.fas

  # Comparison of IR FASTAs via application 'mummer', function 'nucmer'
  nucmer -maxmatch -c 100 -p $ACCN $ACCN.IRa.fas $ACCN.IRb.fas
  show-coords -r -c -l $ACCN.delta > $ACCN.coords
  show-snps $ACCN.delta > $ACCN.snps
  #show-tiling $ACCN.delta > $ACCN.tiling

  # Generate side-by-side comparison
  show-aligns $ACCN.delta "$ACCN"_IRa "$ACCN"_IRb > $ACCN.alignviz

  # File hygiene
  rm $ACCN.IR?.fas
  rm $ACCN.delta
  rm $ACCN.coords
  rm $ACCN.tiling

done
```
