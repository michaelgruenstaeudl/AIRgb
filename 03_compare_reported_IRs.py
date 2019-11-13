#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
OBJECTIVE:
    This script compares the IRs of a plastid genome and reports the number and percentage of differences, if any, between them. Specifically, the script takes a pair of IR sequences (i.e., IRa and IRb; both as separate FASTA files) and compares the two sequences via the application 'mummer' (function 'nucmer') to generate numerical indices that are indicative of the sequence similarlity.

    The output shall be a table of plastid genome records (one record per row) that lists numerically if the inverted repeats in that record are identical.


DESIGN:
    * Like in the other scripts, the evaluation of the records is conducted one by one, not all simultaneously.

    * The IR file pair is located in a spearate folder per plastid genome. The nucmer output shall be saved to the same folder.

    * A visual comparison of the IRs is not desired. Instead, the similarity of the IRs shall be inferred in a numerical fashion.

TO DO:

    * In order to know which accession numbers (and, by extension, which folders) to process, the script shall read in, and loop over, the output table of Script02. Specifically, the script shall process all those accession numbers of the first column ("ACCESSION") which have the entry "yes" for both column #2 ("IRa_REPORTED") and #6 ("IRb_REPORTED").

    * The results of this script are saved as separate files to the individual folders (named by their accession numbers) and appended as new columns to the output table of Script02 (which is read in as input to this script).

NOTES:
    * Foo bar baz
'''

#####################
# IMPORT OPERATIONS #
#####################
#import xml.etree.ElementTree as ET
import os.path, subprocess, shutil
import pandas as pd
import argparse, tarfile
import coloredlogs, logging

###############
# AUTHOR INFO #
###############
__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>, '\
             'Tilman Mehl <tilmanmehl@zedat.fu-berlin.de>'
__copyright__ = 'Copyright (C) 2019 Michael Gruenstaeudl and Tilman Mehl'
__info__ = 'Compare IRs for a series of IR FASTA files'
__version__ = '2019.11.13.1130'

#############
# DEBUGGING #
#############
#import ipdb
# ipdb.set_trace()

#############
# FUNCTIONS #
#############

def main(args):

  # STEP 1. Set up logger
    log = logging.getLogger(__name__)
    coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level='DEBUG', logger=log)



    ## TO DO: Please have this script read in a user-supplied table (i.e., the output table of script 02) and loop over all those accession numbers that have the entry "yes" for columns #2 ("IRa_REPORTED") and #6 ("IRb_REPORTED").

  # STEP 2. Loop though provided folders
    folders = [os.path.abspath(x) for x in args.input]



    ## TO DO: Please have this script append additional columns to the user-supplied table (i.e., the output table of script 02), unless these columns are already in existence (in which case the script simply continues). The column titles shall be: "MUMMER_ABS_SNPS", "MUMMER_ABS_INDELS", "MUMMER_SIMIL_SCORE", and "CONFIRM_BY_CMP".

    main_dir = os.getcwd()
    # Check if outfile exists, write headers
    # Raise an error if outfile exists but does not begin with the required headers
    if not os.path.isfile(args.outfn):
        with open(args.outfn, "w") as outfile:
            outfile.write("ACCESSION\tIR_COUNT\tLEN_A\tLEN_B\tLEN_DIFF\tSNP_COUNT\n")
    else:
        with open(args.outfn, "r") as outfile:
            if not outfile.readline() == "ACCESSION\tIR_COUNT\tLEN_A\tLEN_B\tLEN_DIFF\tSNP_COUNT\n":
                raise Exception('Malformed output file!')



    ## TO DO: Loop over the accession numbers that fullfil the above requirements.

    counter = 0
    for folder in folders:
        counter += 1



        log.info("Processing accession " + str(counter) + "/" + str(len(folders)))
        # Init values that will be written to table
        accession = os.path.basename(folder)
        ir_count = 0

        len_a = None        ## TO DO: Variable to be deleted.
        len_b = None        ## TO DO: Variable to be deleted.
        len_diff = None     ## TO DO: Variable to be deleted.
        
        snp_count = None
        # Change to directory containing sequence files
        os.chdir(folder)
        # Check if two IRs exist for this accession
        if os.path.isfile(accession + "_IRa.fasta") and os.path.isfile(accession + "_IRb_revComp.fasta"):
            log.info("Found both IRs for accession " + accession)
            ir_count = 2



            ## TO DO: Please delete the calculation of length difference here. This calculation is better done in R during the visualization process.

            # STEP 1. Inferring length difference between IRs
            with open(accession + "_IRa.fasta", "r") as ira:            ## TO DO: to be deleted
                len_a = len(ira.readlines()[1])                         ## TO DO: to be deleted
            with open(accession + "_IRb_revComp.fasta", "r") as irb:    ## TO DO: to be deleted
                len_b = len(irb.readlines()[1])                         ## TO DO: to be deleted
            len_diff = abs(len_a - len_b)                               ## TO DO: to be deleted



            ## TO DO: Please move the legacy processing steps with MUMMER into its own function (e.g., "operations_mummer_legacy(accession)" ) and put the corresponding bash code into its preamble.
            """
            def operations_mummer_legacy(accession):
            '''
            # Comparison of IR FASTAs via application 'mummer', function 'nucmer'
            nucmer --maxmatch -c 100 -p $ACCN ${ACCN}_IRa.fasta ${ACCN}_IRb_revComp.fasta -p ${ACCN}_MUMMER 1>${ACCN}_MUMMER_nucmer.log 2>&1
            show-coords -r -c -l ${ACCN}_MUMMER.delta > ${ACCN}_MUMMER.coords
            show-snps ${ACCN}_MUMMER.delta > ${ACCN}_MUMMER.snps
            show-tiling ${ACCN}_MUMMER.delta > ${ACCN}_MUMMER.tiling
            show-aligns ${ACCN}_MUMMER.delta "$ACCN"_IRa "$ACCN"_IRb_revComp > ${ACCN}_MUMMER.alignviz
            '''
                # code here
            """

            # STEP 2. Compare IRs via MUMMER
            log.info("Comparing IRs via MUMMER (function nucmer).")

            nucmerargs = ["nucmer", "--maxmatch", "-c", "100", "-p", accession+"_MUMMER", accession+"_IRa.fasta", accession+"_IRb_revComp.fasta"]
            with open(accession + "_MUMMER_nucmer.log", "w") as nucmer_log:
                nucmer = subprocess.Popen(nucmerargs, stdout=nucmer_log, sterr=subprocess.stdout)
                nucmer.wait()

            coordargs = ["show-coords", "-r", "-c", "-l", accession+"_MUMMER.delta"]
            with open(accession + "_MUMMER.coords", "w") as accession_coords:
                show_coords = subprocess.Popen(coordargs, stdout=accession_coords)
                show_coords.wait()

            snpsargs = ["show-snps", accession+"_MUMMER.delta"]
            with open(accession+"_MUMMER.snps", "w") as accession_snps:
                show_snps = subprocess.Popen(snpsargs, stdout=accession_snps)
                show_snps.wait()

            tilingargs = ["show-tiling", accession+"_MUMMER.delta"]
            with open(accession+"_MUMMER.tiling", "w") as accession_tiling:
                show_tiling = subprocess.Popen(snpsargs, stdout=accession_tiling)
                show_tiling.wait()

            alignsargs = ["show-aligns", accession+"_MUMMER.delta", accession+"_IRa", accession+"_IRb_revComp"]
            with open(accession+"_MUMMER.alignviz","w") as accession_aligns:
                show_aligns = subprocess.Popen(alignsargs, stdout=accession_aligns)
                show_aligns.wait()




            ## TO DO: Please move the new processing steps with MUMMER into its own function (e.g., "operations_mummer(accession)" ) and put the corresponding bash code into its preamble.
            """
            def operations_mummer(accession):
            '''
            # Comparison of IR FASTAs via application 'mummer', function 'dnadiff'
            dnadiff ${ACCN}_IRa.fasta ${ACCN}_IRb_revComp.fasta -p ${ACCN}_MUMMER 1>${ACCN}_MUMMER_dnadiff.log 2>&1
            '''
                # code here
            """

            # STEP 3. Compare IRs via MUMMER
            log.info("Comparing IRs via MUMMER (function dnadiff).")

            dnadiffargs = ["dnadiff", accession+"_IRa.fasta", accession+"_IRb_revComp.fasta", "-p", accession+"_MUMMER"]
            with open(accession + "_MUMMER_dnadiff.log", "w") as dnadiff_log:
                dnadiff = subprocess.Popen(dnadiffargs, stdout=dnadiff_log, sterr=subprocess.stdout)
                dnadiff.wait()




            ## TO DO: Add code here that extracts the relevant info from the output file "accession_MUMMER.report", which was generated in Step 3. Specifically, the information of the lines "AvgIdentity", "TotalSNPs" and "TotalIndels" shall be parsed out, checked for consistency (i.e., in each of the relevant lines, the two values must be identical), and saved to the correct columns (i.e., "MUMMER_SIMIL_SCORE", "MUMMER_ABS_SNPS", "MUMMER_ABS_INDELS") of the input(=output) table.
            
            # STEP 4. Parse relevant info from file "accession_MUMMER.report", check for internal consistency, and save to output file
            # code here #




            ## TO DO: Please move the processing steps with CMP into its own function (e.g., "operations_cmp(accession)" ) and put the corresponding bash code into its preamble.
            """
            def operations_cmp(accession):
            '''
            # Align the IRs
            cat "$ACCN"_IR*.fas > tmp;
            clustalo -i tmp > "$ACCN"_clustalo.fas;

            # Deinterleaving alignment (if necessary)
            # perl -MBio::SeqIO -e 'my $seqin = Bio::SeqIO->new(-fh => \*STDIN, -format => 'fasta'); while (my $seq = $seqin->next_seq) { print ">",$seq->id,"\n",$seq->seq,"\n"; }' < "$ACCN"_clustalo.fas > "$ACCN"_clustalo_deint.fas

            # Numerical comparison via command 'cmp'
            IRA_ALN=$(sed -n '2p' "$ACCN"_clustalo_deint.fas);
            IRB_ALN=$(sed -n '4p' "$ACCN"_clustalo_deint.fas);
            cmp -bl <(echo $IRA_ALN) <(echo $IRB_ALN) | awk '{print $1,$3,$5}';
            '''
                # code here
            """

            # STEP 5. Compare IRs via Bash command 'CMP'
            # Step 5.1. Align the IRs
            log.info("Performing alignment of IR FASTAs.")
            ir_filenames = [accession + "_IRa.fasta", accession + "_IRb_revComp.fasta"]
            with open("tmp","w") as tmp_file:
                for fname in ir_filenames:
                    with open(fname,"r") as infile:
                        for line in infile:
                            tmp_file.write(line)
            clustargs = ["clustalo", "-i", "tmp"]
            with open(accession + "_clustalo.fasta", "w") as clust_file:
                clustalo = subprocess.Popen(clustargs,stdout=clust_file)
                clustalo.wait()

            # Step 5.2. Deinterleave the alignment
            with open(accession + "_clustalo_deint.fasta", "w") as outfile:
              with open(accession + "_clustalo.fasta", "r") as infile:
                  outfile.write(infile.readline())
                  for line in infile:
                      line = line.strip()
                      if(len(line) > 0):
                          if line[0] == '>':
                              if "IRa" in line:
                                  outfile.write(line + "\n")
                              else:
                                  outfile.write("\n" + line + "\n")
                          else:
                              outfile.write(line)
            os.remove("tmp")
            os.remove(accession + "_clustalo.fasta")

            # Step 5.3. Numerical comparison via CMP
            log.info("Comparing differences between IR sequences using generated alignment.")
            with open(accession + ".compare","w") as comparefile:
                cmp_awk = subprocess.Popen("cmp -bl <(sed -n '2p' " + accession + "_clustalo_deint.fasta) <(sed -n '2p' " + accession + "_clustalo_deint.fasta) | awk '{print $1,$3,$5}'", shell=True, executable="/bin/bash", stdout=comparefile)
                cmp_awk.wait()
            with open(accession + ".compare","r") as comparefile:
                snp_count = len(comparefile.readlines())



            ## TO DO: Add code here that extracts the relevant info from the output file "accession_compare", which was generated in Step 5. Specifically, the information of ... shall be parsed out and compared to the information received under step 3. If the two values match, please write a "yes" into column "" of the input(=output) table.
            
            # STEP 6. Parse relevant info from file "accession_compare" and, if consistent with the MUMMER output, save to output file
            # code here #




        elif os.path.isfile(accession + "_IRa.fasta") ^ os.path.isfile(accession + "_IRb_revComp.fasta"):
            log.info("Missing one IR for accession " + accession + "! Skipping comparison.")
            ir_count = 1
        else:
            log.info("Missing both IRs for accession " + accession + "! Skipping comparison.")


        ## TO DO: The following lines are probably integrated into the steps above, as necessary.
        log.info("Writing gathered information to output.")
        os.chdir(main_dir)
        with open(args.outfn,"a") as outfile:
            outfile.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (str(accession), str(ir_count), str(len_a), str(len_b), str(len_diff), str(snp_count)))



########
# MAIN #
########

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))

    
     ## TO DO: Please have this script read in a user-supplied table (i.e., the output table of script 02).
    parser.add_argument("--input", "-i", type=str, required=True, nargs='+', help="List of folder file paths containing inverted repeat FASTA files")


    parser.add_argument("--outfn", "-o", type=str, required=True, help="Name of the file the results of the IR testing will be written to")
    args = parser.parse_args()
    main(args)
