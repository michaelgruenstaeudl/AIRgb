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
import os.path, subprocess
import pandas as pd
import argparse
import coloredlogs, logging

###############
# AUTHOR INFO #
###############
__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>, '\
             'Tilman Mehl <tilmanmehl@zedat.fu-berlin.de>'
__copyright__ = 'Copyright (C) 2019 Michael Gruenstaeudl and Tilman Mehl'
__info__ = 'Compare IRs for a series of IR FASTA files'
__version__ = '2019.11.13.1530'

#############
# DEBUGGING #
#############
#import ipdb
# ipdb.set_trace()

#############
# FUNCTIONS #
#############


def operations_mummer_legacy(accession):
    ## TO DO: Please move the legacy processing steps with MUMMER into its own function (e.g., "operations_mummer_legacy(accession)" ) and put the corresponding bash code into its preamble.
    """
    # Comparison of IR FASTAs via application 'mummer', function 'nucmer'
    # Corresponding bash code:

    nucmer --maxmatch -c 100 -p $ACCN ${ACCN}_IRa.fasta ${ACCN}_IRb_revComp.fasta -p ${ACCN}_MUMMER 1>${ACCN}_MUMMER_nucmer.log 2>&1
    show-coords -r -c -l ${ACCN}_MUMMER.delta > ${ACCN}_MUMMER.coords
    show-snps ${ACCN}_MUMMER.delta > ${ACCN}_MUMMER.snps
    show-tiling ${ACCN}_MUMMER.delta > ${ACCN}_MUMMER.tiling
    show-aligns ${ACCN}_MUMMER.delta "$ACCN"_IRa "$ACCN"_IRb_revComp > ${ACCN}_MUMMER.alignviz
    """
    #log.info("Comparing IRs via MUMMER (function nucmer).")

    nucmerargs = ["nucmer", "--maxmatch", "-c", "100", "-p", accession+"_MUMMER", accession+"_IRa.fasta", accession+"_IRb_revComp.fasta"]
    with open(accession + "_MUMMER_nucmer.log", "w") as nucmer_log:
        #nucmer = subprocess.Popen(nucmerargs, stdout=nucmer_log, sterr=subprocess.stdout) #TM: Where are we trying to pipe the error?
        nucmer = subprocess.Popen(nucmerargs, stdout=nucmer_log)
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


def operations_mummer(accession):
    '''
    # Comparison of IR FASTAs via application 'mummer', function 'dnadiff'
    # Corresponding bash code:

    dnadiff ${ACCN}_IRa.fasta ${ACCN}_IRb_revComp.fasta -p ${ACCN}_MUMMER 1>${ACCN}_MUMMER_dnadiff.log 2>&1
    '''
    dnadiffargs = ["dnadiff", accession+"_IRa.fasta", accession+"_IRb_revComp.fasta", "-p", accession+"_MUMMER"]
    with open(accession + "_MUMMER_dnadiff.log", "w") as dnadiff_log:
        #dnadiff = subprocess.Popen(dnadiffargs, stdout=dnadiff_log, sterr=subprocess.stdout)  #TM: Where are we trying to pipe the error?
        dnadiff = subprocess.Popen(dnadiffargs, stdout=dnadiff_log)
        dnadiff.wait()

def operations_cmp(accession, log):
    '''
    # Compare IRs via Bash command 'CMP'
    # Corresponding bash code:

    # Align the IRs
    cat ${ACCN}_IR*.fasta > tmp
    mafft tmp 1>${ACCN}_IR_aligned.fasta 2>${ACCN}_IR_aligned.log
    rm tmp

    # Deinterleaving alignment (if necessary)
    awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < ${ACCN}_IR_aligned.fasta | tail -n +2 > ${ACCN}_IR_aligned_deint.fasta
    rm ${ACCN}_IR_aligned.fasta

    # Numerical comparison via command 'cmp'
    IRA_ALIGNED=$(sed -n '2p' "$ACCN"_IR_aligned_deint.fasta)
    IRB_ALIGNED=$(sed -n '4p' "$ACCN"_IR_aligned_deint.fasta)
    cmp -bl <(echo $IRA_ALIGNED) <(echo $IRB_ALIGNED) > ${ACCN}_CMP.report
    awk '{print $1,$3,$5}' ${ACCN}_CMP.report
    '''

    # Align the IRs
    log.info("Aligning IRs via MAFFT")
    ir_filenames = [accession + "_IRa.fasta", accession + "_IRb_revComp.fasta"]
    with open("tmp","w") as tmp_file:
        for fname in ir_filenames:
            with open(fname,"r") as infile:
                for line in infile:
                    tmp_file.write(line)
    clustargs = ["mafft", "tmp"]
    with open(accession + "_IR_aligned.fasta", "w") as clust_file:
        IR_aligned = subprocess.Popen(clustargs,stdout=clust_file)
        IR_aligned.wait()

    # Deinterleave the alignment
    with open(accession + "_IR_aligned_deint.fasta", "w") as outfile:
      with open(accession + "_IR_aligned.fasta", "r") as infile:
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
    os.remove(accession+"_IR_aligned.fasta")

    # Numerical comparison via CMP
    log.info("Comparing aligned IRs via Bash function CMP.")
    with open(accession+"_CMP.report", "w") as comparefile:
        cmp_awk = subprocess.Popen("cmp -bl <(sed -n '2p' " + accession + "_IR_aligned_deint.fasta) <(sed -n '4p' " + accession + "_IR_aligned_deint.fasta) | awk '{print $1,$3,$5}'", shell=True, executable="/bin/bash", stdout=comparefile)
        cmp_awk.wait()
    with open(accession+"_CMP.report", "r") as comparefile:
        diff_count = len(comparefile.readlines())

    return diff_count

def main(args):

  # STEP 1. Set up logger
    log = logging.getLogger(__name__)
    coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level='DEBUG', logger=log)

  # STEP 2. Read table from Script02 and append columns that will be filled in this script
    try:
        IRinfo_file = os.path.abspath(args.infn)
        IRinfo_table = pd.read_csv(IRinfo_file, index_col=0, sep='\t', encoding="utf-8")
        relevant_accessions = list(IRinfo_table[(IRinfo_table['IRa_REPORTED']=='yes') & (IRinfo_table['IRb_REPORTED']=='yes')].index)
    except Exception as err:
        log.exception("Error reading from %s: %s" % (str(IRinfo_file), str(err)))
        raise

    if os.path.isdir(args.datadir):
        folders = [os.path.abspath(os.path.join(args.datadir, x)) for x in os.listdir(args.datadir) if x in relevant_accessions]
    else:
        log.error(args.datadir + " is not a valid directory!")
        raise Exception(args.datadir + " is not a valid directory!")

    if len(folders) < len(relevant_accessions):
        log.warning("Missing data folder for at least one provided accession.")

    added_columns = ["MUMMER_SNP_COUNT", "MUMMER_INDEL_COUNT", "MUMMER_SIMIL_SCORE", "CMP_DIFF_COUNT","CONGRUENCE_MUMMER_CMP"]
    if not any(col in list(IRinfo_table.columns) for col in added_columns):
        IRinfo_table = IRinfo_table.reindex(columns = list(IRinfo_table.columns) + added_columns)
    # Note: pd.Int64Dtype() provides a NaN-able integer data type (using regular int results in an error when trying to cast a non-finite value)
    IRinfo_table = IRinfo_table.astype({"MUMMER_SNP_COUNT": pd.Int64Dtype(), "MUMMER_INDEL_COUNT": pd.Int64Dtype(), "MUMMER_SIMIL_SCORE": float, "CMP_DIFF_COUNT": pd.Int64Dtype(), "CONGRUENCE_MUMMER_CMP": str})

    main_dir = os.getcwd()

  # STEP 3. Loop through the relevant accessions
    counter = 0
    for folder in folders:
        counter += 1

        log.info("Processing accession " + str(counter) + "/" + str(len(folders)))
        accession = os.path.basename(folder)

        # Change to directory containing sequence files
        os.chdir(folder)
        # Check if two IRs exist for this accession
        if os.path.isfile(accession + "_IRa.fasta") and os.path.isfile(accession + "_IRb_revComp.fasta"):
            log.info("Found both IRs for accession " + accession)

          # STEP 4. Compare IRs via MUMMER
            log.info("Comparing IRs via MUMMER (function dnadiff).")
            operations_mummer(accession)

          # STEP 5. Parse relevant info from file "accession_MUMMER.report", check for internal consistency, and save to output file
            try:
                with open(accession+"_MUMMER.report","r") as mumreport:
                    is_first_avgidentity = True
                    for line in mumreport.readlines():
                        if line.startswith("AvgIdentity") and is_first_avgidentity:
                            is_first_avgidentity = False
                            if line.split()[1] == line.split()[2]:
                                IRinfo_table.at[accession, "MUMMER_SIMIL_SCORE"] = float(line.split()[1])
                            else:
                                raise Exception("Values in row AvgIdentity differ!")
                        elif line.startswith("TotalSNPs"):
                            if line.split()[1] == line.split()[2]:
                                IRinfo_table.at[accession, "MUMMER_SNP_COUNT"] = int(line.split()[1])
                            else:
                                raise Exception("Values in row TotalSNPs differ!")
                        elif line.startswith("TotalIndels"):
                            if line.split()[1] == line.split()[2]:
                                IRinfo_table.at[accession, "MUMMER_INDEL_COUNT"] = int(line.split()[1])
                            else:
                                raise Exception("Values in row TotalIndels differ!")
            except Exception as err:
                log.exception("Error while parsing MUMMER report: %s\n Skipping this accession." % (err))
                os.chdir(main_dir)
                continue

          # STEP 6. Compare IRs via Bash command 'CMP'
            IRinfo_table.at[accession, "CMP_DIFF_COUNT"] = int(operations_cmp(accession, log))

          # STEP 7. Parse relevant info from file "accession_compare" and, if consistent with the MUMMER output, save to output file
            if IRinfo_table.at[accession, "CMP_DIFF_COUNT"] == IRinfo_table.at[accession, "MUMMER_SNP_COUNT"]:
                IRinfo_table.at[accession, "CONGRUENCE_MUMMER_CMP"] = "yes"
            else:
                IRinfo_table.at[accession, "CONGRUENCE_MUMMER_CMP"] = "no"
            log.info("Writing gathered information to output for accession " + accession)
            IRinfo_table.to_csv(IRinfo_file, sep='\t', header=True)
        else:
            log.warning("At least one IR fasta file is missing. Skipping this accession.")

        os.chdir(main_dir)



########
# MAIN #
########

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
     ## TO DO: Please have this script read in a user-supplied table (i.e., the output table of script 02).
    parser.add_argument("--infn", "-i", type=str, required=True, help="Name of the file that contains IR information for accessions")
    parser.add_argument("--datadir", "-d", type=str, required=True, help="path to data directory")
    #parser.add_argument("--outfn", "-o", type=str, required=True, help="Name of the file the results of the IR testing will be written to")
    args = parser.parse_args()
    main(args)
