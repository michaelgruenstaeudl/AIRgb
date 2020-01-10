#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
OBJECTIVE:
    This script takes a complete plastid genome sequence (in FASTA format) and re-calculates (re-infers) the position (and, thus, the length) of the IRs.

    The output shall be a table of plastid genome records (one record per row) that lists the originally inferred IR length and this newly calculated IR length so that a comparison is possible.

TO DO:
    * Once the IRs are re-calculated, the originally inferred IR length and the newly calculated IR length shall be compared in order to see if previous studies have - on average - overestimated or underestimated the IR length.

    * If differences between the originally inferred IR length and the newly calculated IR length are discovered, it will be interesting to see on which side of the IRs (the side that faces the LSC or the side that faces the SSC) the original inference was incorrect (i.e., on which side a bias in the original inference happened).

    * The following batch code shall be used to re-calculate the compare the IRs (i.e., the IR file pair) of each record:
        ```
        # Self-blasting of the plastid genome sequence in order to infer the IR length
        blastn -db chloroplastGenome.fasta -query chloroplastGenome.fasta -outfmt 7 -strand 'both' | awk '{ if ($4 > 10000 && $4 < 50000) print $4, $7, $8, $9, $10}'
        ```

    * Also do the inference of the IR via MUMMER (specifically, a self-comparison via dnadiff) so that the IR boundaries as inferred via self-BLASTING are confirmed (i.e., similar to the internal confirmation check of the total number of sequence differences via CMP).

DESIGN:
    * Like in the other scripts, the evaluation of the records is conducted one by one, not all simultaneously.

    * The plastid genome sequence of each record will be bundled with the GB file of that record in a record-specific gzip file. Hence, before that record is evaluated here, that specific gzip-file must be unpacked; after the evaluation, the gzip-file must be reconstituted.

NOTES:
    * Start positions for IRa and IRb are sometimes switched around (probably because they were improperly recorded in previous scripts due to missing unique identification)
'''

#####################
# IMPORT OPERATIONS #
#####################
#import xml.etree.ElementTree as ET
import os.path, subprocess
import argparse, tarfile
import pandas as pd
import coloredlogs, logging

###############
# AUTHOR INFO #
###############
__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>, '\
             'Tilman Mehl <tilmanmehl@zedat.fu-berlin.de>'
__copyright__ = 'Copyright (C) 2019 Michael Gruenstaeudl and Tilman Mehl'
__info__ = 'Re-calculates the position and length of the IRs for each '\
           'plastid genome record'
__version__ = '2019.10.17.1530'

#############
# DEBUGGING #
#############
import ipdb
# ipdb.set_trace()


#############
# FUNCTIONS #
#############
def main(args):
    # STEP 1. Set up logger
    log = logging.getLogger(__name__)
    coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level='DEBUG', logger=log)

    # STEP 2. Read table from Script02/Script03 and append columns that will be filled in this script
    try:
        IRinfo_file = os.path.abspath(args.repfn)
        IRinfo_table = pd.read_csv(IRinfo_file, index_col=0, sep='\t', encoding="utf-8")
        relevant_accessions = list(IRinfo_table[(IRinfo_table['IRa_REPORTED']=='yes') & (IRinfo_table['IRb_REPORTED']=='yes')].index)
    except Exception as err:
        log.exception("Error reading from %s: %s" % (str(IRinfo_file), str(err)))
        raise

    if os.path.isdir(args.data):
        folders = [os.path.abspath(os.path.join(args.data, x)) for x in os.listdir(args.data) if x in relevant_accessions]
    else:
        log.error(args.data + " is not a valid directory!")
        raise Exception(args.data + " is not a valid directory!")

    if len(folders) < len(relevant_accessions):
        log.warning("Missing data folder for at least one provided accession.")

    added_columns = ["IRa_CALCULATED", "IRa_CALCULATED_START", "IRa_CALCULATED_END", "IRa_CALCULATED_LENGTH", "IRa_START_COMPARED_OFFSET", "IRa_END_COMPARED_OFFSET", "IRa_LENGTH_COMPARED_DIFFERENCE", "IRb_CALCULATED", "IRb_CALCULATED_START", "IRb_CALCULATED_END", "IRb_CALCULATED_LENGTH", "IRb_START_COMPARED_OFFSET", "IRb_END_COMPARED_OFFSET", "IRb_LENGTH_COMPARED_DIFFERENCE"]
    if not any(col in list(IRinfo_table.columns) for col in added_columns):
        IRinfo_table = IRinfo_table.reindex(columns = list(IRinfo_table.columns) + added_columns)

    #IRinfo_table = IRinfo_table.astype({"IRa_CALCULATED": str, "IRa_CALCULATED_START": pd.Int64Dtype(), "IRa_CALCULATED_END": pd.Int64Dtype(), "IRa_CALCULATED_LENGTH": pd.Int64Dtype(), "IRa_START_COMPARED_OFFSET": pd.Int64Dtype(), "IRa_END_COMPARED_OFFSET": pd.Int64Dtype(), "IRa_LENGTH_COMPARED_DIFFERENCE": pd.Int64Dtype(), "IRb_CALCULATED": str, "IRb_CALCULATED_START": pd.Int64Dtype(), "IRb_CALCULATED_END": pd.Int64Dtype(), "IRb_CALCULATED_LENGTH": pd.Int64Dtype(), "IRb_START_COMPARED_OFFSET": pd.Int64Dtype(), "IRb_END_COMPARED_OFFSET": pd.Int64Dtype(), "IRb_LENGTH_COMPARED_DIFFERENCE": pd.Int64Dtype(), "MUMMER_SNP_COUNT": pd.Int64Dtype(), "MUMMER_INDEL_COUNT": pd.Int64Dtype(), "MUMMER_SIMIL_SCORE": float, "CMP_DIFF_COUNT": pd.Int64Dtype(), "CONGRUENCE_MUMMER_CMP": str})

    IRinfo_table = IRinfo_table.astype({"IRa_CALCULATED": str, "IRb_CALCULATED": str})

    '''
    #folders = [os.path.abspath(x) for x in args.data]
    # Check if outfile exists, write headers
    # Raise an error if outfile exists but does not begin with the required headers
    outfn = os.path.abspath(args.outfn)
    if not os.path.isfile(outfn):
        with open(outfn, "w") as outfile:
            outfile.write("ACCESSION\tSTART_A_ORIG\tSTART_A\tLEN_A_ORIG\tLEN_A\tLEN_A_DIFF\tSTART_B_ORIG\tSTART_B\tLEN_B_ORIG\tLEN_B\tLEN_B_DIFF\n")
    else:
        with open(outfn, "r") as outfile:
            if not outfile.readline() == "ACCESSION\tSTART_A_ORIG\tSTART_A\tLEN_A_ORIG\tLEN_A\tLEN_A_DIFF\tSTART_B_ORIG\tSTART_B\tLEN_B_ORIG\tLEN_B\tLEN_B_DIFF\n":
                raise Exception('Malformed output file!')

    reportfile = os.path.abspath(args.repfn)
    if os.path.isfile(reportfile):
        with open(reportfile) as rep:
            header = rep.readline()
            if not header == "\t".join(["ACCESSION", "IRa_REPORTED", "IRa_REPORTED_START", "IRa_REPORTED_END", "IRa_REPORTED_LENGTH", "IRb_REPORTED", "IRb_REPORTED_START", "IRb_REPORTED_END", "IRb_REPORTED_LENGTH"]) + "\n":
                raise Exception("Missing or malformed header in file " + args.repfn)
        ir_reported = pd.read_csv(reportfile, sep='\t', index_col=0, encoding="utf-8")
    else:
        raise Exception("File not found: " + args.repfn)
    '''

    main_dir = os.getcwd()
    for folder in folders:
        # Init values that will be written to table
        accession = os.path.basename(folder)

        '''
        start_a_orig = None
        start_b_orig = None
        len_a_orig = 0
        len_b_orig = 0
        len_a = 0
        len_b = 0
        start_a = None
        start_b = None

        try:
            ir_acc = ir_reported.loc[accession]
        except:
            log.warning("Could not find any reported IR information for accession " + accession + ". Skipping this accession.")
            continue

        if not ir_acc is None:
            try:
                start_a_orig = int(ir_acc["IRa_REPORTED_START"])
            except:
                start_a_orig = ir_acc["IRa_REPORTED_START"]
            try:
                len_a_orig = int(ir_acc["IRa_REPORTED_LENGTH"])
            except:
                len_a_orig = 0
            try:
                start_b_orig = int(ir_acc["IRb_REPORTED_START"])
            except:
                start_b_orig = ir_acc["IRb_REPORTED_START"]
            try:
                len_b_orig = int(ir_acc["IRb_REPORTED_LENGTH"])
            except:
                len_b_orig = 0
        '''
        # Change to directory containing sequence files
        os.chdir(folder)
        # Calculate IR positions and length
        try:
            mkblastargs = ["makeblastdb", "-in", accession + "_completeSeq.fasta", "-parse_seqids", "-title", accession, "-dbtype", "nucl"]
            mkblastdb_subp = subprocess.Popen(mkblastargs)
            mkblastdb_subp.wait()
        except Exception as err:
            log.exception("Error creating blastdb: %s\nSkipping this accession." % (str(err)))
            continue
        try:
            blastargs = ["blastn", "-db", str(accession) + "_completeSeq.fasta", "-query", accession + "_completeSeq.fasta", "-outfmt", "7", "-strand", "both"]
            blast_subp = subprocess.Popen(blastargs, stdout=subprocess.PIPE)
            awkargs = ["awk", "{if ($4 > 10000 && $4 < 50000) print $4, $7, $8, $9, $10}"]
            awk_subp = subprocess.Popen(awkargs, stdin=blast_subp.stdout, stdout=subprocess.PIPE)
            out, err = awk_subp.communicate()
            if len(out.splitlines()) == 2:
                # TODO: Look up and verify if the first line is indeed always IRb and the second line IRa
                ira_info = out.splitlines()[1].split()
                irb_info = out.splitlines()[0].split()

                IRinfo_table.at[accession, "IRa_CALCULATED"] = "yes"
                IRinfo_table.at[accession, "IRb_CALCULATED"] = "yes"

                IRinfo_table.at[accession, "IRa_CALCULATED_START"] = int(ira_info[1])
                IRinfo_table.at[accession, "IRb_CALCULATED_START"] = int(irb_info[1])

                ''' # TODO: Check which index is for end position
                IRinfo_table.at[accession, "IRa_CALCULATED_END"] = int(ira_info[2])
                IRinfo_table.at[accession, "IRb_CALCULATED_END"] = int(irb_info[2])
                '''
                IRinfo_table.at[accession, "IRa_CALCULATED_LENGTH"] = int(ira_info[0])
                IRinfo_table.at[accession, "IRb_CALCULATED_LENGTH"] = int(irb_info[0])

                IRinfo_table.at[accession, "IRa_START_COMPARED_OFFSET"] = int(IRinfo_table.at[accession, "IRa_REPORTED_START"]) - int(IRinfo_table.at[accession, "IRa_CALCULATED_START"])
                IRinfo_table.at[accession, "IRb_START_COMPARED_OFFSET"] = int(IRinfo_table.at[accession, "IRb_REPORTED_START"]) -int(IRinfo_table.at[accession, "IRb_CALCULATED_START"])

                IRinfo_table.at[accession, "IRa_END_COMPARED_OFFSET"] = int(IRinfo_table.at[accession, "IRa_REPORTED_END"]) - int(IRinfo_table.at[accession, "IRa_CALCULATED_END"])
                IRinfo_table.at[accession, "IRb_END_COMPARED_OFFSET"] = int(IRinfo_table.at[accession, "IRb_REPORTED_END"]) - int(IRinfo_table.at[accession, "IRb_CALCULATED_END"])

                IRinfo_table.at[accession, "IRa_LENGTH_COMPARED_DIFFERENCE"] = int(IRinfo_table.at[accession, "IRa_REPORTED_LENGTH"]) - int(IRinfo_table.at[accession, "IRa_CALCULATED_LENGTH"])
                IRinfo_table.at[accession, "IRb_LENGTH_COMPARED_DIFFERENCE"] = int(IRinfo_table.at[accession, "IRb_REPORTED_LENGTH"]) - int(IRinfo_table.at[accession, "IRb_CALCULATED_LENGTH"])

            else:
                log.warning("Could not calculate IRs for accession " + accession + "." + "\n".join(out.splitlines()))

                IRinfo_table.at[accession, "IRa_CALCULATED"] = "no"
                IRinfo_table.at[accession, "IRb_CALCULATED"] = "no"

                IRinfo_table.at[accession, "IRa_CALCULATED_START"] = "n.a."
                IRinfo_table.at[accession, "IRb_CALCULATED_START"] = "n.a."

                IRinfo_table.at[accession, "IRa_CALCULATED_END"] = "n.a."
                IRinfo_table.at[accession, "IRb_CALCULATED_END"] = "n.a."

                IRinfo_table.at[accession, "IRa_CALCULATED_LENGTH"] = "n.a."
                IRinfo_table.at[accession, "IRb_CALCULATED_LENGTH"] = "n.a."

                IRinfo_table.at[accession, "IRa_START_COMPARED_OFFSET"] = "n.a."
                IRinfo_table.at[accession, "IRb_START_COMPARED_OFFSET"] = "n.a."

                IRinfo_table.at[accession, "IRa_END_COMPARED_OFFSET"] = "n.a."
                IRinfo_table.at[accession, "IRb_END_COMPARED_OFFSET"] = "n.a."

                IRinfo_table.at[accession, "IRa_LENGTH_COMPARED_DIFFERENCE"] = "n.a."
                IRinfo_table.at[accession, "IRb_LENGTH_COMPARED_DIFFERENCE"] = "n.a."
                #raise Exception("Inverted repeat count not == 2")
        except Exception as err:
            log.exception("Error while calculating IRs. %s\n Skipping this accession." % (str(err)))
            continue
        # Write data to outfn
        IRinfo_table.to_csv(IRinfo_file, sep='\t', header=True)
        os.chdir(main_dir)


########
# MAIN #
########

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    # parser.add_argument("--outfn", "-o", type=str, required=True, help="path to output file that contains comparing information on reported vs. inferred IR positions and length")
    parser.add_argument("--repfn", "-r", type=str, required=True, help="path to file that contains information on reported IR positions and length")
    parser.add_argument("--data", "-d", type=str, required=True, help="Data folder containing subfolders with FASTA files")
    args = parser.parse_args()
    main(args)
