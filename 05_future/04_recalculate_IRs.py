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
import os.path, subprocess
import argparse
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
__version__ = '2020.01.23.1500'

#############
# DEBUGGING #
#############
import ipdb
# ipdb.set_trace()


#############
# FUNCTIONS #
#############

def coerceToExactLocation(location):
    exactLocation = None
    if '<' in str(location) or '>' in str(location):
        exactLocation = str(location)[1:]
    else:
        exactLocation = str(location)
    return exactLocation

def main(args):
    # STEP 1. Set up logger
    log = logging.getLogger(__name__)
    coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level='DEBUG', logger=log)

    # STEP 2. Read table from Script02/Script03 and append columns that will be filled in this script
    try:
        IRinfo_file = os.path.abspath(args.repfn)
        IRinfo_table = pd.read_csv(IRinfo_file, index_col=0, sep='\t', encoding="utf-8", na_values="n.a.")
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

    # TODO: IRb_REPORTED_END contains entries that cannot be converted to numeric. Change previous scripts so only numeric data is accepted and saved
    # IRinfo_table = IRinfo_table.astype({"IRa_REPORTED": str, "IRa_REPORTED_LENGTH": pd.Int64Dtype(), "IRb_REPORTED": str, "IRb_REPORTED_LENGTH": pd.Int64Dtype(), "IRa_CALCULATED": str, "IRa_CALCULATED_START": pd.Int64Dtype(), "IRa_CALCULATED_END": pd.Int64Dtype(), "IRa_CALCULATED_LENGTH": pd.Int64Dtype(), "IRa_START_COMPARED_OFFSET": pd.Int64Dtype(), "IRa_END_COMPARED_OFFSET": pd.Int64Dtype(), "IRa_LENGTH_COMPARED_DIFFERENCE": pd.Int64Dtype(), "IRb_CALCULATED": str, "IRb_CALCULATED_START": pd.Int64Dtype(), "IRb_CALCULATED_END": pd.Int64Dtype(), "IRb_CALCULATED_LENGTH": pd.Int64Dtype(), "IRb_START_COMPARED_OFFSET": pd.Int64Dtype(), "IRb_END_COMPARED_OFFSET": pd.Int64Dtype(), "IRb_LENGTH_COMPARED_DIFFERENCE": pd.Int64Dtype(), "MUMMER_SNP_COUNT": pd.Int64Dtype(), "MUMMER_INDEL_COUNT": pd.Int64Dtype(), "MUMMER_SIMIL_SCORE": float, "CMP_DIFF_COUNT": pd.Int64Dtype(), "CONGRUENCE_MUMMER_CMP": str})
    IRinfo_table = IRinfo_table.astype({"IRa_REPORTED": str, "IRb_REPORTED": str, "IRa_CALCULATED": str, "IRb_CALCULATED": str})

    main_dir = os.getcwd()
    for folder in folders:
        # Init values that will be written to table
        accession = os.path.basename(folder)

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
            result_lines = out.splitlines()
            # Note: BLAST sometimes finds additional regions in the sequence that match the length requirements filtered for in awk. We only want the IRs, and therefore need to pick out the two regions with matching length
            if len(result_lines) > 2:
                temp_lines = []
                for i in range(len(result_lines)-1):
                    for j in range(i+1,len(result_lines)):
                        if result_lines[i].split()[1] == result_lines[j].split()[1]:
                            temp_lines.append(result_lines[i])
                            temp_lines.append(result_lines[j])
                            break
                result_lines = temp_lines
            if len(result_lines) == 2:
                # Compare the start positions of the found regions. By default, we assume IRb is located before IRa in the sequence
                if result_lines[0].split()[1] > result_lines[1].split()[1]:
                    ira_info = result_lines[1].split()
                    irb_info = result_lines[0].split()
                else:
                    ira_info = result_lines[0].split()
                    irb_info = result_lines[1].split()

                # Assign calculated info
                IRinfo_table.at[accession, "IRa_CALCULATED"] = "yes"
                IRinfo_table.at[accession, "IRb_CALCULATED"] = "yes"

                IRinfo_table.at[accession, "IRa_CALCULATED_START"] = int(ira_info[1])
                IRinfo_table.at[accession, "IRb_CALCULATED_START"] = int(irb_info[1])

                IRinfo_table.at[accession, "IRa_CALCULATED_END"] = int(ira_info[2])
                IRinfo_table.at[accession, "IRb_CALCULATED_END"] = int(irb_info[2])

                IRinfo_table.at[accession, "IRa_CALCULATED_LENGTH"] = int(ira_info[0])
                IRinfo_table.at[accession, "IRb_CALCULATED_LENGTH"] = int(irb_info[0])

                IRinfo_table.at[accession, "IRa_START_COMPARED_OFFSET"] = int(float(coerceToExactLocation(IRinfo_table.at[accession, "IRa_REPORTED_START"])) - float(coerceToExactLocation(IRinfo_table.at[accession, "IRa_CALCULATED_START"])))
                IRinfo_table.at[accession, "IRb_START_COMPARED_OFFSET"] = int(float(coerceToExactLocation(IRinfo_table.at[accession, "IRb_REPORTED_START"])) - float(coerceToExactLocation(IRinfo_table.at[accession, "IRb_CALCULATED_START"])))

                IRinfo_table.at[accession, "IRa_END_COMPARED_OFFSET"] = int(float(coerceToExactLocation(IRinfo_table.at[accession, "IRa_REPORTED_END"])) - float(coerceToExactLocation(IRinfo_table.at[accession, "IRa_CALCULATED_END"])))
                IRinfo_table.at[accession, "IRb_END_COMPARED_OFFSET"] = int(float(coerceToExactLocation(IRinfo_table.at[accession, "IRb_REPORTED_END"])) - float(coerceToExactLocation(IRinfo_table.at[accession, "IRb_CALCULATED_END"])))

                IRinfo_table.at[accession, "IRa_LENGTH_COMPARED_DIFFERENCE"] = int(float(coerceToExactLocation(IRinfo_table.at[accession, "IRa_REPORTED_LENGTH"])) - float(coerceToExactLocation(IRinfo_table.at[accession, "IRa_CALCULATED_LENGTH"])))
                IRinfo_table.at[accession, "IRb_LENGTH_COMPARED_DIFFERENCE"] = int(float(coerceToExactLocation(IRinfo_table.at[accession, "IRb_REPORTED_LENGTH"])) - float(coerceToExactLocation(IRinfo_table.at[accession, "IRb_CALCULATED_LENGTH"])))

            else:
                log.warning("Could not calculate IRs for accession " + accession + "." + "\n".join([str(line).strip() for line in result_lines]))

                IRinfo_table.at[accession, "IRa_CALCULATED"] = "no"
                IRinfo_table.at[accession, "IRb_CALCULATED"] = "no"

        except Exception as err:
            log.exception("Error while calculating IRs. %s\n Skipping this accession." % (str(err)))
            continue
        # Write data to outfn
        IRinfo_table.to_csv(IRinfo_file, sep='\t', header=True, na_rep="n.a.")
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
