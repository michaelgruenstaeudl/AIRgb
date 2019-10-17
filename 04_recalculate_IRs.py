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
DESIGN:
    * Like in the other scripts, the evaluation of the records is conducted one by one, not all simultaneously.

    * The plastid genome sequence of each record will be bundled with the GB file of that record in a record-specific gzip file. Hence, before that record is evaluated here, that specific gzip-file must be unpacked; after the evaluation, the gzip-file must be reconstituted.

NOTES:
    * Foo bar baz
'''

#####################
# IMPORT OPERATIONS #
#####################
#import xml.etree.ElementTree as ET
import os.path, subprocess
import argparse, tarfile
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

    # STEP 2. Loop though provided archives
    archives = [os.path.abspath(x) for x in args.input]
    main_dir = os.getcwd()
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

    for archive in archives:
        # Init values that will be written to table
        accession = os.path.basename(archive).split('.')[0]
        start_a_orig = None
        start_b_orig = None
        len_a_orig = 0
        len_b_orig = 0
        len_a = 0
        len_b = 0
        start_a = None
        start_b = None

        # Extract sequence file from archive
        tar = tarfile.open(archive, "r:gz")
        tar.extractall()
        tar.close()
        acc_path = os.path.abspath(accession)
        # Change to directory containing sequence files
        os.chdir(acc_path)
        # Read in original IR start positions and lengths
        if os.path.isfile(accession + "_irpos"):
            with open(accession + "_irpos","r") as posfile:
                for line in posfile.readlines():
                    if "A:" in line:
                        start_a_orig = line.split(':')[1].split('\t')[0]
                        len_a_orig = line.split(':')[1].split('\t')[2]
                    elif "B:" in line:
                        start_b_orig = line.split(':')[1].split('\t')[0]
                        len_b_orig = line.split(':')[1].split('\t')[2]
                    else:
                        raise Exception("Unexpected line in position file.")
        # Calculate IR positions and length
        # TODO: need to execute makeblastdb before blastn can operate
        blastargs = ["blastn", "-db", str(accession) + ".fasta", "-query", accession + ".fasta", "-outfmt", "7", "-strand", "both"]
        blast_subp = subprocess.Popen(blastargs, stdout=subprocess.PIPE)
        awkargs = ["awk", "{if ($4 > 10000 && $4 < 50000) print $4, $7, $8, $9, $10}"]
        awk_subp = subprocess.Popen(awkargs, stdin=blast_subp.stdout, stdout=subprocess.PIPE)
        out, err = awk_subp.communicate()
        if len(out.splitlines()) == 2:
            ira_info = out.splitlines()[0].split()
            irb_info = out.splitlines()[1].split()
            len_a = ira_info[0]
            len_b = irb_info[0]
            start_a = ira_info[1]
            start_b = irb_info[1]
        else:
            raise Exception("Inverted repeat count not == 2")
        # Write data to outfn
        with open(outfn, "a") as outfile:
            outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (str(accession), str(start_a_orig), str(start_a), str(len_a_orig), str(len_a), str(abs(len_a - len_a_orig)), str(start_b_orig), str(start_b), str(len_b_orig), str(len_b), str(abs(len_b - len_b_orig))))
        os.chdir(main_dir)
        # Bundle and compress data
        tar = tarfile.open(acc_path + ".tar.gz", "w:gz")
        tar.add(acc_path, os.path.basename(acc_path))
        tar.close()


########
# MAIN #
########

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    parser.add_argument("--outfn", "-o", type=str, required=True, help="path to output file")
    parser.add_argument("--input", "-i", type=str, required=True, nargs='+', help="List of archive file paths containing FASTA files")
    args = parser.parse_args()
    main(args)
