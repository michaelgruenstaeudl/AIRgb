#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
OBJECTIVE:
    This script takes a set of IR file pairs (in FASTA format) and - for each pair - compares the IRs via application 'mummer' (function 'nucmer') and generates numerical indices for the sequence similarlity.

    The output shall be a table of plastid genome records (one record per row) that lists numerically if the inverted repeats in that record are identical.

TO DO:
    * A visual comparison of the IRs is not desired. Instead, the similarity of the IRs shall be inferred in a numerical fashion.

    * The numerical comparison shall indicate in separate columns: (a) if and how many SNPs exist between the two IRs, and (b) if length differences exist between the IRs and what length difference there is.

    * The following batch code shall be used to compare the IRs (i.e., the IR file pair) of each record:
        ```
        # Define which record to work on
        ACCN=#Define here

        # Comparison of IR FASTAs via application 'mummer', function 'nucmer'
        nucmer --maxmatch -c 100 -p $ACCN $ACCN.IRa.fas $ACCN.IRb.fas
        show-coords -r -c -l $ACCN.delta > $ACCN.coords
        show-snps $ACCN.delta > $ACCN.snps
        show-tiling $ACCN.delta > $ACCN.tiling

        # Generate side-by-side comparison
        show-aligns $ACCN.delta "$ACCN"_IRa "$ACCN"_IRb > $ACCN.alignviz

        # Align the IRs
        cat "$ACCN"_IR*.fas > tmp;
        clustalo -i tmp > "$ACCN"_clustalo.fas;

        # Deinterleaving alignment (if necessary)
        # perl -MBio::SeqIO -e 'my $seqin = Bio::SeqIO->new(-fh => \*STDIN, -format => 'fasta'); while (my $seq = $seqin->next_seq) { print ">",$seq->id,"\n",$seq->seq,"\n"; }' < "$ACCN"_clustalo.fas > "$ACCN"_clustalo_deint.fas

        # Numerical comparison via command 'cmp'
        IRA_ALN=$(sed -n '2p' "$ACCN"_clustalo_deint.fas);
        IRB_ALN=$(sed -n '4p' "$ACCN"_clustalo_deint.fas);
        cmp -bl <(echo $IRA_ALN) <(echo $IRB_ALN) | awk '{print $1,$3,$5}';
        ```

DESIGN:
    * Like in the other scripts, the evaluation of the records is conducted one by one, not all simultaneously.

    * The IR file pair of each record will be bundled with the GB file of that record in a record-specific gzip file. Hence, before that record is evaluated here, that specific gzip-file must be unpacked; after the evaluation, the gzip-file must be reconstituted.

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
__version__ = '2019.10.14.1600'

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

  # STEP 2. Loop though provided archives
    archives = [os.path.abspath(x) for x in args.input]
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

    counter = 0
    for archive in archives:
        counter += 1
        log.info("Processing accession " + counter + "/" + str(len(archives)))
        # Init values that will be written to table
        accession = os.path.basename(archive).split('.')[0]
        ir_count = 0
        len_a = None
        len_b = None
        len_diff = None
        snp_count = None
        # Extract sequence file from archive
        tar = tarfile.open(archive, "r:gz")
        tar.extractall()
        tar.close()
        acc_path = os.path.abspath(accession)
        # Change to directory containing sequence files
        os.chdir(acc_path)
        # Check if two IRs exist for this accession
        if os.path.isfile(accession + "_IRa.fasta") and os.path.isfile(accession + "_IRb_revComp.fasta"):
            log.info("Found both IRs for accession " + accession)
            ir_count = 2
            with open(accession + "_IRa.fasta", "r") as ira:
                len_a = len(ira.readlines()[1])
            with open(accession + "_IRb_revComp.fasta", "r") as irb:
                len_b = len(irb.readlines()[1])
            len_diff = abs(len_a - len_b)

            # Process IRs with mummer
            log.info("Comparing IR FASTAs using mummer's nucmer function.")
            nucmerargs = ["nucmer", "--maxmatch", "-c", "100", "-p", accession, accession + "_IRa.fasta", accession + "_IRb_revComp.fasta"]
            nucmer = subprocess.Popen(nucmerargs)
            nucmer.wait()
            coordargs = ["show-coords", "-r", "-c", "-l", accession + ".delta"]
            with open(accession + ".coords", "w") as accession_coords:
                show_coords = subprocess.Popen(coordargs, stdout=accession_coords)
                show_coords.wait()
            snpsargs = ["show-snps", accession + ".delta"]
            with open(accession + ".snps", "w") as accession_snps:
                show_snps = subprocess.Popen(snpsargs, stdout=accession_snps)
                show_snps.wait()
            tilingargs = ["show-tiling", accession + ".delta"]
            with open(accession + ".tiling", "w") as accession_tiling:
                show_tiling = subprocess.Popen(snpsargs, stdout=accession_tiling)
                show_tiling.wait()

            # Generate side-by-side comparison
            alignsargs = ["show-aligns", accession + ".delta", accession + "_IRa", accession + "_IRb_revComp"]
            with open(accession + ".alignviz","w") as accession_aligns:
                show_aligns = subprocess.Popen(alignsargs, stdout=accession_aligns)
                show_aligns.wait()

            # Align the IRs
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

            # Deinterleave created alignment FASTA file
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
            # Remove interleaved and tmp file
            os.remove("tmp")
            os.remove(accession + "_clustalo.fasta")

            # Numerical comparison via command 'cmp'
            log.info("Comparing differences between IR sequences using generated alignment.")
            with open(accession + ".compare","w") as comparefile:
                cmp_awk = subprocess.Popen("cmp -bl <(sed -n '2p' " + accession + "_clustalo_deint.fasta) <(sed -n '2p' " + accession + "_clustalo_deint.fasta) | awk '{print $1,$3,$5}'", shell=True, executable="/bin/bash", stdout=comparefile)
                cmp_awk.wait()
            with open(accession + ".compare","r") as comparefile:
                snp_count = len(comparefile.readlines())
            # Append created files to archive
            log.info("Archiving generated files.")
            tar = tarfile.open(archive, "w:gz")
            tar.add(acc_path, accession)
            tar.close()
        elif os.path.isfile(accession + "_IRa.fasta") ^ os.path.isfile(accession + "_IRb_revComp.fasta"):
            log.info("Missing one IR for accession " + accession + "! Skipping comparison.")
            ir_count = 1
        else:
            log.info("Missing both IRs for accession " + accession + "! Skipping comparison.")

        log.info("Writing gathered information to output.")
        os.chdir(main_dir)
        with open(args.outfn,"a") as outfile:
            outfile.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (str(accession), str(ir_count), str(len_a), str(len_b), str(len_diff), str(snp_count)))



########
# MAIN #
########

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    parser.add_argument("--input", "-i", type=str, required=True, nargs='+', help="List of archive file paths containing inverted repeat FASTA files")
    parser.add_argument("--outfn", "-o", type=str, required=True, help="Name of the file the results of the IR testing will be written to")
    args = parser.parse_args()
    main(args)
