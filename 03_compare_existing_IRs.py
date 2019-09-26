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
import os.path, subprocess
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
__version__ = '2018.09.18.1200'

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
    for archive in args.input:
        accession = os.path.basename(archive).split('.')[0]
        # Extract sequence file from archive
        tar = tarfile.open(archive, "r:gz")
        tar.extractall()
        tar.close()
        # Change to directory containing sequence files
        os.chdir(accession)
        # Process IRs with mummer
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




########
# MAIN #
########

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    parser.add_argument("--input", "-i", type=str, required=True, nargs=+, help="List of archive file names containing inverted repeat FASTA files")
    parser.add_argument("-outfn", "-o", type=str, required=True, help="Name of the file the results of the IR testing will be written to")
    args = parser.parse_args()
    main(args)
