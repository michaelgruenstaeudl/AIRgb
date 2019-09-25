# -*- coding: utf-8 -*-
'''
OBJECTIVE:
    This script takes a list of GenBank accession numbers and - for each accession number - downloads the record from GenBank, parses it via Biopython, extracts several aspects relevant for subsequent analyses and saves these aspects in separate files.

    The output shall be a set of files for each GenBank record:
    (a) the GenBank record in GB format
    (b) the full plastid genome sequence in FASTA format
    (c) the sequence of the IRa in FASTA format
    (d) the reverse-complemented sequence of the IRb in FASTA format

    The output is the bundled in a single gzip file.

TO DO:
    * The inverted repeats (i.e. IRa and IRb) of a plastid genome record may be labelled in different ways, depending on the record. This script shall be flexible enough to identify the different naming conventions, yet always extract only a single IR pair per record.

    * Here is starting code that you can use:
    ```
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    import os
    import tarfile

    # First, download the GenBank record currently under study.
    # Then, do the following:

    inFn = # path to, and filename of, the input record
    inDir = # where the GB record is located (so that the output files are saved right next to them)
    outFn_fullSeq = # output file containing full sequence in FASTA format
    outFn_IRa = #
    outFn_IRbRC = #

    # Make output folder for record
    outDir = os.path.join(inDir, 'output')
    if not os.path.exists(outDir):
        os.makedirs(outDir)

    # PLEASE NOTE: A lot of the code used hereafter is exeplified in the example "Random subsequences" on https://biopython.org/wiki/SeqIO

    # Read record
    rec = SeqIO.read(inFn, "genbank")
    accn = rec.id.split('.')[0]

    # Save full sequence as FASTA file
    with open(outFn_fullSeq, "w") as outHdl_fullSeq:
        SeqIO.write(rec, outHdl_fullSeq, "fasta") # But sequence must NOT be interleaved! (I know that AlignIO.write can write non-interleaved FASTA files, but I am not sure if SeqIO.write can do that. Another option would be to write the sequence to a string handle and the save the FASTA file via the regular Python write function (see "Writing to a string" in https://biopython.org/wiki/SeqIO).)

    # Identify all IRs and raise exception if not exatcly two IRs detected
    all_misc_features = [feature for feature in rec.features if feature.type=='misc_feature'] # Note: IRs are not always identified by the feature tpye "misc-feature"; there are other commonly used feature types for IRs. COnversely, there are other elements of a plastid genome that are also called "misc_feature". In short, a better identification of the IRs is necessary here!
        # Important: Raise exception if not exactly two IRs detected

    # Extract the IR sequences; reverse completement the IRb
    IRa_seq = all_misc_features[0].extract(rec).seq
    IRbRC_seq = all_misc_features[1].extract(rec).seq.reverse_complement()

    IRa_rec = SeqRecord(IRa_seq, accn+'_IRa', '', '')  # Note: Since SeqIO.write can only write record objects, I am converting the sequence object into a record object here
    IRbRC_rec = SeqRecord(IRbRC_seq, accn+'_IRb_revComp', '', '')

    # Save IRa_seq and IRb_seq_revComp as separate, non-interleaved FASTA files to outDir
    with open(outFn_IRa, "w") as outHdl_IRa, with open(outFn_IRbRC, "w") as outHdl_IRbRC:
        SeqIO.write(IRa_rec, outHdl_IRa, "fasta")
        SeqIO.write(IRbRC_rec, outHdl_IRbRC, "fasta")

    # Gzip outDir (which contains the record and all newly generated files)
    tar = tarfile.open(inDir+".tar.gz", "w:gz")
    tar.add(inDir, arcname="TarName")
    tar.close()
    ```

DESIGN:
    * Like in the other scripts, the evaluation of the records is conducted one by one, not all simultaneously.

    * The output files of each record shall be bundled together as a record-specific gzip file.

NOTES:
    * Foo bar baz
'''

#####################
# IMPORT OPERATIONS #
#####################
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import xml.etree.ElementTree as ET
import os, subprocess, argparse
import tarfile, coloredlogs, logging, time
###############
# AUTHOR INFO #
###############
__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>, '\
             'Tilman Mehl <tilmanmehl@zedat.fu-berlin.de>'
__copyright__ = 'Copyright (C) 2019 Michael Gruenstaeudl and Tilman Mehl'
__info__ = 'Compare IRs for a series of IR FASTA files'
__version__ = '2019.09.23.1700'

#############
# DEBUGGING #
#############
import ipdb
# ipdb.set_trace()

#############
# FUNCTIONS #
#############

# Saves fetched GenBank flatfile with accession number "id" to outdir and returns the path and filename of the file
def fetchGBflatfile(outdir, id):
    with open(os.path.join(outdir,str(id) + ".gb"), "w") as gbFile:
        efetchargs = ["efetch", "-db", "nucleotide", "-format", "gb", "-id", str(id)]
        efetch = subprocess.Popen(efetchargs, stdout=gbFile)
        efetch.wait()
    return os.path.join(outdir,str(id) + ".gb")

# Identifies and returns IRa and IRb as SeqRecord objects
def getInvertedRepeats(rec):
    # Hierarchy of preference for identifying inverted repeats:
    # 1. feature is of type "repeat_region" AND
    #   1.1. qualifier rpt_type == "inverted" AND
    #       1.1.1 qualifier note contains either "inverted repeat a", "inverted repeat b", "ira", or "irb"
    #       1.1.2 qualifier note is empty
    #   1.2 qualifier rpt_type is empty AND
    #       1.2.1 qualifier note contains either "inverted repeat a", "inverted repeat b", "ira", or "irb"
    #       1.2.2 qualifier note contains either "inverted repeat a", "inverted repeat b", "ira", or "irb"
    #       1.2.3 qualifier note contains either ("inverted" and "repeat") or "IR"
    # 2. feature is of type "misc_feature" AND
    #   2.1. qualifier note contains either "inverted repeat a", "inverted repeat b", "ira", or "irb"
    #   2.2. qualifier note contains either ("inverted" and "repeat") or "IR"

    IRa = None
    IRb = None
    # Checking repeat_regions
    all_repeats = [feature for feature in rec.features if feature.type=='repeat_region']
    for repeat in all_repeats:
        if "rpt_type" in repeat.qualifiers:
            if repeat.qualifiers["rpt_type"].lower() == "inverted":
                if "note" in repeat.qualifiers:
                    if "ira" in repeat.qualifiers["note"][0].lower() or "inverted repeat a" in repeat.qualifiers["note"][0].lower():
                        IRa = repeat.extract(rec)
                    elif "irb" in repeat.qualifiers["note"][0].lower() or "inverted repeat b" in repeat.qualifiers["note"][0].lower():
                        IRb = repeat.extract(rec)
                elif IRa is None:
                    IRa = repeat.extract(rec)
                elif IRb is None:
                    IRb = repeat.extract(rec)
        elif "note" in repeat.qualifiers:
            if "ira" in repeat.qualifiers["note"][0].lower() or "inverted repeat a" in repeat.qualifiers["note"][0].lower():
                IRa = repeat.extract(rec)
            elif "irb" in repeat.qualifiers["note"][0].lower() or "inverted repeat b" in repeat.qualifiers["note"][0].lower():
                IRb = repeat.extract(rec)
            elif ("inverted" in repeat.qualifiers["note"][0].lower() and "repeat" in repeat.qualifiers["note"][0].lower()) or "IR" in misc_feature.qualifiers["note"][0]:
                if IRa is None:
                    IRa = repeat.extract(rec)
                elif IRb is None:
                    IRb = repeat.extract(rec)

    # Only check misc_features if neither inverted repeat was found through the repeat_region qualifier
    # If one inverted repeat was tagged as repeat_region, the other probably would have been tagged the same
    if IRa is None and IRb is None:
        all_misc_features = [feature for feature in rec.features if feature.type=='misc_feature']
        for misc_feature in all_misc_features:
            if "ira" in misc_feature.qualifiers["note"][0].lower() or "inverted repeat a" in misc_feature.qualifiers["note"][0].lower():
                IRa = misc_feature.extract(rec)
            elif "irb" in misc_feature.qualifiers["note"][0].lower() or "inverted repeat b" in misc_feature.qualifiers["note"][0].lower():
                IRb = misc_feature.extract(rec)
            elif ("inverted" in misc_feature.qualifiers["note"][0].lower() and "repeat" in misc_feature.qualifiers["note"][0].lower()) or "IR" in misc_feature.qualifiers["note"][0]:
                if IRa is None:
                    IRa = misc_feature.extract(rec)
                elif IRb is None:
                    IRb = misc_feature.extract(rec)

    return IRa, IRb.reverse_complement()



def main(args):

  # STEP 1. Set up logger
    log = logging.getLogger(__name__)
    coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level='DEBUG', logger=log)

  # STEP 2. Read in accession numbers (if necessary)
    accNumbers = []
    if args.infile:
        with open(args.infile,"r") as infile:
            accNumbers = infile.readlines()
    else:
        accNumbers = args.list

  # STEP 3. Loop though accession numbers and save the relevant data
    for accession in accNumbers:
        # Create output folder(s) for accession data if it does not exist yet.
        outputFolder = os.path.join(args.outdir,str(accession))
        if not os.path.exists(outputFolder):
            os.makedirs(outputFolder)
        # STEP 3.1. Fetch GenBank flatfile for accession and save it
        log.info("Saving GenBank flat file for accession " + str(accession))
        gbFn = fetchGBflatfile(outputFolder, accession)
        rec = SeqIO.read(gbFn, "genbank")
        # STEP 3.2. Write sequence in FASTA format
        log.info("Writing sequence as FASTA for accession " + str(accession))
        with open(os.path.join(outputFolder, accession + ".fasta"),"w") as fastaOut:
            fastaOut.write(">" + str(rec.id) + " " + str(rec.description) + "\n")
            fastaOut.write(str(rec.seq))
        # STEP 3.3. Extract inverted repeat sequences from full sequence and write them in FASTA format
        IRa_seq = None
        IRbRC_seq = None
        IRa_seq, IRbRC_seq = getInvertedRepeats(rec)
        if not (IRa_seq is None or IRbRC_seq is None):
            log.info("Found both inverted repeats for accession " + str(accession))
            with open(os.path.join(outputFolder, accession + "_IRa.fasta"),"w") as IRa_fasta:
                IRa_fasta.write(">" + str(accession) + "_IRa\n")
                IRa_fasta.write(str(IRa_seq.seq))
            with open(os.path.join(outputFolder, accession + "_IRb_revComp.fasta"),"w") as IRb_fasta:
                IRb_fasta.write(">" + str(accession) + "_IRb_revComp\n")
                IRb_fasta.write(str(IRbRC_seq.seq))
        elif not IRa_seq is None and IRbRC_seq is None:
            log.info("Only one inverted repeat found for accession " + str(accession))
            with open(os.path.join(outputFolder, accession + "_IRa.fasta"),"w") as IRa_fasta:
                IRa_fasta.write(">" + str(accession) + "_IRa\n")
                IRa_fasta.write(str(IRa_seq.seq))
        elif IRa_seq is None and not IRbRC_seq is None:
            log.info("Only one inverted repeat found for accession " + str(accession))
            IRbRC_rec = SeqRecord(IRbRC_seq, str(accession) +'_IRb_revComp', '', '')
            with open(os.path.join(outputFolder, accession + "_IRb_revComp.fasta"),"w") as IRb_fasta:
                IRb_fasta.write(">" + str(accession) + "_IRb_revComp\n")
                IRb_fasta.write(str(IRbRC_seq.seq))
        else:
            log.info("No inverted repeats found for accession " + str(accession))
        # STEP 3.4. Bundle and compress accession data
        tar = tarfile.open(outputFolder + ".tar.gz", "w:gz")
        tar.add(outputFolder)
        tar.close()


########
# MAIN #
########

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    inputType = parser.add_mutually_exclusive_group(required=True)
    inputType.add_argument("--infile", "-i", type=str, help="File with list of NCBI nucleotide accession numbers")
    inputType.add_argument("--list", "-l", type=str, nargs='+', help="List of NCBI nucleotide accession numbers")
    parser.add_argument("--outdir", "-o", type=str, required=True, help="path to output directory")
    args = parser.parse_args()
    #rec = SeqIO.read("/home/tilmanmehl/Documents/plastid_summary/getIRtest/KX832333/KX832333.gb", "genbank")
    #print(rec.id)
    main(args)
