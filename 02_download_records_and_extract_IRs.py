#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
OBJECTIVE:
    This script takes a list of GenBank accession numbers and - for each accession number - downloads the record from GenBank, parses it via Biopython, extracts the inverted repeat (IR) regions and saves both their reported position as well as their reported sequences to files.

    The output shall be a set of files for each GenBank record:
    (a) the GenBank record in GB format
    (b) the full plastid genome sequence in FASTA format
    (c) the reported positions of IRa and IRb as a table
    (d) the reported sequence of IRa and IRb (reverse-complemented), both in FASTA format

    The downloaded GB record is gzipped and saved in a folder titled "records".
    The other files generated (i.e., the three FASTA files and the table) are saved in a folder called "data" and, within that, their own subfolder (named with the accession number)

TO DO:
    # Extract the IR sequences; reverse completement the IRb
    IRa_seq = all_misc_features[0].extract(rec).seq
    IRbRC_seq = all_misc_features[1].extract(rec).seq.reverse_complement()
    # Note: The code regarding reverse complementing (or not) is correct, if str(IRa) == str(IRb) for a typical record. (It's the exceptions to this rule that we are curious about in this research project.)

DESIGN:
    * Like in the other scripts, the processing of the individual records is conducted one by one, not all simultaneously.
    
    * The IRs (i.e. IRa and IRb) of a plastid genome record may be labelled in different ways, depending on the record. This script shall be flexible enough to identify the different naming conventions, yet always extract only a single IR pair per record.

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
__version__ = '2019.10.30.1130'

#############
# DEBUGGING #
#############
import ipdb
# ipdb.set_trace()

#############
# FUNCTIONS #
#############


def fetchGBflatfile(outdir, id):
    ''' Saves fetched GenBank flatfile with accession number "id" to 
        outdir and returns the path and filename of the file '''
        
    ## TO DO: Please implement raising an exception if a record is not successfully downloaded; then please catch the exception (as well as do a warning message via the logger) in main(), so that looping over the different records is not interrupted but just results in the warning message.
        
    with open(os.path.join(outdir,str(id) + ".gb"), "w") as gbFile:
        efetchargs = ["efetch", "-db", "nucleotide", "-format", "gb", "-id", str(id)]
        efetch = subprocess.Popen(efetchargs, stdout=gbFile)
        efetch.wait()
    return os.path.join(outdir, str(id)+".gb")


def getInvertedRepeats(rec):
    ''' Identifies the IR regions from a GB record and returns IRa and 
        IRb as SeqFeature objects '''
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
    
    # Parsing out all repeat_regions of record
    all_repeat_features = [feature for feature in rec.features if feature.type=='repeat_region']
    ## TO DO: Produce info via logger and raise exception if "all_repeat_features" were empty (very rare but theoretically possible!)

    # Looping through repeat regions and attempting to identify IRs
    for repeat_feature in all_repeat_features:
        if "rpt_type" in repeat_feature.qualifiers:
            if repeat_feature.qualifiers["rpt_type"][0].lower() == "inverted":
                if "note" in repeat_feature.qualifiers:
                    if "ira" in repeat_feature.qualifiers["note"][0].lower() or "inverted repeat a" in repeat_feature.qualifiers["note"][0].lower():
                        IRa = repeat_feature
                    elif "irb" in repeat_feature.qualifiers["note"][0].lower() or "inverted repeat b" in repeat_feature.qualifiers["note"][0].lower():
                        IRb = repeat_feature
                elif IRa is None:  ## TO DO: Does this line represent "if qualifier note is empty"?
                    IRa = repeat_feature
                elif IRb is None:  ## TO DO: Does this line represent "if qualifier note is empty"?
                    IRb = repeat_feature
        elif "note" in repeat_feature.qualifiers:
            if "ira" in repeat_feature.qualifiers["note"][0].lower() or "inverted repeat a" in repeat_feature.qualifiers["note"][0].lower():
                IRa = repeat_feature
            elif "irb" in repeat_feature.qualifiers["note"][0].lower() or "inverted repeat b" in repeat_feature.qualifiers["note"][0].lower():
                IRb = repeat_feature
            elif ("inverted" in repeat_feature.qualifiers["note"][0].lower() and "repeat" in repeat_feature.qualifiers["note"][0].lower()) or "IR" in misc_feature.qualifiers["note"][0]:
                if IRa is None:
                    IRa = repeat_feature
                elif IRb is None:
                    IRb = repeat_feature

    # Only check misc_features if neither inverted repeat was found through the repeat_region qualifier
    # If one inverted repeat was tagged as repeat_region, the other probably would have been tagged the same
    if IRa is None and IRb is None:
        all_misc_features = [feature for feature in rec.features if feature.type=='misc_feature']
        for misc_feature in all_misc_features:
            if "ira" in misc_feature.qualifiers["note"][0].lower() or "inverted repeat a" in misc_feature.qualifiers["note"][0].lower():
                IRa = misc_feature
            elif "irb" in misc_feature.qualifiers["note"][0].lower() or "inverted repeat b" in misc_feature.qualifiers["note"][0].lower():
                IRb = misc_feature
            elif ("inverted" in misc_feature.qualifiers["note"][0].lower() and "repeat" in misc_feature.qualifiers["note"][0].lower()) or "IR" in misc_feature.qualifiers["note"][0]:
                if IRa is None:
                    IRa = misc_feature
                elif IRb is None:
                    IRb = misc_feature

    # Biopython automatically seems to extract IRb in a reverse complement fashion
    # This was true for a record (NC_043815) that had two "repeat_region" features with "rpt_type=inverted", one of which had its sequence marked "complement"
    # Will have to test if the behaviour changes for misc_feature records or records where rpt_type=inverted is omitted
    
    # MG: The code regarding reverse complementing (or not) is correct, if str(IRa) == str(IRb) for a typical record. (It's the exceptions to this rule that we are curious about in this research project.)
    return IRa, IRb


def writeReportedIRpos(filename, IRa_feature, IRb_feature):
    ''' Writes the reported start and end positions, as well as the 
        lengths, of the IRs in a given SeqFeature object as a table '''

    with open(filename, "w") as outfile:
        outfile.write("IR\tStart\tEnd\tLength\n")
        if IRa_feature:
            outfile.write("Reported IRa:" + "\t" + str(int(IRa_feature.location.start)) + "\t" + str(int(IRa_feature.location.end)) + "\t" + str(abs(int(IRa_feature.location.start) - int(IRa_feature.location.end))) + "\n")
            
        ## TO DO: if not IRa_feature: write "not identified" for each value

        if IRb_feature:
            outfile.write("Reported IRb:" + "\t" + str(int(IRb_feature.location.start)) + "\t" + str(int(IRb_feature.location.end)) + "\t" + str(abs(int(IRb_feature.location.start) - int(IRb_feature.location.end))))

        ## TO DO: if not IRb_feature: write "not identified" for each value


def main(args):

  # STEP 1. Set up logger
    log = logging.getLogger(__name__)
    coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level='DEBUG', logger=log)

  # STEP 2. Read in accession numbers to loop over
    accNumbers = []
    if args.infile:
        with open(args.infile, "r") as infile:
            
            ## TO DO: Parsing needs to be improved. Only the second column of the infile (disregarding the title line) contains the accession numbers.
            accNumbers = infile.read().splitlines()
    else:
        accNumbers = args.list

    if accNumbers:
        # Create records and data folder (unless already existant)
        if not os.path.exists(args.recordsdir):
            os.makedirs(args.recordsdir)
        if not os.path.exists(args.datadir):
            os.makedirs(args.datadir)

  # STEP 3. Loop over accession numbers and conduct the parsing
    for accession in accNumbers:
        # Create accession subfolder
        
        ## TO DO: The accession folder must not yet exist at this point; if it does, raise a warning and skip that particular accession
        accessionFolder = os.path.join(args.recordsdir, str(accession))
        if not os.path.exists(accessionFolder):
            os.makedirs(accessionFolder)
            
        ## TO DO: Check if recordID is identical to str(accession); it should be; if not, raise a warning and skip that particular accession
        recordID = str(rec.id).split('.')[0]


        # STEP 3.1. Fetch GenBank flatfile for accession and save it
        log.info("Saving GenBank flat file for accession `%s`." % (str(accession)))

        ## TO DO: Implement catching an exception from fetchGBflatfile if a record is not successfully downloaded; in such cases, the accession is simply skipped, with a warning message logged
        gbFn = fetchGBflatfile(accessionFolder, accession)

        ## TO DO: Implement catching an exception from SeqIO.read if a record is not successfully parsed (e.g., because it is malformed); in such cases, the accession is simply skipped, with a warning message logged
        rec = SeqIO.read(gbFn, "genbank")


        # STEP 3.2. Write full sequence in FASTA format
        log.info("Writing sequence as FASTA for accession `%s`." % (str(accession)))
        with open(os.path.join(accessionFolder, str(accession)+"_completeSeq.fasta"), "w") as fastaOut:
            fastaOut.write(">" + str(accession)+"_completeSequence" + "\n")
            fastaOut.write(str(rec.seq) + "\n")


        # STEP 3.3. Extract the reported IR positions and sequences; write both to file, if present
        #           Specifically, write the reported IR positions as a table and write IR sequences in FASTA format
        IRa_feature = None
        IRbRC_feature = None


        # Step 3.3.1. Parse reported IR positions and check if both IRa and IRb are present

        ## TO DO: Implement catching an exception from getInvertedRepeats so that the loop is not interrupted; in case of an exception, the accession is simply skipped, with a warning message logged
        IRa_feature, IRbRC_feature = getInvertedRepeats(rec)
        if not (IRa_feature is None or IRbRC_feature is None):
            log.info("Both IRs (IRa and IRb) detected in accession `%s`." % (str(accession)))
        elif not IRa_feature is None and IRbRC_feature is None:
            log.info("Only IRa detected in accession `%s`." % (str(accession)))
        elif IRa_feature is None and not IRbRC_feature is None:
            log.info("Only IRb detected in accession `%s`." % (str(accession)))
        else:
            log.info("No IRs detected in accession `%s`." % (str(accession)))


        # Step 3.3.2. Write the reported IR positions as a table
        tableFn = os.path.join(accessionFolder, str(accession)+"_ReportedIRpositions.tsv")
        writeReportedIRpos(tableFn, IRa_feature, IRbRC_feature)


        # Step 3.3.3. Write write IR sequences in FASTA format
        
        ## TO DO: Move the following lines into their own function, which should be preferentially named "writeReportedIRseqs"
        if not (IRa_feature is None or IRbRC_feature is None):
            with open(os.path.join(accessionFolder, accession + "_IRa.fasta"),"w") as IRa_fasta:
                IRa_fasta.write(">" + str(accession) + "_IRa\n")
                IRa_fasta.write(str(IRa_feature.extract(rec).seq) + "\n")
            with open(os.path.join(accessionFolder, accession + "_IRb_revComp.fasta"),"w") as IRb_fasta:
                IRb_fasta.write(">" + str(accession) + "_IRb_revComp\n")
                IRb_fasta.write(str(IRbRC_feature.extract(rec).seq) + "\n")
        elif not IRa_feature is None and IRbRC_feature is None:
            with open(os.path.join(accessionFolder, accession + "_IRa.fasta"),"w") as IRa_fasta:
                IRa_fasta.write(">" + str(accession) + "_IRa\n")
                IRa_fasta.write(str(IRa_feature.extract(rec).seq) + "\n")
        elif IRa_feature is None and not IRbRC_feature is None:
            IRbRC_rec = SeqRecord(IRbRC_feature, str(accession) +'_IRb_revComp', '', '')
            with open(os.path.join(accessionFolder, accession + "_IRb_revComp.fasta"),"w") as IRb_fasta:
                IRb_fasta.write(">" + str(accession) + "_IRb_revComp\n")
                IRb_fasta.write(str(IRbRC_feature.extract(rec).seq) + "\n")


        ## TO DO: Only the downloaded GB record is gzipped and saved in a folder titled "records". The other files generated (i.e., the three FASTA files and the table) are saved in a folder called "data" and, within that, their own subfolder (named with the accession number)
        ## Desired output per record: 
        # records/NC_027250.gz
        # data/NC_027250/NC_027250_IRa.fasta
        # data/NC_027250/NC_027250_IRb_revComp.fasta
        # data/NC_027250/NC_027250_completeSeq.fasta
        # data/NC_027250/NC_027250_ReportedIRpositions.tsv

        # STEP 3.4. Bundle and compress accession data
        tar = tarfile.open(accessionFolder + ".tar.gz", "w:gz")
        tar.add(accessionFolder, os.path.basename(accessionFolder))
        tar.close()


########
# MAIN #
########

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    inputType = parser.add_mutually_exclusive_group(required=True)
    inputType.add_argument("--infile", "-i", type=str, help="File with list of NCBI nucleotide accession numbers")
    inputType.add_argument("--list", "-l", type=str, nargs='+', help="List of NCBI nucleotide accession numbers")
    parser.add_argument("--recordsdir", "-r", type=str, required=False, default="./records/", help="path to records directory")
    parser.add_argument("--datadir", "-a", type=str, required=False, default="./data/", help="path to data directory")
    args = parser.parse_args()
    main(args)
