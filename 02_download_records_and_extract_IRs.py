#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
OBJECTIVE:
    This script takes a list of GenBank accession numbers and - for each accession number - downloads the record from GenBank, parses it via Biopython, extracts the inverted repeat (IR) regions and saves both their reported position as well as their reported sequences to files. The IR is hereby identified either explicitly (i.e., via repeat annotations with appropriate keywords) or implicitly (i.e., via annotations of the large and small single copy region and then taking the complement set thereof).

    The output shall be a set of files for each GenBank record:
    (a) the GenBank record in GB format
    (b) the full plastid genome sequence in FASTA format
    (c) the reported positions of IRa and IRb as a table
    (d) the reported sequence of IRa and IRb (reverse-complemented), both in FASTA format

    The downloaded GB record is gzipped and saved in a folder titled "records".
    The other files generated (i.e., the three FASTA files and the table) are saved in a folder called "data" and, within that, their own subfolder (named with the accession number)

DESIGN:
    * Like in the other scripts, the processing of the individual records is conducted one by one, not all simultaneously.

    * The IRs (i.e. IRa and IRb) of a plastid genome record may be labelled in different ways, depending on the record. This script shall be flexible enough to identify the different naming conventions, yet always extract only a single IR pair per record.

    * The output files of each record shall be bundled together as a record-specific gzip file.

TO DO:

    * In function "getInvertedRepeats()", let us please add another step which - if IRa and IRb both remain None after STEP 2 - attempts to infer the position of the IR implicitly via the  large (LSC) and small single copy (SSC) regions. SPecifically, that step shall extract the positions of the LSC and the SSC and calculate the IRs as the complement set of the LSC and the SSC. Details are below in the code.

    * In STEP 3.6. function "main()", let us please add the execution of a function that explicitly checks if the IRb must be reverse complemented or not. Currently, the FASTA output of the IRb is inconsistent, as an explicit reverse complementing of the IRb appears necessary in some but not all records. Maybe this presumably inconsistent behaviour is coupled with the presence of a qualifier named "rpt_type=inverted" in "repeat_region" features (as indicated by Tilman - see lines 165 ff. herein). Details are below in the code.

NOTES:
    * none for now
'''

#####################
# IMPORT OPERATIONS #
#####################
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from fuzzywuzzy import fuzz
import pandas as pd
import os, subprocess, argparse
import tarfile, coloredlogs, logging

###############
# AUTHOR INFO #
###############
__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>, '\
             'Tilman Mehl <tilmanmehl@zedat.fu-berlin.de>'
__copyright__ = 'Copyright (C) 2019 Michael Gruenstaeudl and Tilman Mehl'
__info__ = 'Compare IRs for a series of IR FASTA files'
__version__ = '2019.11.12.1900'

#############
# DEBUGGING #
#############
import ipdb
# ipdb.set_trace()

#############
# FUNCTIONS #
#############


def fetchGBflatfile(outdir, id, log):
    ''' Saves fetched GenBank flatfile with accession number "id" to
        outdir and returns the path and filename of the file '''

    gbFile = os.path.join(outdir, str(id) + ".gb")
    with open(gbFile, "w") as outfile:
        efetchargs = ["efetch", "-db", "nucleotide", "-format", "gb", "-id", str(id)]
        efetch = subprocess.Popen(efetchargs, stdout=outfile)
        efetch.wait()
    if not os.path.isfile(gbFile):
        raise Exception("Error retrieving GenBank flatfile of accession " + str(id))
    elif os.path.getsize(gbFile) == 0:
        raise Exception("Error retrieving GenBank flatfile of accession " + str(id))

    return gbFile


def getInvertedRepeats(rec):
    ''' Identifies the IR regions from a GenBank record and returns IRa and
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

    # STEP 1. Parsing out all repeat_regions, loop through them, and attempt to identify IRs
    all_repeat_features = [feature for feature in rec.features if feature.type=='repeat_region']
    for repeat_feature in all_repeat_features:
        if "rpt_type" in repeat_feature.qualifiers:
            if repeat_feature.qualifiers["rpt_type"][0].lower() == "inverted":
                if "note" in repeat_feature.qualifiers:
                    # If the "note" qualifier contains explicit mention of which IR (a/b) we're looking at, assign it to the appropriate variable
                    if "ira" in repeat_feature.qualifiers["note"][0].lower() or "inverted repeat a" in repeat_feature.qualifiers["note"][0].lower():
                        IRa = repeat_feature
                    elif "irb" in repeat_feature.qualifiers["note"][0].lower() or "inverted repeat b" in repeat_feature.qualifiers["note"][0].lower():
                        IRb = repeat_feature
                    # If the "note" qualifier holds no information on which IR we're looking at, assign the repeat feature to one of the variables that hasn't been initialized yet.
                    elif IRa is None:
                        IRa = repeat_feature
                    elif IRb is None:
                        IRb = repeat_feature
                # If the "note" qualifier does not exist, assign the repeat feature to one of the variables that hasn't been initialized yet.
                elif IRa is None:  # If a note qualifier does not exist to inform us if an IR is IRa or IRb, then the identified IR is by default assigned to IRa first, and the second one to IRb
                    IRa = repeat_feature
                elif IRb is None:
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
            else:
                log.info("Found a repeat region without further identifying information. Ignoring this feature.")

    # STEP 2. Parsing out all misc_features, loop through them, and attempt to identify IRs
    # Only check misc_features if neither inverted repeat was found through the repeat_region qualifier; if one inverted repeat was tagged as repeat_region, the other probably would have been tagged the same
    if IRa is None and IRb is None:
        all_misc_features = [feature for feature in rec.features if feature.type=='misc_feature']
        if len(all_repeat_features) == 0 and len(all_misc_features) == 0:
            raise Exception("Record does not contain any features with which the IR are typically marked (i.e., repeat_region, misc_feature).")
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

    # STEP 3. Inferring the position of the IR implicitly by extracting the positions of the large (LSC) and small single copy (SSC) regions and calculating the IRs as the complement set thereof.

    if IRa is None and IRb is None:
        ssc = None
        lsc = None
        all_misc_features = [feature for feature in rec.features if feature.type=='misc_feature']
        for misc_feature in all_misc_features:
            if "ssc" in misc_feature.qualifiers["note"][0].lower() or "small single copy" in misc_feature.qualifiers["note"][0].lower():
                ssc = misc_feature
            if "lsc" in misc_feature.qualifiers["note"][0].lower() or "large single copy" in misc_feature.qualifiers["note"][0].lower():
                lsc = misc_feature
        if lsc and ssc:
            # TM: Do we need to add +1 to the positions?
            SSC_start = int(ssc.location.start)
            SSC_end = int(ssc.location.end)
            LSC_start = int(lsc.location.start)
            LSC_end = int(lsc.location.end)
            IRa = SeqFeature(FeatureLocation(SSC_end, LSC_start), type="misc_feature", strand=-1)
            IRa.qualifiers["note"] = "inverted repeat a"
            IRb = SeqFeature(FeatureLocation(LSC_end, SSC_start), type="misc_feature", strand=-1)
            IRb.qualifiers["note"] = "inverted repeat b"
        else:
            raise Exception("No annotations for SSC or LSC found! Cannot infer IR positions.")

    return IRa, IRb


def writeReportedIRpos(filename, IRinfo_table, accession, IRa_feature, IRb_feature):
    ''' Writes the reported start and end positions, as well as the
        lengths, of the IRs in a given SeqFeature object as a table '''

    if IRa_feature:
        IRinfo_table.loc[accession]["IRa_REPORTED_START"] = str(int(IRa_feature.location.start))
        IRinfo_table.loc[accession]["IRa_REPORTED_END"] = str(int(IRa_feature.location.end))
        IRinfo_table.loc[accession]["IRa_REPORTED_LENGTH"] = str(abs(int(IRa_feature.location.start) - int(IRa_feature.location.end)))
    else:
        IRinfo_table.loc[accession]["IRa_REPORTED_START"] = "not identified"
        IRinfo_table.loc[accession]["IRa_REPORTED_END"] = "not identified"
        IRinfo_table.loc[accession]["IRa_REPORTED_LENGTH"] = "not identified"

    if IRb_feature:
        IRinfo_table.loc[accession]["IRb_REPORTED_START"] = str(int(IRb_feature.location.start))
        IRinfo_table.loc[accession]["IRb_REPORTED_END"] = str(int(IRb_feature.location.end))
        IRinfo_table.loc[accession]["IRb_REPORTED_LENGTH"] = str(abs(int(IRb_feature.location.start) - int(IRb_feature.location.end)))
    else:
        IRinfo_table.loc[accession]["IRb_REPORTED_START"] = "not identified"
        IRinfo_table.loc[accession]["IRb_REPORTED_END"] = "not identified"
        IRinfo_table.loc[accession]["IRb_REPORTED_LENGTH"] = "not identified"

    with open(filename, "w") as outfile:
        IRinfo_table.to_csv(outfile, sep='\t', header=True)



def writeReportedIRseqs(output_folder, rec, accession, IRa_feature, IRbRC_feature):
    if not (IRa_feature is None or IRbRC_feature is None):
        with open(os.path.join(output_folder, accession + "_IRa.fasta"),"w") as IRa_fasta:
            IRa_fasta.write(">" + str(accession) + "_IRa\n")
            IRa_fasta.write(str(IRa_feature.extract(rec).seq) + "\n")
        with open(os.path.join(output_folder, accession + "_IRb_revComp.fasta"),"w") as IRb_fasta:
            IRb_fasta.write(">" + str(accession) + "_IRb_revComp\n")
            IRb_fasta.write(str(IRbRC_feature.extract(rec).seq) + "\n")
    elif not IRa_feature is None and IRbRC_feature is None:
        with open(os.path.join(output_folder, accession + "_IRa.fasta"),"w") as IRa_fasta:
            IRa_fasta.write(">" + str(accession) + "_IRa\n")
            IRa_fasta.write(str(IRa_feature.extract(rec).seq) + "\n")
    elif IRa_feature is None and not IRbRC_feature is None:
        IRbRC_rec = SeqRecord(IRbRC_feature, str(accession) +'_IRb_revComp', '', '')
        with open(os.path.join(output_folder, accession + "_IRb_revComp.fasta"),"w") as IRb_fasta:
            IRb_fasta.write(">" + str(accession) + "_IRb_revComp\n")
            IRb_fasta.write(str(IRbRC_feature.extract(rec).seq) + "\n")


def main(args):

  # STEP 1. Set up logger
    log = logging.getLogger(__name__)
    coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level='DEBUG', logger=log)

  # STEP 2. Read in accession numbers to loop over
    accNumbers = []
    if args.infn:
        with open(args.infn, "r") as infile:
            accNumbers = [line.split('\t')[1] for line in infile.read().splitlines()[1::]]
    else:
        accNumbers = args.list

    if accNumbers:
        # Create records and data folder (unless already existant)
        if not os.path.exists(args.recordsdir):
            os.makedirs(args.recordsdir)
        if not os.path.exists(args.datadir):
            os.makedirs(args.datadir)


    columns = ["ACCESSION", "IRa_REPORTED", "IRa_REPORTED_START", "IRa_REPORTED_END", "IRa_REPORTED_LENGTH", "IRb_REPORTED", "IRb_REPORTED_START", "IRb_REPORTED_END", "IRb_REPORTED_LENGTH"]
    if os.path.isfile(args.outfn):
        with open(args.outfn) as outfile:
            header = outfile.readline()
            if not header == "\t".join(columns) + "\n":
                raise Exception("Missing or malformed header in provided output file!")
        IRinfo_table = pd.read_csv(args.outfn, sep='\t', index_col=0, encoding="utf-8")
    else:
        IRinfo_table = pd.DataFrame(columns=columns[1::]) # "ACCESSION" is omitted as column name, because it will be set as name of the index column
        IRinfo_table.index.name = "ACCESSION"

  # STEP 3. Loop over accession numbers and conduct the parsing
    for accession in accNumbers:

        # STEP 3.1. Create accession subfolder
        accessionFolder = os.path.join(args.datadir, str(accession))
        if not os.path.exists(accessionFolder):
            os.makedirs(accessionFolder)
        else:
            log.warning("Folder for accession `%s` already exists. Skipping this accession." % (str(accession)))
            continue

        # STEP 3.2. Fetch GenBank flatfile for accession and save it
        log.info("Saving GenBank flat file for accession `%s`." % (str(accession)))
        try:
            gbFn = fetchGBflatfile(accessionFolder, accession, log)
        except:
            log.warning("Error retrieving accession `%s`. Skipping this accession." % (str(accession)))
            os.rmdir(accessionFolder)
            continue

        try:

            # STEP 3.3. Parse GenBank flatfile
            try:
                rec = SeqIO.read(gbFn, "genbank")
            except Exception as err:
                raise Exception("Error reading record of accession `%s`: %s. Skipping this accession." % (str(accession), str(err)))
                continue

            # STEP 3.4. Internal check
            recordID = str(rec.id).split('.')[0]
            if not recordID == str(accession): ## Note: This internal check ensures that we are actually dealing with the record that was intended to be downloaded via efetch.
                log.warning("Accession number mismatch. Expected: `%s`. Retrieved: `%s`. Skipping this accession." % (str(accession), recordID))
                continue

            # STEP 3.5. Write full sequence in FASTA format
            log.info("Writing sequence as FASTA for accession `%s`." % (str(accession)))
            with open(os.path.join(accessionFolder, str(accession)+"_completeSeq.fasta"), "w") as fastaOut:
                fastaOut.write(">" + str(accession)+"_completeSequence" + "\n")
                fastaOut.write(str(rec.seq) + "\n")


            # STEP 3.6. Extract the reported IR positions and sequences; write both to file, if present
            #           Specifically, write the reported IR positions as a table and write IR sequences in FASTA format
            IRa_feature = None
            IRb_feature = None
            IRbRC_feature = None
            if not str(accession) in IRinfo_table.index:
                IRinfo_table = IRinfo_table.append(pd.Series(name=str(accession)))
            try:
                IRa_feature, IRb_feature = getInvertedRepeats(rec)

                # TO DO: Let us add a function that explicitly checks if the IRb_feature must be reverse complemented or not. The current FASTA output is too inconsistent (i.e., sometimes the reverse-complementing appears to have been done, sometimes not). An explicit check, followed by a reverse-complementing operation (where required) appears indicated. The situation is complicated by the fact that a certain level of non-identity (i.e., up to 10% of all nucleotides) must remain permissible, even if the IRb has been correctly reverse-complemented. Hence, I suggest designing a function that takes the IRa and the IRb as received from function "getInvertedRepeats()" and then conducts an approximate string matching to see if a reverse-complementing is indicated. There are a series of Python libraries for approximate string matching available.

                # Example code (just an idea!)
                '''
                pip install fuzzywuzzy
                pip install python-Levenshtein  # optional, I think; to make it faster

                from fuzzywuzzy import fuzz

                score_noRC = fuzz.ratio(IRa_feature.extract(rec).seq, IRb_feature.extract(rec).seq)
                score_RC = fuzz.ratio(IRa_feature.extract(rec).seq, IRb_feature.extract(rec).seq.reverse_complement())

                if score_noRC > score_RC:
                    return IRb_feature.extract(rec).seq
                else:
                    return IRb_feature.extract(rec).seq.reverse_complement()
                '''

                if IRa_feature and IRb_feature:
                    score_noRC = fuzz.ratio(IRa_feature.extract(rec).seq, IRb_feature.extract(rec).seq)
                    score_RC = fuzz.ratio(IRa_feature.extract(rec).seq, IRb_feature.extract(rec).seq.reverse_complement())
                    if score_noRC > score_RC:
                        IRbRC_feature = IRb_feature
                    else:
                        IRbRC_feature = IRb_feature
                        # TM: Switching strands should be all we need. Not sure if the features we detect have their strand property set though - need to test
                        IRbRC_feature.strand = IRb_feature.strand * -1



                if not (IRa_feature is None or IRbRC_feature is None):
                    log.info("Both IRs (IRa and IRb) detected in accession `%s`." % (str(accession)))
                    IRinfo_table.loc[str(accession)]["IRa_REPORTED"] = "yes"
                    IRinfo_table.loc[str(accession)]["IRb_REPORTED"] = "yes"
                elif not IRa_feature is None and IRbRC_feature is None:
                    log.info("Only IRa detected in accession `%s`." % (str(accession)))
                    IRinfo_table.loc[str(accession)]["IRa_REPORTED"] = "yes"
                    IRinfo_table.loc[str(accession)]["IRb_REPORTED"] = "no"
                elif IRa_feature is None and not IRbRC_feature is None:
                    log.info("Only IRb detected in accession `%s`." % (str(accession)))
                    IRinfo_table.loc[str(accession)]["IRa_REPORTED"] = "no"
                    IRinfo_table.loc[str(accession)]["IRb_REPORTED"] = "yes"
                else:
                    log.info("No IRs detected in accession `%s`." % (str(accession)))
                    IRinfo_table.loc[str(accession)]["IRa_REPORTED"] = "no"
                    IRinfo_table.loc[str(accession)]["IRb_REPORTED"] = "no"
            except Exception as err:
                IRinfo_table.loc[str(accession)]["IRa_REPORTED"] = "no"
                IRinfo_table.loc[str(accession)]["IRb_REPORTED"] = "no"
                writeReportedIRpos(args.outfn, IRinfo_table, str(accession), IRa_feature, IRbRC_feature)
                raise Exception("Error while extracting IRs for accession `%s`: %s Skipping further processing of this accession." % (str(accession), str(err)))
                continue

            # STEP 3.7. Write the reported IR positions as a table
            writeReportedIRpos(args.outfn, IRinfo_table, str(accession), IRa_feature, IRbRC_feature)

            # STEP 3.8. Write IR sequences in FASTA format
            writeReportedIRseqs(accessionFolder, rec, accession, IRa_feature, IRbRC_feature)

        except Exception as err:
            log.warning(str(err))

        # STEP 3.9. Compress the GenBank flatfile and delete uncompressed flatfile
        finally:
            tar = tarfile.open(os.path.join(args.recordsdir, accession + ".tar.gz"), "w:gz")
            tar.add(gbFn, os.path.basename(gbFn))
            tar.close()
            os.remove(gbFn)


########
# MAIN #
########

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    EitherOr = parser.add_mutually_exclusive_group(required=True)
    EitherOr.add_argument("--infn", "-i", type=str, help="path to input file; input is a summary table of NCBI accessions (tab-delimited, accession numbers in second column)")
    EitherOr.add_argument("--inlist", "-l", type=str, nargs='+', help="List of NCBI nucleotide accession numbers")
    parser.add_argument("--outfn", "-o", type=str, required=True, help="path to output file that contains information on IR positions and length")
    parser.add_argument("--recordsdir", "-r", type=str, required=False, default="./records/", help="path to records directory")
    parser.add_argument("--datadir", "-d", type=str, required=False, default="./data/", help="path to data directory")
    args = parser.parse_args()
    main(args)
