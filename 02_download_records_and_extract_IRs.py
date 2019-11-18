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
    * The IRs (i.e. IRa and IRb) of a plastid genome record may be annotated in different ways (i.e., via different features and feature qualifiers), depending on the record. This script is flexible enough to identify the different naming conventions, yet always extract only a single IR pair per record. 
    
    * This script can also infer the position of the IR implicitly via the position of large (LSC) and small single copy (SSC) regions. Specifically, it extracts the positions of the LSC and the SSC and calculate the IRs as the complement set of the LSC and the SSC, if no explicit annotation information for the IR is present.

    * The plastid genome records on GenBank are inconsistent in their annotations of the reading direction of the IR. Specifically, the annotations sometimes indicate that one of the IRs is on the opposite strand than its counterpart (i.e., is reverse-complemented compared to its counterpart), sometimes not. In order to avoid issues when extracting the IRs in FASTA-format, an explicit check of the IRs per genome, followed by a reverse-complementing operation (where required), is indicated. The situation is complicated by the fact that a certain level of non-identity (i.e., up to 10% of all nucleotides) must remain permissible between the IRs, because the identification of such a non-identity is the overall objective of this investigation. Consequently, this script explicitly checks for each GenBank record (in function "getInvertedRepeats()") if one of the IR features of a record must be reverse complemented or not. Specifically, the string conducts approximate string comparisons and, based on the results, switches the strand information for one of the IRs where necessary.

    * The output files of each record shall be bundled together as a record-specific gzip file.

TO DO:
    * see code locations indicated with "# TO DO"

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
__version__ = '2019.11.18.1900'

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
        
        # TO DO: The inference of the IR via feature "repeat_region" with the qualifier "rpt_type" and the qualifier value "inverted" shall only be permissible if the length of the feature is at least 1000bp. (Otherwise, small sequence curiosities that were annotated as 'inverted repeats' in very early plastid genomes would be counted as genuine IRs.) Please have this script write a log.info-line in such a case that an 'inverted repeat' less than 1000bp was detected but not counted (due to its short length).
        
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
                print("Found a repeat region without further identifying information. Ignoring this feature.")
                #log.info("Found a repeat region without further identifying information. Ignoring this feature.")

    # STEP 2. Parsing out all misc_features, loop through them, and attempt to identify IRs
    # Only check misc_features if neither inverted repeat was found through the repeat_region qualifier; if one inverted repeat was tagged as repeat_region, the other probably would have been tagged the same
    
    if IRa is None and IRb is None:
        all_misc_features = [feature for feature in rec.features if (feature.type=='misc_feature' and 'pseudo' not in feature.qualifiers)]
        all_qualifiers = [misc_feature.qualifiers for misc_feature in all_misc_features]
        if len(all_repeat_features) == 0 and len(all_misc_features) == 0:
            raise Exception("Record does not contain any features which the IR are typically marked with (i.e., feature `repeat_region`, `misc_feature`).")
        if "note" not in [q.keys() for q in all_qualifiers][0]:
            raise Exception("Record does not contain any qualifiers for feature `misc_feature` which the IR are typically named with (i.e., qualifier `note`).")
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


    ## TO DO: In very early plastid genome records, the authors only listed the junction sites between the single-copy and the IR regions in a GenBank record. See, for example, record NC_001320 which provides three of the four junction sites, namely LSC-IRb, IRb-SSC, and SSC-IRa. The fourth junction site (IRa-LSC) is always the last nucleotide in the sequence and must, thus, not be given explicitly:
    '''
        misc_feature    80592..80593
                     /note="[JLB]; Junction LSC-IR"
        misc_feature    101391..101392
                     /note="[JSB]; Junction IR-SSC"
        misc_feature    113726..113727
                     /note="[JSA]; Junction SSC-IR"
    '''
    # Here's another example, this time from record NC_002202:
    '''
        misc_feature    82719..82720
                     /note="JLB, junction LSC-IRB"
        misc_feature    107792..107793
                     /note="JSB, junction IRB-SSC"
        misc_feature    125652..125653
                     /note="JSA, junction SSC-IRA"
    '''
    # Here's yet another example that illustrates the diversity of how these junction sites were listed:
    '''
        misc_feature    1..159886
                     /note="IR-LSC boundary"
        misc_feature    88152..88153
                     /note="LSC-IR boundary"
        misc_feature    114536..114537
                     /note="IR-SSC boundary"
        misc_feature    133500..133501
                     /note="SSC-IR boundary"
    '''
    # In such cases, it is necessary to have a separate block of code to determine the correct IR position. Please design and add such code below! The current code already identifies the IRs, but mistakes the junction sites as the full length and, thus, lists the IRs as having a length of only one or a few bp.
    
    # STEP 3. Inferring the position of the IR if only the junction sites are provided
    if (IRa is not None and len(IRa.extract(rec).seq)<100) or (IRb is not None and len(IRb.extract(rec).seq)<100):  # IMPORTANT: Please keep "or" in this line, as the circumstance that the IRa-LSC junction site is between the last and the first bp of the plastome gives the IRa sometimes the appearance that it has the length of the entrire plastome (and is, thus >100bp).
        raise Exception("IRa (length: `%s`) or IRb (length: `%s`) are far too short to be genuine IRs" % (len(IRa.extract(rec).seq), len(IRb.extract(rec).seq))) # remove once better code is implemented


    # STEP 4. Inferring the position of the IR implicitly by extracting the positions of the large (LSC) and small single copy (SSC) regions and calculating the IRs as the complement set thereof.
    if IRa is None and IRb is None:
        ssc = None
        lsc = None
        all_misc_features = [feature for feature in rec.features if (feature.type=='misc_feature' and 'pseudo' not in feature.qualifiers)]
        all_qualifiers = [misc_feature.qualifiers for misc_feature in all_misc_features]
        if len(all_misc_features) == 0:
            raise Exception("Record does not contain any features which the single-copy regions are typically marked with (i.e., feature `misc_feature`).")
        if "note" not in [q.keys() for q in all_qualifiers][0]:
            raise Exception("Record does not contain any qualifiers for feature `misc_feature` which the single-copy regions are typically named with (i.e., qualifier `note`).")
        for misc_feature in all_misc_features:
            if "ssc" in misc_feature.qualifiers["note"][0].lower() or "small single copy" in misc_feature.qualifiers["note"][0].lower():
                ssc = misc_feature
            if "lsc" in misc_feature.qualifiers["note"][0].lower() or "large single copy" in misc_feature.qualifiers["note"][0].lower():
                lsc = misc_feature
        if lsc and ssc:
            SSC_start = int(ssc.location.start)
            SSC_end = int(ssc.location.end)
            LSC_start = int(lsc.location.start)
            LSC_end = int(lsc.location.end)
            IRa = SeqFeature(FeatureLocation(SSC_end+1, LSC_start-1), type="misc_feature", strand=-1)
            IRa.qualifiers["note"] = "inverted repeat a"
            IRb = SeqFeature(FeatureLocation(LSC_end+1, SSC_start-1), type="misc_feature", strand=-1)
            IRb.qualifiers["note"] = "inverted repeat b"
        else:
            raise Exception("Record does not contain the information necessary to infer the position of either the IR or the single-copy regions.")

    return IRa, IRb


def writeReportedIRpos(filename, IRinfo_table, accession, IRa_feature, IRb_feature):
    ''' Writes the reported start and end positions, as well as the
        lengths, of the IRs in a given SeqFeature object as a table '''

    if IRa_feature:
        IRinfo_table.loc[accession]["IRa_REPORTED_START"] = str(int(IRa_feature.location.start))
        IRinfo_table.loc[accession]["IRa_REPORTED_END"] = str(int(IRa_feature.location.end))
        IRinfo_table.loc[accession]["IRa_REPORTED_LENGTH"] = str(abs(int(IRa_feature.location.start) - int(IRa_feature.location.end)))
    else:
        IRinfo_table.loc[accession]["IRa_REPORTED_START"] = "n.a."
        IRinfo_table.loc[accession]["IRa_REPORTED_END"] = "n.a."
        IRinfo_table.loc[accession]["IRa_REPORTED_LENGTH"] = "n.a."

    if IRb_feature:
        IRinfo_table.loc[accession]["IRb_REPORTED_START"] = str(int(IRb_feature.location.start))
        IRinfo_table.loc[accession]["IRb_REPORTED_END"] = str(int(IRb_feature.location.end))
        IRinfo_table.loc[accession]["IRb_REPORTED_LENGTH"] = str(abs(int(IRb_feature.location.start) - int(IRb_feature.location.end)))
    else:
        IRinfo_table.loc[accession]["IRb_REPORTED_START"] = "n.a."
        IRinfo_table.loc[accession]["IRb_REPORTED_END"] = "n.a."
        IRinfo_table.loc[accession]["IRb_REPORTED_LENGTH"] = "n.a."

    with open(filename, "w") as outfile:
        IRinfo_table.to_csv(outfile, sep='\t', header=True)



def writeReportedIRseqs(output_folder, rec, accession, IRa_feature, IRbRC_feature):
    ''' Writes the reported IRs to file in FASTA format '''

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
        accNumbers = args.inlist

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
                IRa_feature, IRbRC_feature = getInvertedRepeats(rec)
                if IRa_feature and IRbRC_feature:
                    score_noRC = fuzz.ratio(IRa_feature.extract(rec).seq, IRbRC_feature.extract(rec).seq)
                    score_RC = fuzz.ratio(IRa_feature.extract(rec).seq, IRbRC_feature.extract(rec).seq.reverse_complement())
                    if score_noRC < score_RC:
                        # Switching strands should be sufficient
                        IRbRC_feature.strand = IRbRC_feature.strand * -1

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
