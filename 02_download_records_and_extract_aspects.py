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
#import xml.etree.ElementTree as ET
import os.path, subprocess, calendar
import pandas as pd
import argparse, sys
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
import ipdb
# ipdb.set_trace()

"""

#############
# FUNCTIONS #
#############

def getNewUIDs(query, outfn):
    '''
    Gets a list of all UIDs found by the search query via ESearch and EFetch.
    
    Args:
        query (str): a string that specifies the query term for NCBI
        outfn (str): name of the output file
    
    Returns:
        A set of UIDs that aren't yet included in the output file.
    '''

  # STEP 1. Get list of UIDs
    esearchargs = ["esearch", "-db", "nucleotide", "-query", query] # Diesen Aufruf evtl durch Entrezpy ersetzen? Damit bekäme ich ja sofort alle UIDs. ## MG: That's a minor issue. Don't bother with it for the moment.
    esearch = subprocess.Popen(esearchargs, stdout=subprocess.PIPE)
    efetchargs = ["efetch", "-db", "nucleotide", "-format", "uid"]
    efetch = subprocess.Popen(efetchargs, stdin=esearch.stdout, stdout=subprocess.PIPE)
    out, err = efetch.communicate()

  # STEP 2. Exclude all previously processed UIDs (i.e., those already in outfile) from input.
    uids = set(map(int,out.splitlines()))
    plastidTable = pd.read_csv(outfn, sep="\t")
    try:
        complement_set = uids - set(plastidTable['UID'])
    except Exception as e:
        log.info(("Calculation of complement set not successful:\n%s" % (e.message)))
    return complement_set


def getEntryInfo(uid):
    '''
    Gets the GenBank flatfile for given UID in XML-format, then parses XML content for relevant data and returns it as a list.
    
    Args:
        uid (str): A unique identifier (actually a long integer that only has relevance in NCBI internally)
    
    Returns:
        - accession number
        - sequence version number
        - organism name
        - sequence length
        - date the record went online
        - the first reference for this sequence including:
            -- name of authors
            -- name of the publication
            -- full citation of the publication
    '''

  # STEP 1. For a given UID, get the record summary in XML format
    esummaryargs = ["esummary", "-db", "nucleotide", "-format", "gb", "-mode", "xml", "-id", str(uid)]
    esummary = subprocess.Popen(esummaryargs, stdout=subprocess.PIPE)
    out, err = esummary.communicate()

  # STEP 2. Parse out the relevant info from XML-formatted record summary
    root = ET.fromstring(out)
    seqDaten = root.find("GBSeq")
    fields = []
    fields.append(seqDaten.find("GBSeq_primary-accession").text)
    fields.append((seqDaten.find("GBSeq_accession-version").text).split('.')[1])
    fields.append(seqDaten.find("GBSeq_organism").text)
    fields.append(seqDaten.find("GBSeq_length").text)
    # Parse and format the date that the record was first online
    month_map = {"JAN":"01", "FEB":"02", "MAR":"03", "APR":"04", "MAY":"05", "JUN":"06", "JUL":"07", "AUG":"08", "SEP":"09", "OCT":"10", "NOV":"11", "DEC":"12"}
    create_date = seqDaten.find("GBSeq_create-date").text.split('-')
    fields.append(create_date[2] + "-" + month_map[create_date[1]] + "-" + create_date[0])
    # Parse all info related to the authors and the publication
    topReference = seqDaten.find("GBSeq_references").find("GBReference")
    authors = topReference.find("GBReference_authors").findall("GBAuthor")
    authstring = ""
    for author in authors:
        authstring = authstring + author.text + ", "
    authstring = authstring[:-2]
    fields.append(authstring)
    fields.append(topReference.find("GBReference_title").text)
    fields.append(topReference.find("GBReference_journal").text)
    return fields
"""


def main(outfn, query):

  # STEP 1. Set up logger
    log = logging.getLogger(__name__)
    coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level='DEBUG', logger=log)
    ''' # Reactivate the following lines if coloredlogs turns out to be disadvantageous.
    log = logging.getLogger()
    log.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter("[%(levelname)s] - %(message)s - %(asctime)s", "%Y-%m-%d %H:%M:%S")
    ch.setFormatter(formatter)
    log.addHandler(ch)
    '''
    
"""
  # STEP 2. Check if output file already exists
    if not os.path.isfile(outfn):
        with open(outfn, "w") as summaryFile:
            summaryFile.write("UID\tACCESSION\tVERSION\tORGANISM\tSEQ_LEN\tCREATE_DATE\tAUTHORS\tTITLE\tREFERENCE\n")

  # STEP 3. Get new UIDs, parse out relevant info and save as lines in output file.
  
    # Load previously processed data from outfile
    log.info(("Obtaining previously processed data from %s" % (str(outfn))))
    plastid_summary = pd.read_csv(outfn, sep='\t', index_col=0, encoding='utf-8')
    
    # Re-initialize outfile (as new outfile)
    log.info(("Re-initializing file %s" % (str(outfn))))
    with open(outfn, "a") as summaryFile:

        # Write previously processed data to new outfile and then drop it
        log.info(("Writing previously processed data to %s" % (str(outfn))))
        plastid_summary.to_csv(summaryFile, sep='\t', header=False)
        plastid_summary.drop(plastid_summary.index, inplace=True)

        '''
        TODO:
        When a large masterlist already exists, sometimes the following exception is thrown. This exception is not trown when no masterlist is present.

        Traceback (most recent call last):
          File "01_generate_plastome_availability_table.py", line 169, in <module>
            main(args.outfn, args.query)
          File "01_generate_plastome_availability_table.py", line 140, in main
            plastid_summary.drop(plastid_summary.index, inplace=True)
          File "/home/michael_science/.local/lib/python2.7/site-packages/pandas/core/generic.py", line 1417, in drop
            indexer = -axis.isin(labels)
        TypeError: The numpy boolean negative, the `-` operator, is not supported, use the `~` operator or the logical_not function instead.
        '''

        # Iteratively write new data to new outfile
        for uid in getNewUIDs(query, outfn):
            log.info(("Reading and parsing UID %s, writing to %s" % (str(uid), str(outfn))))
            plastid_summary.loc[uid] = getEntryInfo(uid)
            #log.info(("Appending summary of UID %s to %s"  % (str(uid), outfn)))
            plastid_summary.to_csv(summaryFile, sep='\t', header=False)
            plastid_summary.drop([uid], inplace=True)
"""

########
# MAIN #
########

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    parser.add_argument("--outfn", "-o", type=str, required=True, help="path to output file")
    parser.add_argument("--query", "-q", type=str, required=False, default="Magnoliophyta[ORGN] AND 00000100000[SLEN] : 00000200000[SLEN] AND complete genome[TITLE] AND (chloroplast[TITLE] OR plastid[TITLE])", help="(Optional) Entrez query that will replace the standard query")
    args = parser.parse_args()
    main(args.outfn, args.query)

"""
