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
__info__ = 'Re-calculates the position and length of the IRs for each '\
           'plastid genome record'
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
