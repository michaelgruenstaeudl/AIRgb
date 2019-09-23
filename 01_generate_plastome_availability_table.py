# -*- coding: utf-8 -*-
'''
OBJECTIVE:
    This script generates a table of plastid genome records that are currently available on NCBI. It simultaneously ensures that no plastid genome is counted twice (issue about regular vs. RefSeq NC_ records).

    The output is a table in which each row contains the parsed information of a single record. Each row contains nine, tab-separated columns in the following order:
    1. the unique identifier,
    2a. the accession number,
    2b. synonyms of the accession number (i.e., issue about regular vs. RefSeq NC_ records),
    3. the sequence version number,
    4. the organism name,
    5. the sequence length,
    6. the date the record went online,
    7. the authors (uppermost AUTHORS line in GB-file),
    8. the name of the publication (uppermost TITLE line in GB-file), and
    9. the full citation of the publication (see uppermost JOURNAL line in GB-file)


TO DO:

    * Can you please fix the issue regarding line: "plastid_summary.drop(plastid_summary.index, inplace=True)" See my explanations there.

    * The current script version appears to duplicate the existing data (probably due to the failure in line "plastid_summary.drop(plastid_summary.index, inplace=True)"). You can see this when looking at the output file (e.g., look for duplicates of the first uid). Beware: This issue is not apparent from the log.

    * To better handle the above two isues (and similar issues), we should modify the code so that the uids are processed in a static order (as opposed to the randomness currently experienced). For example, it would be great if the uids are always processed with the olderst one starting first! (For that, the command starting with "efetchargs" should generate a date-sorted output. Is this possible?) In general, the random retrieval/processing of the uids is a weak spot in our script that should be fixed, if possible.

    * The script shall ensure that no plastid genome is counted twice (issue about regular vs. RefSeq NC_ records). If a dual counting is present, the COMMENT line of a GB-file would contain the information which other record the reference sequence is identical to.

    * Upon initialization, print out how many plastid genome entries are on GenBank (e.g., in function getNewUIDs). Then print out how many are already in the masterlist (e.g., in function main after reading in the existing masterlist).


DESIGN:

    There are thousands of plastid genome sequences on GenBank. The parsing of the records is, thus, conducted one by one, not all simultaneously. Specifically, a list of unique identifiers is first obtained and then this list is looped over.


NOTES:

    * For testing purposes (i.e., to work only on a handful of records), the start sequence length can be increased to 190000 (`00000190000[SLEN]`).

    * Searches in the sequence databases of NCBI (nucleotide, protein, EST, GSS) allow the usage of [these fields](https://www.ncbi.nlm.nih.gov/books/NBK49540/).

    * The searches automated here can be done manually in a Linux shell:
    ### Generating uidlist
    ```
    esearch -db nucleotide -query \
    "Magnoliophyta[ORGN] AND \
    00000100000[SLEN] : 00000200000[SLEN] AND \
    complete genome[TITLE] AND\
    (chloroplast[TITLE] OR plastid[TITLE]) \
    " | efetch -db nucleotide -format uid > uidlist.txt
    ```

    ### Parsing accession number for each UID
    ```
    for i in $(cat uidlist.txt); do
      ACCN=$(esummary -db nucleotide -id $i | xmllint --xpath 'string(//Caption)' -);
      echo "$ACCN"
    done;
    ```
'''

#####################
# IMPORT OPERATIONS #
#####################
import xml.etree.ElementTree as ET
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
__info__ = 'Collect summary information on all plastid sequences stored ' \
           'in NCBI GenBank'
__version__ = '2018.09.17.1900'

#############
# DEBUGGING #
#############
import ipdb
# ipdb.set_trace()

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

  # STEP 2. Check if output file already exists
    if not os.path.isfile(outfn):
        with open(outfn, "w") as summaryFile:
            summaryFile.write("UID\tACCESSION\tVERSION\tORGANISM\tSEQ_LEN\tCREATE_DATE\tAUTHORS\tTITLE\tREFERENCE\n")

  # STEP 3. Get new UIDs, parse out relevant info and save as lines in output file.

    # Load headers from outfile
    plastid_summary = pd.read_csv(outfn, nrows=0, sep='\t', index_col=0, encoding='utf-8')

    # Re-initialize outfile (as new outfile)
    #log.info(("Re-initializing file %s" % (str(outfn))))
    with open(outfn, "a") as summaryFile: # Note: This is to append, not to write anew!

        # Iteratively write new data to new outfile
        for uid in getNewUIDs(query, outfn):
            log.info(("Reading and parsing UID %s, writing to %s" % (str(uid), str(outfn))))
            plastid_summary.loc[uid] = getEntryInfo(uid)
            #log.info(("Appending summary of UID %s to %s"  % (str(uid), outfn)))
            plastid_summary.to_csv(summaryFile, sep='\t', header=False)
            plastid_summary.drop([uid], inplace=True)

########
# MAIN #
########

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    parser.add_argument("--outfn", "-o", type=str, required=True, help="path to output file")
    parser.add_argument("--query", "-q", type=str, required=False, default="Magnoliophyta[ORGN] AND 00000100000[SLEN] : 00000200000[SLEN] AND complete genome[TITLE] AND (chloroplast[TITLE] OR plastid[TITLE])", help="(Optional) Entrez query that will replace the standard query")
    args = parser.parse_args()
    main(args.outfn, args.query)
