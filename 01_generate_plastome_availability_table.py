#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
OBJECTIVE:
    This script generates a table of plastid genome records that are currently available on NCBI. It simultaneously ensures that no plastid genome is counted twice (issue about regular vs. RefSeq NC_ records).

    The output is a table in which each row contains the parsed information of a single record. Each row contains nine, tab-separated columns in the following order:
    01. the unique identifier,
    02a. the accession number,
    02b. synonyms of the accession number (i.e., issue about regular vs. RefSeq NC_ records),
    03. the sequence version number,
    04. the organism name,
    05. the sequence length,
    06. the date the record went online,
    07. the authors (uppermost AUTHORS line in GB-file),
    08. the name of the publication (uppermost TITLE line in GB-file), and
    09. the full citation of the publication (see uppermost JOURNAL line in GB-file)

    ! Still to do: !
    10. Any note if a REFSEQ accession number exists and to which regular accession number it is equal to.


TO DO:

    * The script shall ensure that no plastid genome is counted twice (issue about regular vs. RefSeq NC_ records). If a dual counting is present, the COMMENT line of a GB-file would contain the information which other record the reference sequence is identical to.

    * Once the script works well, let us parse the date of the oldest existing UID (from the already existing output file) and limit the esearch to searching NCBI for records only after that date (use esearch option "-mindate"; see "esearch -h" for more info).

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
__version__ = '2019.10.09.1330'

#############
# DEBUGGING #
#############
import ipdb
# ipdb.set_trace()

#############
# FUNCTIONS #
#############

def getQueriedUIDs(query, outfn):
    '''
    Gets a list of all UIDs found by the search query via ESearch and EFetch.

    Args:
        query (str): a string that specifies the query term for NCBI

    Returns:
        A list of UIDs as integers that is returned by the query.
    '''

  # STEP 1. Get list of UIDs
    esearchargs = ['esearch', '-db', 'nucleotide', '-sort', '"Date Released"', '-query', query]
    esearch = subprocess.Popen(esearchargs, stdout=subprocess.PIPE)

    ### TO DO: Due to the date sorting in the esearchargs, the UIDs are sorted with the newest first. However, we want the oldest first. This could be done after esearch or after efetch, but may have different effects!

    efetchargs = ["efetch", "-db", "nucleotide", "-format", "uid"]
    efetch = subprocess.Popen(efetchargs, stdin=esearch.stdout, stdout=subprocess.PIPE)
    out, err = efetch.communicate()

    return map(int, out.splitlines().reverse())

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
    accession = seqDaten.find("GBSeq_primary-accession").text
    fields.append(accession)
    fields.append((seqDaten.find("GBSeq_accession-version").text).split('.')[1])
    fields.append(seqDaten.find("GBSeq_organism").text)
    fields.append(seqDaten.find("GBSeq_length").text)

    # Parse and format the date that the record was first online
    month_map = {"JAN":"01", "FEB":"02", "MAR":"03", "APR":"04", "MAY":"05", "JUN":"06", "JUL":"07", "AUG":"08", "SEP":"09", "OCT":"10", "NOV":"11", "DEC":"12"}
    create_date = seqDaten.find("GBSeq_create-date").text.split('-')
    fields.append(create_date[2] + "-" + month_map[create_date[1]] + "-" + create_date[0])

    # Parse all info related to the authors and the publication
    references = seqDaten.find("GBSeq_references").findall("GBReference")
    authstring = ""
    title = ""
    citation = ""

    # Properly format the author output
    for ref in references:
        # Look for a reference that has authors (not all entries have a reference with authors)
        authors = ref.find("GBReference_authors")
        if authors:
            title = ref.find("GBReference_title").text
            citation = ref.find("GBReference_journal").text
            authors = ref.find("GBReference_authors").findall("GBAuthor")
            for author in authors:
                authstring = authstring + author.text.replace(","," ") + ", " #Not sure if this does what it is supposed to do? -> Testing
                #authstring = authstring[:-2]
            break
    fields.append(authstring)
    fields.append(title)
    fields.append(citation)

    # Parse comment field for RefSeq
    note = ""
    origseq = None
    if accession[:3] = "NC_":
        comments = root.find("GBSeq_comment").text
        for comment in comments.split(";"):
            if "PROVISIONAL REFSEQ: This record has not yet been subject to final NCBI review" in comment:
                origseq = comment.split(" ")[-1][:-1]
                note = "The reference sequence %s is identical to record %s." % (accession, origseq)

    return fields, origseq

def constructBlacklist(df):
    '''
        returns a set of sequence accessions for which a RefSeq exists

        Args:
            df (pandas data frame): a table containing plastid summaries
    '''
    blacklist = set()
    for note in df['NOTE']:
        if "The reference sequence NC_" in note:
            blacklist.add(note.split(" ")[-1][:-1])
    return blacklist

def main(outfn, query):

  # STEP 1. Set up logger
    log = logging.getLogger(__name__)
    coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level='DEBUG', logger=log)

  # STEP 2. Check if output file already exists, read existing UIDs
    UIDs_alreadyProcessed = set()
    if os.path.isfile(outfn):
        with open(outfn, "r") as outputFile:
            UIDs_alreadyProcessed = set(map(int, [row.split('\t')[0] for row in outputFile]))
            log.info(("Summary file `%s` already exists. %s UIDs read." % (str(outfn), str(len(UIDs_alreadyProcessed)))))
    else:
        with open(outfn, "w") as outputFile:
            outputFile.write("UID\tACCESSION\tVERSION\tORGANISM\tSEQ_LEN\tCREATE_DATE\tAUTHORS\tTITLE\tREFERENCE\nNOTE\n")
            log.info(("Summary file `%s` does not exist; generating new file. Thus, no UIDs read." % (str(outfn))))

  # STEP 3. Get all existing UIDs and calculate which to be processed
    UIDs_allExisting = getQueriedUIDs(query)
    log.info(("Number of unique UIDs currently on NCBI: `%s`" % (str(len(UIDs_allExisting)))))

    UIDs_notYetProcessed = set()
    if len(UIDs_alreadyProcessed) > 0:
        plastidTable = pd.read_csv(outfn, sep="\t")
        try:
            UIDs_notYetProcessed = UIDs_allExisting - UIDs_alreadyProcessed
        except Exception as e:
            log.info(("Calculation of complement set not successful:\n%s" % (e.message)))

    log.info(("Number of UIDs to be processed: `%s`" % (str(len(UIDs_notYetProcessed)))))

  # STEP 4. Obtain info from new UIDs and save as lines in output file.
    if len(UIDs_notYetProcessed) > 0:
      # Load format of the summary file
        outputHandle = pd.read_csv(outfn, nrows=0, sep='\t', index_col=0, encoding='utf-8')
      # Append to outfile
        blacklist = set()
        with open(outfn, "a") as outputFile:
            for uid in UIDs_notYetProcessed:
                log.info(("Reading and parsing UID %s, writing to %s" % (str(uid), str(outfn))))

                entry, origseq = getEntryInfo(uid) # Note: This is where the heavy-lifting is done!
                if origseq:
                    blacklist.add(origseq)
                outputHandle.loc[uid] = entry
                outputHandle.to_csv(outputFile, sep='\t', header=False)
                outputHandle.drop([uid], inplace=True)

  # STEP 5. Go through entire list and remove REFSEQ/regular duplicates
  ### TO DO: Read in the entire outputFile again and loop through the blacklist tuples. For every tuple, see if both (!) the REFSEQ accession number as well as the regular accession number are in the outFile. If yes, remove the line of the the regular accession number because it is a duplicate of the REFSEQ accession number line.
    # -> don't need to check if both numbers are in outFile. blacklist will only contain the regular accession number if the RefSeq accession was just processed (i.e. both numbers are in the outFile)
    # -> in case the regular number is NOT in the outFile, nothing needs to be done.
        outputHandle = pd.read_csv(outfn, sep='\t', index_col=1, encoding='utf-8') # Using ACCESSION column as index to make accessing it more efficient
        for entry in blacklist:
            # Attempting to drop regular duplicate. If the duplicate doesn't exist, nothing happens (except for a log message)
            try:
                outputHandle.drop([entry], inplace=True)
            except:
                log.info("Could not find accession %s when trying to remove it." % str(entry))
        with open(outfn, "w") as outputFile:
            # Replace existing outputFile content with updated list
            outputHandle.to_csv(outputFile, sep='\t', header=True)


########
# MAIN #
########

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    parser.add_argument("-o", "--outfn", type=str, required=True, help="path to output file")
    parser.add_argument("-q", "--query", type=str, required=False, default="Magnoliophyta[ORGN] AND 00000180000[SLEN] : 00000200000[SLEN] AND complete genome[TITLE] AND (chloroplast[TITLE] OR plastid[TITLE])", help="(Optional) Entrez query that will replace the standard query")
    args = parser.parse_args()
    main(args.outfn, args.query)
