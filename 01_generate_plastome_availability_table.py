#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
OBJECTIVE:
    This script generates a table of plastid genome records that are currently available on NCBI. It simultaneously ensures that no plastid genome is counted twice (issue about regular vs. RefSeq NC_ records). If an output file already exists, this script will identify the most recent date of any of the processed records and search NCBI only for younger records.

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
    10. Any note if a REFSEQ accession number exists and to which regular accession number it is equal to.
    11. the taxonomy

DESIGN:
    There are thousands of plastid genome sequences on GenBank. The parsing of the records is, thus, conducted one by one, not all simultaneously. Specifically, a list of unique identifiers is first obtained and then this list is looped over.

TO DO:
    * none for now

NOTES:
    * none for now

'''

#####################
# IMPORT OPERATIONS #
#####################
import xml.etree.ElementTree as ET
import os.path, subprocess, calendar
import pandas as pd
import argparse, sys
import coloredlogs, logging
from datetime import date, datetime

###############
# AUTHOR INFO #
###############
__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>, '\
             'Tilman Mehl <tilmanmehl@zedat.fu-berlin.de>'
__copyright__ = 'Copyright (C) 2019 Michael Gruenstaeudl and Tilman Mehl'
__info__ = 'Collect summary information on all plastid sequences stored ' \
           'in NCBI GenBank'
__version__ = '2020.02.18.1230'

#############
# DEBUGGING #
#############
import ipdb
# ipdb.set_trace()

#############
# FUNCTIONS #
#############

def getQueriedUIDs(query, mindate):
    '''
    Gets a list of all UIDs found by the search query via ESearch and EFetch.

    Args:
        query (str): a string that specifies the query term for NCBI
        mindate (date): a date that specifies the beginning of the date range from which to fetch records
    Returns:
        A list of UIDs as integers that is returned by the query.
    '''

  # STEP 1. Get list of UIDs
    if mindate:
        # TO DO: resolve issue around choice of mindate/maxdate:
        # Assuming maxdate = date.now() (date at which the script is run):
        # if mindate = newest_date (in table) it will fetch records that already exist in the table
            # -> since there is a very slight chance that refseq and original were both created on the same day it will add the original on the next run but won't remove it since the refseq that adds the original to the blacklist won't be processed
        # if mindate = newest_date+1 (day after newest date) it will skip records that were created on newest_date but after the script ran.
        # Possible solutions:
            # Rewrite blacklist creation so it parses ALL existing records in the table for refseqs (not just the ones that are added by the current script run), then removes the originals
                # -> decreases performance (maybe negligible?)
            # Set maxdate = date.now()-1 (day before the script is run) and mindate = newest_date+1 (day after newest_date)
                # -> no risk of getting existing entries or skipping entries
                # -> script does not fetch the very newest records that were added at the date of execution (is working with 1 day old data detrimental?)

        # Answer by MG: You are making a very good point! I would rather err on the side of poor performance and have a cumulative blacklist output file that is evaluated in its entirety every time the script is executed. That way we don't have to deal with a maxdate, which would only make things more complicated.

        esearchargs = ['esearch', '-db', 'nucleotide', '-sort', '"Date Released"', '-mindate', mindate.strftime("%Y/%m/%d"), '-maxdate', date.today().strftime("%Y/%m/%d"), '-query', query]
    else:
        esearchargs = ['esearch', '-db', 'nucleotide', '-sort', '"Date Released"', '-query', query]
    esearch = subprocess.Popen(esearchargs, stdout=subprocess.PIPE)
    efetchargs = ["efetch", "-db", "nucleotide", "-format", "uid"]
    efetch = subprocess.Popen(efetchargs, stdin=esearch.stdout, stdout=subprocess.PIPE)
    out, err = efetch.communicate()
    return list(map(int, out.splitlines()))[::-1]


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
    # No need to add UID as a separate field here since it will be added by pandas as an index later
    root = ET.fromstring(out)
    uidData = root.find("GBSeq")
    fields = []
    accession = uidData.find("GBSeq_primary-accession").text
    fields.append(accession)
    fields.append((uidData.find("GBSeq_accession-version").text).split('.')[1])
    fields.append(uidData.find("GBSeq_organism").text)
    fields.append(uidData.find("GBSeq_length").text)

    # Parse and format the date that the record was first online
    month_map = {"JAN":"01", "FEB":"02", "MAR":"03", "APR":"04", "MAY":"05", "JUN":"06", "JUL":"07", "AUG":"08", "SEP":"09", "OCT":"10", "NOV":"11", "DEC":"12"}
    create_date = uidData.find("GBSeq_create-date").text.split('-')
    fields.append(create_date[2] + "-" + month_map[create_date[1]] + "-" + create_date[0])

    # Parse all info related to the authors and the publication
    references = uidData.find("GBSeq_references").findall("GBReference")
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
                authstring = authstring + author.text + ","
            authstring = authstring[:-2]
            break
    fields.append(authstring)
    fields.append(title)
    fields.append(citation)

    # Parse comment field for RefSeq
    note = ""
    duplseq = None
    if accession[:3] == "NC_":
        comments = uidData.find("GBSeq_comment").text
        for comment in comments.split(";"):
            keyw1 = "PROVISIONAL REFSEQ: This record has not yet been subject to final NCBI review"
            keyw2 = "The reference sequence is identical to"
            if keyw1 in comment or keyw2 in comment:
                duplseq = comment.split(" ")[-1][:-1]
                note = "The REFSEQ accession `%s` is identical to accession `%s`." % (accession, duplseq)
    fields.append(note)

    fields.append(uidData.find("GBSeq_taxonomy").text)

    return fields, accession, duplseq


def main(outfn, query, update_only):

  # STEP 1. Set up logger
    log = logging.getLogger(__name__)
    coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level='DEBUG', logger=log)


  # STEP 2. Check if output file already exists, read existing UIDs, infer mindate
    UIDs_alreadyProcessed = []
    mindate = None
    if os.path.isfile(outfn) and sum(1 for line in open(outfn) if line.rstrip()) > 1: # If outfile exists and if outfile contains more than a single line (i.e., the title line)
        with open(outfn, "r") as outputFile:
            UIDs_alreadyProcessed = [row.split('\t')[0] for row in outputFile] # Read UID column
            UIDs_alreadyProcessed = list(map(int, UIDs_alreadyProcessed[1:])) # Discard header and convert UIDs to integer values
            log.info("Summary file `%s` already exists. Number of UIDs read: %s" % (str(outfn), str(len(UIDs_alreadyProcessed))))
            if(update_only):
                try: # Trying to set creation date of newest record as starting date for search query. If the file exists but there are no records in it (just the headers), it will skip.
                    outputFile.seek(0) # Return to beginning of file. Otherwise the next line will always throw an exception
                    mindate = datetime.strptime(outputFile.readlines()[-1].split('\t')[5], '%Y-%m-%d')
                    log.info("NOTE: Only records more recent than `%s` are being looked for." % (str(mindate)))
                except:
                    log.error("Could not read newest date from existing summary file.")
                    raise Exception
    else:
        with open(outfn, "w") as outputFile:
            log.info(("Summary file `%s` does not exist; generating new file. Thus, no UIDs read." % (str(outfn))))
            outputFile.write("UID\tACCESSION\tVERSION\tORGANISM\tSEQ_LEN\tCREATE_DATE\tAUTHORS\tTITLE\tREFERENCE\tNOTE\tTAXONOMY\n")


  # STEP 3. Get all existing UIDs and calculate which to be processed
    try:
        UIDs_allExisting = getQueriedUIDs(query, mindate)
    except Exception as err:
        log.exception("Error while retrieving UID list: " + str(err))
    log.info(("Number of UIDs on NCBI: %s" % (str(len(UIDs_allExisting)))))
    UIDs_notYetProcessed = []
    if len(UIDs_alreadyProcessed) > 0:
        try:
            UIDs_notYetProcessed = [uid for uid in UIDs_allExisting if uid not in UIDs_alreadyProcessed]
        except Exception as e:
            log.error(("Calculation of complement set of UIDs not successful: %s" % (e.message)))
            raise Exception
    else:
        UIDs_notYetProcessed = UIDs_allExisting
    log.info(("Number of UIDs to be processed: %s" % (str(len(UIDs_notYetProcessed)))))


  # STEP 4. Set up output handle, obtain info from new UIDs, and save them as lines in output file.
    blacklistfn = os.path.join(os.path.dirname(outfn), os.path.basename(outfn) + ".duplicates")
    if len(UIDs_notYetProcessed) > 0:
      # Load format of the summary file
      # This only throws an exception for me when outfn exists and is completely empty (i.e. not even headers)
      # if it has headers but is otherwise empty it should work fine
        outputHandle = pd.read_csv(outfn, nrows=0, sep='\t', index_col=0, encoding='utf-8')
      # Append to outfile, but also save duplicate ones to blacklist
        with open(outfn, "a") as outputFile:
            for uid in UIDs_notYetProcessed:
                log.info(("Reading and parsing UID `%s`, writing to `%s`." % (str(uid), str(outfn))))
                # Parse data
                try:
                    parsedUIDdata, accession, duplseq = getEntryInfo(uid) # Note: This is where the heavy-lifting is done!
                except Exception as err:
                    log.exception("Error retrieving info for UID " + str(uid) + ": " + str(err) + "\nSkipping this accession.")
                    continue
                # Write data to summary file
                outputHandle.loc[uid] = parsedUIDdata
                outputHandle.to_csv(outputFile, sep='\t', header=False)
                outputHandle.drop([uid], inplace=True)
                # Write refseq and duplseq to blacklist file
                if duplseq:
                    with open(blacklistfn, "a") as blfile:
                        blfile.write("%s\t%s\n" % (str(accession), str(duplseq)))

  # STEP 5. Go through entire list and remove accession numbers that are duplicates of REFSEQs
    if os.path.isfile(blacklistfn) and sum(1 for line in open(blacklistfn) if line.rstrip()) > 0:
        with open(blacklistfn, "r") as blfile:
            blacklist = [tuple(line.strip("\n").split("\t")) for line in blfile.readlines()]
        outputHandle = pd.read_csv(outfn, sep='\t', index_col=0, encoding='utf-8')

        for refseq, duplseq in blacklist:
            try: # Attempting to drop regular accession duplicate.
                outputHandle.drop(outputHandle.loc[outputHandle['ACCESSION'] == duplseq].index, inplace=True)
                log.info("Removed accession `%s` from `%s` because it is a duplicate of REFSEQ `%s`." % (duplseq, str(os.path.basename(outfn)), refseq))
            except: # If the duplicate doesn't exist, nothing happens (except for a log message)
                log.info("Could not find accession `%s` when trying to remove it." % str(duplseq)) # Reverted this to log level "info" since there is nothing wrong with a duplicate not existing (anymore)
        with open(outfn, "w") as outputFile:
            outputHandle.to_csv(outputFile, sep='\t', header=True) # Replace existing outputFile content with updated list


########
# MAIN #
########

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    parser.add_argument("-o", "--outfn", type=str, required=True, help="path to output file")
    parser.add_argument("-q", "--query", type=str, required=False, default="Magnoliophyta[ORGN] AND 170000:210000[SLEN] AND complete genome[TITLE] AND (chloroplast[TITLE] OR plastid[TITLE]) NOT unverified[TITLE] NOT partial[TITLE] AND 2000/01/01:2019/12/31[PDAT]", help="(Optional) Entrez query that will replace the standard query")
    parser.add_argument("-u", "--update_only", action="store_true", required=False, default=False, help="Will only add entries with more recent creation date than the most recent existing entry.")
    args = parser.parse_args()
    main(os.path.abspath(args.outfn), args.query, args.update_only)
