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
    esearchargs = ['esearch', '-db', 'nucleotide', '-sort', '"Date Released"', '-query', query]
    esearch = subprocess.Popen(esearchargs, stdout=subprocess.PIPE)

    ### TO DO: Due to the date sorting in the esearchargs, the UIDs are sorted with the newest first. However, we want the oldest first. This could be done after esearch or after efetch, but may have different effects!

    efetchargs = ["efetch", "-db", "nucleotide", "-format", "uid"]
    efetch = subprocess.Popen(efetchargs, stdin=esearch.stdout, stdout=subprocess.PIPE)
    out, err = efetch.communicate()
    
    ### TO DO: Maybe I can just reverse the out list? Does that mess things up?
    ### return out.splitlines().reverse()


  ## TO DO: Please move STEP2 to main()
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
                authstring = authstring + author.text.replace(","," ") + ", "
                #authstring = authstring[:-2]
            break
    fields.append(authstring)
    fields.append(title)
    fields.append(citation)

    ### TO DO: Please also parse out the field "COMMENT" and check if there is the following information:
    '''
    LOCUS       NC_027250             195251 bp    DNA     circular PLN 18-JUN-2015
    DEFINITION  Carex siderosticta plastid, complete genome.
    ...
    COMMENT     PROVISIONAL REFSEQ: This record has not yet been subject to final
                NCBI review. The reference sequence is identical to KP751906.
    '''
    ### If a REFSEQ (accession number always starts with "NC_") has been created, then the regular accession (here: KP751906) is a duplicate and must be removed from the output list. To generate a safe removal mechanism, please save both the REFSEQ accession number (here: NC_027250) and the regular accession (here: KP751906) as a tuple and return to main(). In main, we must collect these tuples in a blacklist and remove any duplicates.
    
    ### TO DO: If a case of REFSEQ/regular accession number exists, this shall be written to fields as a note (fields.append(note)) with the following language: "The reference sequence NC_027250 is identical to record KP751906."

    return fields


def main(outfn, query):

  # STEP 1. Set up logger
    log = logging.getLogger(__name__)
    coloredlogs.install(fmt='%(asctime)s [%(levelname)s] %(message)s', level='DEBUG', logger=log)

  # STEP 2. Check if output file already exists, read existing UIDs
    if os.path.isfile(outfn):
        with open(outfn, "r") as outputFile:
            pass
            ### TO DO: Parse out the UIDs of the outputFile, save to list `UIDs_alreadyProcessed`
            ### log.info(("Summary file `%s` already exists. %s UIDs read." % (str(outfn), str(len(UIDs_alreadyProcessed)))))

    if not os.path.isfile(outfn):
        with open(outfn, "w") as outputFile:
            outputFile.write("UID\tACCESSION\tVERSION\tORGANISM\tSEQ_LEN\tCREATE_DATE\tAUTHORS\tTITLE\tREFERENCE\n")
            log.info(("Summary file `%s` does not exist; generating new file. Thus, no UIDs read." % (str(outfn))))
            UIDs_alreadyProcessed = []

  # STEP 3. Get all existing UIDs and calculate which to be processed
    UIDs_allExisting = getNewUIDs(query, outfn)
    log.info(("Number of unique UIDs corrently on NCBI: `%s`" % (str(len(UIDs_allExisting)))))

    ### TO DO: Move here from getNewUIDs(). Specifically, remove the entries of UIDs_alreadyProcessed from the list UIDs_allExisting, saving the remainder as list `UIDs_notYetProcessed`
    ### log.info(("Number of UIDs to be processed: `%s`" % (str(len(UIDs_notYetProcessed)))))
    
  # STEP 4. Obtain info from new UIDs and save as lines in output file.
    # Load format of the summary file
    outputHandle = pd.read_csv(outfn, nrows=0, sep='\t', index_col=0, encoding='utf-8')
    # Append to outfile
    with open(outfn, "a") as outputFile:
        blacklist = []
        for uid in UIDs_notYetProcessed:
            log.info(("Reading and parsing UID %s, writing to %s" % (str(uid), str(outfn))))
            outputHandle.loc[uid] = getEntryInfo(uid)  # Note: This is where the heavy-lifting is done!

            ### TO DO: Next to the full UID info, we may (not always) also receive a tuple per UID to be saved in a blacklist. This blacklist should be appended to after every new UID being processed. This blacklist is used later in step 5.
            ### blacklist.append(getEntryInfo(uid)[2]) # Or something like that!

            outputHandle.to_csv(outputFile, sep='\t', header=False)
            outputHandle.drop([uid], inplace=True)

  # STEP 5. Go through entire list and remove REFSEQ/regular duplicates
    ### TO DO: Read in the entire outputFile again and loop through the blacklist tuples. For every tuple, see if both (!) the REFSEQ accession number as well as the regular accession number are in the outFile. If yes, remove the line of the the regular accession number because it is a duplicate of the REFSEQ accession number line.


########
# MAIN #
########

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    parser.add_argument("-o", "--outfn", type=str, required=True, help="path to output file")
    parser.add_argument("-q", "--query", type=str, required=False, default="Magnoliophyta[ORGN] AND 00000180000[SLEN] : 00000200000[SLEN] AND complete genome[TITLE] AND (chloroplast[TITLE] OR plastid[TITLE])", help="(Optional) Entrez query that will replace the standard query")
    args = parser.parse_args()
    main(args.outfn, args.query)
