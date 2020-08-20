import re
import os.path
import argparse
import coloredlogs, logging
import time
import sys
#import PlastomeIntegrityChecks as pic
import fetchpubmed
import entrezpy.conduit
from ete3 import NCBITaxa
from pathlib import Path
from datetime import datetime

# For suppressing console output
import io
from contextlib import redirect_stdout

###############
# AUTHOR INFO #
###############
__author__ = 'Michael Gruenstaeudl <m.gruenstaeudl@fu-berlin.de>, '\
             'Tilman Mehl <tilmanmehl@zedat.fu-berlin.de>'
__copyright__ = 'Copyright (C) 2019 Michael Gruenstaeudl and Tilman Mehl'
__info__ = 'Create or append a list of species names that are proven to lack one or more inverted repeats in their plastid genome'
__version__ = '2020.08.20.1800'

#############
# DEBUGGING #
#############
import ipdb
# ipdb.set_trace()

def fetch_pubmed_articles(mail, query):
    articles = None
    cond = entrezpy.conduit.Conduit(mail)
    fetch_pubmed = cond.new_pipeline()
    sid = fetch_pubmed.add_search({'db': 'pubmed', 'term': query, 'rettype': 'count'})
    fid = fetch_pubmed.add_fetch({'retmode':'xml'}, dependency=sid, analyzer=fetchpubmed.PubMedAnalyzer())
    a = cond.run(fetch_pubmed)
    result = a.get_result()
    if result.size() >= 1:
        articles = result.pubmed_records
    return articles

'''
def fetch_pmc_articles(mail, query):
    articles = None
    cond = entrezpy.conduit.Conduit(mail)
    fetch_pmc = cond.new_pipeline()
    sid = fetch_pmc.add_search({'db': 'pmc', 'term': query, 'rettype': 'count'})
    fid = fetch_pmc.add_fetch({'retmode':'xml'}, dependency=sid, analyzer=fetchpmc.PMCAnalyzer())
    a = cond.run(fetch_pmc)
    result = a.get_result()
    if result.size() >= 1:
        articles = result.pmc_records
    return articles
'''

def get_abstract_text(article):
    text = ""
    for elem in article.xml.findall(".//AbstractText"):
        text += elem.text + "\n"
    return text

def get_article_keywords(article):
    keywords = []
    for elem in article.xml.findall(".//Keyword"):
        keywords.append(elem.text)
    return keywords

def get_species_from_pubmed_article(article, ncbi):
    '''
    Parses a pubmed article for species names
    '''
    species = set()
    abstract_text = get_abstract_text(article)
    keywords = get_article_keywords(article)
    for keyword in keywords:
        abstract_text += keyword + "\n"
    ncbi_query_results = ncbi.get_name_translator(construct_species_query(abstract_text))
    for name, id in ncbi_query_results.items():
        # We're only interested in species-rank taxons
        if ncbi.get_rank(id)[id[0]] == "species":
            species.add(name)

    return species

def construct_species_query(text):
    '''
    Constructs a string for every two adjacent words in the text and returns them as a list.
    '''
    species_queries = []
    taxonomic_expressions = ["sp.", "x", "var.", "subsp.", "f."]
    regexp = r'[^a-zA-Z0-9.\-:\s]'
    words = text.split()
    for i in range(0,len(words)-2):
        # Since not all species names are simply two words, we need to account for these cases.
        # Below is an attempt to account for the cases that were detected by a brief look at existing names
        # We're also removing all characters that cannot(should not?) occur in taxonomy (e.g. a species name might be in parentheses)
        if words[i+1] in taxonomic_expressions:
            # "sp." is preceded by one word and is usually followed by one word, but we also found a case with two following words
            if words[i+1] == "sp.":
                if i < len(words)-3:
                    species_queries.append(re.sub(regexp, '', words[i] + " " + words[i+1] + " " + words[i+2]))
                    if i < len(words)-4:
                        species_queries.append(re.sub(regexp, '', words[i] + " " + words[i+1] + " " + words[i+2] + " " + words[i+3]))
            if words[i+1] == "var." or words[i+1] == "f.":
                # "var." and "f." are preceded by two words and followed by one word
                if i < len(words)-3 and i > 0:
                    species_queries.append(re.sub(regexp, '', words[i-1] + " " + words[i] + " " + words[i+1] + " " + words[i+2]))
            if words[i+1] == "subsp.":
                # "subsp." is preceded by two words and is usually followed by one word, but we also found a case with two following words
                if i < len(words)-3 and i > 0:
                    species_queries.append(re.sub(regexp, '', words[i-1] + " " + words[i] + " " + words[i+1] + " " + words[i+2]))
                    if i < len(words)-4:
                        species_queries.append(re.sub(regexp, '', words[i-1] + " " + words[i] + " " + words[i+1] + " " + words[i+2] + " " + words[i+3]))
            if words[+1] == "x":
                # "x" is preceded by two words and followed by two words
                if i < len(words)-4 and i > 0:
                    species_queries.append(re.sub(regexp, '', words[i-1] + " " + words[i] + " " + words[i+1] + " " + words[i+2] + " " + words[i+3]))
        else:
            species_queries.append(re.sub(regexp, '', words[i] + " " + words[i+1]))
    return species_queries

def get_irl_clade_species(ncbi):
    species_ids = []
    irl_clade_tree = ncbi.get_topology([ncbi.get_name_translator(['IRL clade'])['IRL clade'][0]])
    for leaf in irl_clade_tree.iter_leaves():
         species_ids.append(int(leaf.name))
    species = set(ncbi.translate_to_names(species_ids))
    return species

def read_blacklist(fp_blacklist):
    '''
    Read a file of blacklisted genera/species.
    Params:
     - fp_blacklist: file path to input file
    '''
    blacklist = set()
    with open(fp_blacklist, "r") as fh_blacklist:
        for line in [l.strip() for l in fh_blacklist.readlines()]:
            if not line.startswith("#"):
                blacklist.add(line)
    return blacklist

def append_blacklist(fp_blacklist, blacklist):
    with open(fp_blacklist, "a") as fh_blacklist:
        fh_blacklist.write("# Update on %s\n" % datetime.now().strftime("%Y-%m-%d, %H:%M"))
        for entry in blacklist:
            fh_blacklist.write(entry + "\n")

def main(args):

    ## STEP 1. Initialize variables
    # TODO: replace static mail address (maybe with a simple username@hostname ?)
    mail = "tilmanmehl@fu-berlin.de"
    ncbi = NCBITaxa()
    blacklist = set()
    blacklist_existing = set()

    ## STEP 2. Read blacklist if the file already exists
    if os.path.isfile(args.file_blacklist):
        print("Reading existing blacklist ...")
        blacklist_existing = read_blacklist(args.file_blacklist)

    ## STEP 3. Assemble species names of IRL clade of Fabaceae
    print("\nFetching species names of taxa in 'IRL clade' of Fabaceae ...")
    irl_clade_species = get_irl_clade_species(ncbi)
    print("  Adding new species names to blacklist ...")
    blacklist = irl_clade_species.difference(blacklist_existing)

    ## STEP 4. Assemble species names of active search
    print("\nSearching for matching PubMed articles ...")
    # Suppressing console output by entrezpy.conduit.Conduit.new_pipeline()
    # Source: https://codingdose.info/2018/03/22/supress-print-output-in-python/
    trap = io.StringIO()
    with redirect_stdout(trap):
        articles = fetch_pubmed_articles(mail, args.query)
    print("  Number of matching articles found: %s" % str(len(articles)))

    print("\nIterating through matching articles ...")
    for article in articles:
        print("\n  Fetching species names of article: %s" % str(article.xml.findall(".//ArticleTitle")[0].text))
        article_species = get_species_from_pubmed_article(article, ncbi)
        print("  Number of species names parsed: %s" % str(len(article_species)))
        blacklist = blacklist.union(article_species.difference(blacklist_existing))

    ## STEP 5. Append only new taxon names to blacklist
    print("\nCalculating and appending species names not previously in blacklist ...")
    append_blacklist(args.file_blacklist, blacklist)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="  --  ".join([__author__, __copyright__, __info__, __version__]))
    parser.add_argument("-f", "--file_blacklist", type=str, required=True, help="path to blacklist file")
    parser.add_argument("-q", "--query", type=str, required=False, default="inverted[TITLE] AND repeat[TITLE] AND loss[TITLE]", help="query used to fetch PMC articles that will be scanned for species with missing IRs")
    args = parser.parse_args()
    main(args)
