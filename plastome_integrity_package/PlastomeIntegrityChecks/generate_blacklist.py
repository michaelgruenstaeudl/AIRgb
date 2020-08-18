import os.path, time
import sys
import PlastomeIntegrityChecks as pic
from ete3 import NCBITaxa
from pathlib import Path

if len(sys.argv) > 0:

    blacklist = set()
    if os.path.isfile(sys.argv[1]):
        print("Reading existing blacklist %s ..." % sys.argv[1])
        with open(sys.argv[1], "r") as fh_blacklist:
            for line in [l.rstrip() for l in fh_blacklist.readlines()]:
                if not line.startswith("#"):
                    blacklist.add(line)

    ncbi = NCBITaxa()
#    if (time.time() - os.path.getmtime(os.path.join(Path.home(), ".etetoolkit/taxa.sqlite"))) > 2628000:
#        print("Local NCBI taxonomy database is older than one month. Updating...")
#        ncbi.update_taxonomy_database()

    print("Getting taxonomy of 'IRL clade'...")
    species_ids = []
    irl_clade_tree = ncbi.get_topology([ncbi.get_name_translator(['IRL clade'])['IRL clade'][0]])
    for leaf in irl_clade_tree.iter_leaves():
         species_ids.append(int(leaf.name))
    species = set(ncbi.translate_to_names(species_ids))
    print("Merging 'IRL clade' species with blacklist...")
    blacklist = blacklist.union(species)

    print("Writing blacklist to %s." % sys.argv[1])
    with open(sys.argv[1], "w") as fh_blacklist:
        for species in blacklist:
            fh_blacklist.write(species + "\n")
