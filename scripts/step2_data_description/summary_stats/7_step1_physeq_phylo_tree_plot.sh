#!/bin/bash 

# this script was used to generate a phylogenetic tree from the taxonomy dump file

# define variables

TXDMP=/u/project/rwayne/software/Crux/crux_db/TAXO
WORK=/u/project/rwayne/meixilin/abiotic_transect
OUTPUT=${WORK}/derive_data/step1_mk_phyloseq

cd ${WORK}
# find all the "phylum" in the taxonomy dump file 
grep -w "phylum" ${TXDMP}/nodes.dmp > ${OUTPUT}/all_phylum_20190119.txt

# notes on the field of all_phylum_20190119.txt file 
# nodes.dmp file consists of taxonomy nodes. The description for each node includes the following
# fields:
#         tax_id                                  -- node id in GenBank taxonomy database
#         parent tax_id                           -- parent node id in GenBank taxonomy database
#         rank                                    -- rank of this node (superkingdom, kingdom, ...)
#         embl code                               -- locus-name prefix; not unique
#         division id                             -- see division.dmp file
#         inherited div flag  (1 or 0)            -- 1 if node inherits division from parent
#         genetic code id                         -- see gencode.dmp file
#         inherited GC  flag  (1 or 0)            -- 1 if node inherits genetic code from parent
#         mitochondrial genetic code id           -- see gencode.dmp file
#         inherited MGC flag  (1 or 0)            -- 1 if node inherits mitochondrial gencode from parent
#         GenBank hidden flag (1 or 0)            -- 1 if name is suppressed in GenBank entry lineage
#         hidden subtree root flag (1 or 0)       -- 1 if this subtree has no sequence data yet
#         comments                                -- free-text comments and citations

# notes on taxonomy files 
# Taxonomy names file (names.dmp):                
#         tax_id                                  -- the id of node associated with this name
#         name_txt                                -- name itself
#         unique name                             -- the unique variant of this name if name not unique
#         name class                              -- (synonym, common name, ...)

# now find the name of all the nodes 

