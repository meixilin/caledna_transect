# notes on creating a tree based phyla

1. download all phyla available at the time when building the crux database. 

```shell
/u/project/rwayne/software/Crux/crux_db/TAXO
```

to the final data folder: 

```shell
./final_data/TAXO
```

2. the taxid for all phylum available were collected and subjected to ncbi common tree generator
https://www.ncbi.nlm.nih.gov/Taxonomy/CommonTree/wwwcmt.cgi

the ranks of 
superkingdom
kingdom
phylum 
were selected and a phylip tree was downloaded. 

3. we compared the taxonomy coverage of different primers by parsing the available variables 

# name changes in current data 

> 09/14/2019 vs jan 2018

| Previous name            | Current name                      | Previous rank   | Current rank |
| ------------------------ | --------------------------------- | --------------- | ------------ |
| Chromerida               | Colpodellida                      | phylum          | no rank      |
| Haplosporidia            | Haplosporida                      | phylum          | order        |
| candidate division WPS-2 | Candidatus Eremiobacteraeota      | phylum          | phylum       |
| candidate division ZB3   | Candidatus Marinamargulisbacteria | phylum          | class        |
| Candidatus Acetothermia  | Candidatus Bipolaricaulota        | phylum          | phylum       |
| Entorrhizomycota         | N.A.                              | phylum as fungi |              |
|                          |                                   |                 |              |
|                          |                                   |                 |              |
|                          |                                   |                 |              |

```bash
linmeixideMBP:~ linmeixi$ diff /Users/linmeixi/UCLA/Lab/abiotic_transect/derive_data/step1_mk_phyloseq/phyla/phyliptree_no0.phy /Users/linmeixi/UCLA/Lab/abiotic_transect/derive_data/step1_mk_phyloseq/phyla/phyliptree_no0_oldname.phy 
9c9
< Colpodellida:4,
---
> Chromerida:4,
45a46
> Entorrhizomycota:4,
63,64c64,65
< 'Candidatus Tectomicrobia':4,
< 'Candidatus Melainabacteria':4,
---
> Candidatus_Tectomicrobia:4,
> Candidatus_Melainabacteria:4,
67c68
< 'candidate division NC10':4,
---
> candidate division NC10:4,
70c71
< 'Candidatus Cloacimonetes':4,
---
> Candidatus_Cloacimonetes:4,
87c88
< 'Candidatus Bipolaricaulota':4,
---
> Candidatus_Acetothermia:4,
```








