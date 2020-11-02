#!/bin/bash
WORK="/Users/linmeixi/Lab/caledna_transect"
cd ${WORK}
RSCRIPT=scripts/step2_data_description/concordance/2_describe_5_dup_sites_20200107.R
OTHERSCRIPT=scripts/step2_data_description/concordance/2_describe_5_dup_sites_full_20200107.R # considers full phylist 
TXLEVEL=("Family" "Genus" "Species")

for kk in ${TXLEVEL[@]}; do
    Rscript --vanilla ${RSCRIPT} "rare" "0" ${kk}
    Rscript --vanilla ${OTHERSCRIPT} "rare" "0" ${kk}
    Rscript --vanilla ${RSCRIPT} "deco" "2" ${kk}
    Rscript --vanilla ${OTHERSCRIPT} "deco" "2" ${kk}
    Rscript --vanilla ${RSCRIPT} "deco" "99" ${kk}
    Rscript --vanilla ${OTHERSCRIPT} "deco" "99" ${kk}
done

