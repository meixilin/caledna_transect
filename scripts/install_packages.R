# install cran based packages --------
if(!require(devtools)) {install.packages("devtools")}
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}
if (!requireNamespace("reshape2", quietly = TRUE)) {install.packages("reshape2")}
if (!requireNamespace("ggpubr", quietly = TRUE)) {install.packages("ggpubr")}
if (!requireNamespace("FSA", quietly = TRUE)) {install.packages("FSA")}
if (!requireNamespace("taxize", quietly = TRUE)) {install.packages("taxize")}
if (!requireNamespace("taxonomizr", quietly = TRUE)) {install.packages("taxonomizr")}
if (!requireNamespace("pls", quietly = TRUE)) {install.packages("pls")}
if (!requireNamespace("zetadiv", quietly = TRUE)) {install.packages("zetadiv")}
if (!requireNamespace("leaflet", quietly = TRUE)) {install.packages("leaflet")}
if (!requireNamespace("adespatial", quietly = TRUE)) {install.packages("adespatial")}
if (!requireNamespace("iNEXT", quietly = TRUE)) {install.packages("iNEXT")}
if (!requireNamespace("pals", quietly = TRUE)) {install.packages("pals")}
if (!requireNamespace("usedist", quietly = TRUE)) {install.packages("usedist")}
if (!requireNamespace("gstat", quietly = TRUE)) {install.packages("gstat")}
if (!requireNamespace("rgdal", quietly = TRUE)) {install.packages("rgdal")}
if (!requireNamespace("rasterVis", quietly = TRUE)) {install.packages("rasterVis")}
if (!requireNamespace("ggcorrplot", quietly = TRUE)) {install.packages("ggcorrplot")}
if (!requireNamespace("geosphere", quietly = TRUE)) {install.packages("geosphere")}
if (!requireNamespace("ggmap", quietly = TRUE)) {install.packages("ggmap")}
if (!requireNamespace("ggedit", quietly = TRUE)) {install.packages("ggedit")}
if (!requireNamespace("rgeos", quietly = TRUE)) {install.packages("rgeos")}
if (!requireNamespace("ggsn", quietly = TRUE)) {install.packages("ggsn")}
if (!requireNamespace("mapview", quietly = TRUE)) {install.packages("mapview")}
if (!requireNamespace("leaflet.extras", quietly = TRUE)) {install.packages("leaflet.extras")}

# install github based packages --------
# ranacapa for downstream analysis: 
require(devtools)
devtools::install_github("gauravsk/ranacapa")
# microbiomeSeq
require(devtools)
devtools::install_github("umerijaz/microbiomeSeq")
# depending on your previous packages, you might need to install these
BiocManager::install("impute")
BiocManager::install("preprocessCore")
BiocManager::install("GO.db")
BiocManager::install("AnnotationDbi")
BiocManager::install("DESeq2")
BiocManager::install("WGCNA")

# install gradient forest --------
if (!requireNamespace("randomForest", quietly = TRUE)) {install.packages("randomForest")}
if (!requireNamespace("gradientForest", quietly = TRUE)) {install.packages("gradientForest", repos="http://R-Forge.R-project.org")}

# # if shown error: 
# make: gfortran: No such file or directory
# make: *** [rfsub.o] Error 1
# ERROR: compilation failed for package ‘extendedForest’

# # install gfortran as recommended in CRAN:https://cran.r-project.org/bin/macosx/tools/
# if needed, you need to modify the Makeconf file to specify the path of gfortran executable installed in this way, For example: 
# /Library/Frameworks/R.framework/Resources/etc/Makeconf line 51, 
# FC = /usr/local/gfortran/bin/gfortran


