# Title: query variable colour 
# Author: Meixi Lin (meixilin@ucla.edu)
# Date: Mon Apr 20 00:55:21 2020

# def functions --------
query_var_color <- function(groupv) {
    load("./derive_data/step4_beta_diver/util/catlist_pal.RData")
    var_pal = catlist_pal[[groupv]]
    return(var_pal)
}

# # example use: 
# dummy = data.frame(c("Others", "Aquatic"), 1:2)
# colnames(dummy) = c("x", "value")
# test = ggplot(dummy, aes(x = x, y = value, color = x)) +
#     geom_point() +
#     scale_color_manual(values = query_var_color('majorhab'))
