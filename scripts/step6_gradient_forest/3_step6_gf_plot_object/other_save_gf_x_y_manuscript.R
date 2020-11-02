setwd("~/Lab/caledna_transect/")
load("/Users/linmeixi/UCLA/Lab/abiotic_transect/derive_data/step4_gradient_forest/2_final/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17.RData")
write.csv(x = gf$X, file = "./derive_data/step4_gradient_forest/2_final/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17_X.csv")
write.csv(x = gf$Y, file = "./derive_data/step4_gradient_forest/2_final/gf_deco_all_Family_Presence_2000_2_0.05_FALSE_17_Y.csv")
