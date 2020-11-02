# define colour scheme 
mycolours <- pals::polychrome(n = 6)[-1]
myname =  c("unknown", "Aquatic", "Herbaceous", "Shrub", "Tree")
names(mycolours) <- myname

rand <- data.frame(myname, 1:5, 1:5)
rand <- rand[c(2:5,1),]
colnames(rand) = c("Habitat", "Longitude", "Latitude")
rand$Habitat <- factor(rand$Habitat, levels = rand$Habitat, ordered = T)


pp <- ggplot(data = rand, aes(x = Longitude, y = Latitude, fill = Habitat)) + 
    geom_point(shape = 23, colour = 'black', size = 2) + 
    scale_fill_manual(values = mycolours) +
    theme_bw()

pp
ggsave(filename = "./plots/step6_gradient_forest/species_response_leg.pdf", plot = pp, width = 4, height = 4, device = "pdf")

# define sample 
type <- data.frame(c("Detected", "Not_detected"), 1:2, 1:2)
colnames(type) = c("Sample", "Longitude", "Latitude")
pp2 <- ggplot(data = type, aes(x = Longitude, y = Latitude, color = Sample, shape = Sample)) + 
    geom_point(size = 2) +
    scale_color_manual(values = c("black", "slategrey")) + 
    scale_shape_manual(values = c(23,16)) + 
    theme_bw()
ggsave(filename = "./plots/step6_gradient_forest/species_response_leg2.pdf", plot = pp2, width = 4, height = 4, device = "pdf")


