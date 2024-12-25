library(ggplot2)

file = "Insert path to phenotype/heritability information file here"

data = read.csv(file, header = F) # Loading in the data
names(data) = c("Tick", "Heritability", "Mean phenotype") # Naming the columns


# coeff is used to scale the axes of the plot appropriately
coeff = max(data$Heritability)/max(data$`Mean phenotype`) 


# Plots heritability and mean phenotype on one graph
ggplot(data, aes(Tick)) + 
  geom_line(aes(y = Heritability/coeff, colour = "He")) + 
  geom_line(aes(y = `Mean phenotype`, colour = "Phe"))+
  scale_color_manual(values = c(He = "black", Phe = "red"))+
  scale_y_continuous(
    name = "Mean Phenotype",
    sec.axis = sec_axis(~.*coeff, name = "Heritability")
  )+
  labs(y = "Mean Phenotype",
       x = "Generations",
  ) +
  theme_bw()+
  theme(axis.title.y.left = element_text(color = "red"),
        legend.position = "none",
        axis.title.x = element_text(size = 18),
        axis.title.y.right = element_text(size = 18),
        axis.title.y = element_text(size = 18))