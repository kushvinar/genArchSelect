library(ggplot2)
library(dplyr)
library(tidyverse)
library(cowplot)

setwd("C:\\Users\\kushv\\genArchSelect\\QTL_finalModel_SURC/output_files/N_10k_L_10e7/withoutScale/")
file = "3264789426475091028QTL_V2_He+Pheno_10k.csv"
AlleleFile = "7503075833640078399_QTL_3_10k.csv"
FitFunc = 'N = 10k'

data = read.csv(file, header = F)
names(data) = c("Tick", "Heritability", "Mean phenotype")

coeff = max(data$Heritability)/max(data$`Mean phenotype`)


# Plots heritability and mean phenotype on one graph
ggplot(data, aes(Tick)) + 
  geom_line(aes(y = Heritability/coeff, colour = "He")) + 
  geom_line(aes(y = `Mean phenotype`, colour = "Phe"))+
  # geom_line(data = as.data.frame(Expected), aes(x = Ticks, y = Expected, colour = 'Exp') )+
  scale_color_manual(values = c(He = "black", Phe = "red", Exp = "Blue"))+
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

# Plotting segregating sites. Use list `AlleleFreqs` from LD_calc.R
ggplot(data.frame(Ticks = as.numeric(Ticks), Sites = unlist(lapply(AlleleFreqs, nrow))), aes(Ticks, Sites))+
  geom_point()+
  labs(title = FitFunc, y = "Number of QTLs") + 
  theme_bw()

# Plotting fixed QTLs from AlleleFreqs
FixedQTLs = data.frame(Ticks = Ticks, FixedQTLs = as.numeric(lapply(AlleleFreqs, function(x){
  return(sum(x$Freq == 1))
})))

ggplot(FixedQTLs, aes(x = Ticks, y = FixedQTLs)) +
  geom_line()+
  theme_bw()

FixedQTLs[FixedQTLs$FixedQTLs != 0, ]
FixedQTLs = FixedQTLs[FixedQTLs$Ticks>= 16500 & FixedQTLs$Ticks<= 17500,]

44325


# Plotting allele Frequencies

CombinedFrequencies = data.table::rbindlist(AlleleFreqs)

CombinedFrequencies %>%
  ggplot( aes(x = Tick,
              y = Freq,
              group = Mutation)) +
  geom_line() +
  theme_bw()

a = CombinedFrequencies %>%
  group_by(Mutation)%>%
  filter(n()>=500)

CombinedFrequencies %>%
  ggplot( aes(x = Tick,
              y = Freq,
              group = Mutation)) +
  geom_line()

# Delta p

data = pivot_wider(CombinedFrequencies, names_from = Tick, values_from = Freq)

sum(data$`64250` == 1, na.rm = T)

data = data[,-c(2,3)]

data = data[(data[,2] != 1), ]

deltaPvalues = data.frame(Mutation = data$Mutation)

for(i in 3:ncol(data)){
  deltaPvalues[,i-1] = data[i]-data[i-1]
}



deltaPvalues_long = pivot_longer(deltaPvalues, !Mutation, names_to = "ticks", values_to = "values")

ggplot(deltaPvalues_long, aes(x = ticks, 
                              y = values)) +
         geom_boxplot(na.rm = T) + 
  theme_bw()+
  theme(axis.text.x = element_text(size = 12,
                                    angle = 90),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 18))+
  labs(x = "Generations",
       y = "Chnage in Allele Frequency")
  


SurcPlotter = function(PhenoFile, AlleleFile, FitFunc = "10k",FileName, Filter = NULL){
  
  
  # Plots heritability and mean phenotype on one graph
  Phenodata = read.csv(PhenoFile, header = F)
  names(Phenodata) = c("Tick", "Heritability", "Mean phenotype")
  
  coeff = max(Phenodata$Heritability)/max(Phenodata$`Mean phenotype`)
  
  if(!is.null(Filter)){
    Phenodata = Phenodata[Filter,]
  }
  
  
  a = ggplot(Phenodata, aes(Tick)) + 
    geom_line(aes(y = Heritability/coeff, colour = "He")) + 
    geom_line(aes(y = `Mean phenotype`, colour = "Phe"))+
    scale_color_manual(values = c(He = "black", Phe = "red"))+
    scale_y_continuous(
      name = "Mean Phenotype",
      sec.axis = sec_axis(~.*coeff, name = "Heritability")
    )+
    labs(y = "Mean Phenotype",
         x = "Generations") +
    theme_bw()+
    theme(axis.title.y.left = element_text(color = "red"),
          legend.position = "none",
          axis.title.x = element_text(size = 18),
          axis.title.y.right = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     size = 12),
          axis.text.y.left = element_text(size = 12),
          axis.text.y.right = element_text(size = 12))
  
  
  
  Alleledata = read.table(AlleleFile, sep = ",")
  colnames(Alleledata) = c("Tick", "Mutation", "Freq")
  
  Ticks = unique(unlist(Alleledata[,1]))
  SplitData = vector(mode = "list", length = length(Ticks))
  
  for (i in 1:length(Ticks)){
    SplitData[[i]] = Alleledata[Alleledata$Tick==Ticks[i],]  
  }
  
  names(SplitData) = Ticks
  
  AlleleFreqs = SplitData
  
  if(!is.null(Filter)){
    AlleleFreqs = AlleleFreqs[Filter]
  }
  
  
  # Plotting allele Frequencies
  
  CombinedFrequencies = data.table::rbindlist(AlleleFreqs)
  
  b = CombinedFrequencies %>%
    ggplot( aes(x = Tick,
                y = Freq,
                group = Mutation)) +
    geom_line() +
    theme_bw() + 
    labs(x = "Generations",
         y = "Allele frequencies")+
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     size = 12),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          axis.text.y = element_text(size = 12))
  
  
  
  
  Wide_Alleledata = pivot_wider(CombinedFrequencies, names_from = Tick, values_from = Freq)
  
  if(any(Wide_Alleledata[,2] == 1, na.rm = T)){
    Wide_Alleledata = Wide_Alleledata[-which(Wide_Alleledata[,2] == 1),]
  }
  print(Wide_Alleledata)  

  
  deltaPvalues = data.frame(Mutation = Wide_Alleledata$Mutation)
  
  if(grepl("1k", FitFunc)){
    generationwidth=10
    col_index = 2
    
    for(i in seq(2+generationwidth, ncol(Wide_Alleledata), generationwidth)){
      deltaPvalues[,col_index] = Wide_Alleledata[i]-Wide_Alleledata[i-generationwidth]
      col_index = col_index + 1
    }
  } else{
    for(i in 3:ncol(Wide_Alleledata)){
      deltaPvalues[,i-1] = Wide_Alleledata[i]-Wide_Alleledata[i-1]
    }
  }
    
  
  
  deltaPvalues_long = pivot_longer(deltaPvalues, !Mutation, names_to = "ticks", values_to = "values")
  
  print(length(unique(deltaPvalues_long$ticks)))
  
  deltaPvalues_long$ticks = as.numeric(deltaPvalues_long$ticks)
  
  deltaPvalues_long = deltaPvalues_long %>%
                        arrange(ticks)
  
  Ticks = unique(deltaPvalues_long$ticks)
  BreakTicks = Ticks[seq(1, length(Ticks)-1, length.out = 8)]
  
    
  c = ggplot(deltaPvalues_long, aes(x = factor(ticks), 
                                 y = values)) +
    geom_boxplot(na.rm = T,
                 outlier.size = 0.6)+
    labs(x = "Generations",
         y = paste("\u394", "Allele Frequency")) +
    scale_x_discrete(breaks = BreakTicks)+
    theme_bw() + 
     theme(axis.text.x  = element_text(angle = 90,
                                      hjust = 1),
           axis.title.x = element_text(size = 18),
           axis.title.y = element_text(size = 18), 
           axis.text.y = element_text(size = 12))
     
  
  
  
  plt = plot_grid(a, b, c, ncol = 3)
  plot(plt)
  
  save_plot(FileName, plot = plt, dpi = 900, base_height = 9, base_width = 16)
}

deltaPvalues_long$ticks[seq(1, length(deltaPvalues_long$ticks), length.out = 8)]
scale_x_discrete(breaks = seq(min(deltaPvalues_long$ticks), max(deltaPvalues_long$ticks), length.out = 7))+
  
  seq(min(deltaPvalues_long$ticks), max(deltaPvalues_long$ticks), length.out = 7)

deltaPvalues_long$ticks[seq(1, length(deltaPvalues_long$ticks), length.out = 8)]
deltaPvalues_long$ticks[c(T,rep(F,Xscale))]


# Code for batchwise plotting
Phe = c("QTL_V2_He+pheno_BurnIn2MilReal_5197473836330009195.csv",
        "QTL_V2_He+pheno_BurnIn2MilReal_5690401104725564838.csv",
        "QTL_V2_He+pheno_BurnIn2MilReal_5838701142826674457.csv",
        "QTL_V2_He+pheno_BurnIn2MilReal_6522260778554470372.csv",
        "QTL_V2_He+pheno_BurnIn2MilReal_6553982241615284217.csv",
        "QTL_V2_He+pheno_BurnIn2MilReal_7069034506788787739.csv",
        "QTL_V2_He+pheno_BurnIn2MilReal_9131510556339839757.csv")

QTLs = c("QTL_V2_5197473836330009195.csv",
         "QTL_V2_5690401104725564838.csv",
         "QTL_V2_5838701142826674457.csv",
         "QTL_V2_6522260778554470372.csv",
         "QTL_V2_6553982241615284217.csv",
         "QTL_V2_7069034506788787739.csv",
         "QTL_V2_9131510556339839757.csv")

files = data.frame(Phe, QTLs)

apply(files, MARGIN = 1, function(x){
  PheFile = x[1]
  FitFunc = '1.0 + (dnorm(x, 5.0, 2)/scale)'
  
  
  
  data = read.csv(PheFile, header = F)
  names(data) = c("Tick", "Heritability", "Mean phenotype")
  
  coeff = max(data$Heritability)/max(data$`Mean phenotype`)
  
  
  # Plots heritability and mean phenotype on one graph
  plot = ggplot(data, aes(Tick)) + 
    geom_line(aes(y = Heritability/coeff, colour = "He")) + 
    geom_line(aes(y = `Mean phenotype`, colour = "Phe"))+
    scale_color_manual(values = c(He = "black", Phe = "red"))+
    scale_y_continuous(
      name = "Mean Phenotype",
      sec.axis = sec_axis(~.*coeff, name = "Heritability")
    )+
    labs(title = PheFile, y = "Mean Phenotype",
    ) +
    theme_bw()
})


apply(files, MARGIN = 1, function(x){
  AlleleFile = x[2]
  
  data = read.table(AlleleFile, sep = ",")
  colnames(data) = c("Tick", "Mutation", "Freq")
  
  Ticks = unique(unlist(data[,1]))
  SplitData = vector(mode = "list", length = length(Ticks))
  
  for (i in 1:length(Ticks)){
    SplitData[[i]] = data[data$Tick==Ticks[i],]  
  }
  
  names(SplitData) = Ticks
  
  AlleleFreqs = SplitData
  
  # Plotting segregating sites. Use list `AlleleFreqs` from LD_calc.R
  ggplot(data.frame(Ticks = Ticks, Sites = unlist(lapply(AlleleFreqs, nrow))), aes(Ticks, Sites))+
    geom_line()+
    labs(title = AlleleFile, y = "Number of QTLs")
  
})

apply(files, MARGIN = 1, function(x){
  AlleleFile = x[2]
  
  data = read.table(AlleleFile, sep = ",")
  colnames(data) = c("Tick", "Mutation", "Freq")
  
  Ticks = unique(unlist(data[,1]))
  SplitData = vector(mode = "list", length = length(Ticks))
  
  for (i in 1:length(Ticks)){
    SplitData[[i]] = data[data$Tick==Ticks[i],]  
  }
  
  names(SplitData) = Ticks
  
  AlleleFreqs = SplitData
  
  # Plotting fixed QTLs from AlleleFreqs
  FixedQTLs = data.frame(Ticks = Ticks, FixedQTLs = as.numeric(lapply(AlleleFreqs, function(x){
    return(sum(x$Freq == 1))
  })))
  
  ggplot(FixedQTLs, aes(x = Ticks, y = FixedQTLs)) +
    geom_line() + labs(title = AlleleFile)
})


# Code to plot Fitness Functions

# scale = dnorm(5.0, 5.0, 1.0)
fitness1 = function(x) 1.0 + dnorm(x, 5, 2)
fitness2 = function(x) 1.0 + dnorm(x, 5, 4)
# fitness3 = function(x) 1.0 + dnorm(x, 5, 0.5)/scale

par(mar = c(5, 6, 4, 1)+.1)
curve(fitness1, 0, 10,
      xlab = "Phenotype",
      ylab = "Relative fitness",
      cex.lab = 2,
      cex.axis = 1.5,
      lwd = 3)
curve(fitness2, 0, 10,
      lwd = 3,
      add = T, col = "red")
# curve(fitness1, 0, 10, add = T, col = "blue")


legend("topright", 
       legend = c("Low Selection", "High Selection"),
       col = c("red", "black"),
       pch = c(1),
       bty = "n",
       pt.cex = 1,
       cex = 2,
       text.col = "black",
       horiz = FALSE,)




  EverythingPlotter = function(PhenoFile, AlleleFile, FitFunc){
  data = read.csv(PhenoFile, header = F)
  names(data) = c("Tick", "Heritability", "Mean phenotype")
  
  coeff = max(data$Heritability)/max(data$`Mean phenotype`)
  
  
  # Plots heritability and mean phenotype on one graph
  a = ggplot(data, aes(Tick)) + 
    geom_line(aes(y = Heritability/coeff, colour = "He")) + 
    geom_line(aes(y = `Mean phenotype`, colour = "Phe"))+
    scale_color_manual(values = c(He = "black", Phe = "red"))+
    scale_y_continuous(
      name = "Mean Phenotype",
      sec.axis = sec_axis(~.*coeff, name = "Heritability")
    )+
    labs(title = FitFunc, y = "Mean Phenotype",
    ) +
    theme_bw()+
    theme(axis.title.y.left  = element_text(color = "red"),
          legend.position = "none")
  
  
  
  
  data = read.table(AlleleFile, sep = ",")
  colnames(data) = c("Tick", "Mutation", "Freq")
  
  Ticks = unique(unlist(data[,1]))
  SplitData = vector(mode = "list", length = length(Ticks))
  
  for (i in 1:length(Ticks)){
    SplitData[[i]] = data[data$Tick==Ticks[i],]  
  }
  
  names(SplitData) = Ticks
  
  AlleleFreqs = SplitData
  
  
  
  # Plotting segregating sites. Use list `AlleleFreqs` from LD_calc.R
  b = ggplot(data.frame(Ticks = Ticks, Sites = unlist(lapply(AlleleFreqs, nrow))), aes(Ticks, Sites))+
    geom_line()+
    labs(title = FitFunc, y = "Number of QTLs") + 
    theme_bw()
  
  # Plotting fixed QTLs from AlleleFreqs
  FixedQTLs = data.frame(Ticks = Ticks, FixedQTLs = as.numeric(lapply(AlleleFreqs, function(x){
    return(sum(x$Freq == 1))
  })))
  
  c = ggplot(FixedQTLs, aes(x = Ticks, y = FixedQTLs)) +
    geom_line()+
    theme_bw()
  

  
  
  # Plotting allele Frequencies
  
  CombinedFrequencies = data.table::rbindlist(AlleleFreqs)
  
  d = CombinedFrequencies %>%
    ggplot( aes(x = Tick,
                y = Freq,
                group = Mutation)) +
    geom_line() +
    theme_bw()
  
  
  plot_grid(a, b, c, d) + coord_fixed()
}




deltaPvalues = data.frame(Mutation = Wide_Alleledata$Mutation)

generationwidth=15
col_index = 2

for(i in seq(2+generationwidth, ncol(Wide_Alleledata), generationwidth)){
  deltaPvalues[,col_index] = Wide_Alleledata[i]-Wide_Alleledata[i-generationwidth]
  col_index = col_index + 1
  }



for (i in seq(13,130,counter) ){
  print(i)
  print(i-counter)
}

# Making fitness function graph
x_vals <- seq(-10, 20, length.out = 500)
fitness_data <- data.frame(
  x = rep(x_vals, 2),
  fitness = c(1.0 + dnorm(x_vals, 5, 2),
              1.0 + dnorm(x_vals, 5, 4)),
  selection = rep(c("High Selection", "Low Selection"), each = length(x_vals))
)

ggplot(fitness_data, aes(x = x, y = fitness, color = selection)) +
  geom_line(size = 1.5) +
  labs(x = "Phenotype", y = "Relative fitness", color = NULL) +
  scale_color_manual(values = c("High Selection" = "black", "Low Selection" = "red")) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "right",  # default is right inside — we'll fix that below
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    plot.margin = margin(20, 80, 20, 20)  # space for external legend
  ) +
  theme_bw()+
  theme(
  
    legend.justification = c("left", "top"),
    plot.margin = margin(20, 100, 20, 20)
  )

ggplot(fitness_data, aes(x = x, y = fitness, color = selection)) +
  geom_line(size = 1.5) +
  labs(x = "Phenotype", y = "Relative fitness") +
  scale_color_manual(values = c("High Selection" = "black", "Low Selection" = "red")) +
  theme_minimal(base_size = 16) +
  theme_bw()+
  theme(
    legend.position = c(0.95, 0.95),  # inside top right
    legend.justification = c("right", "top"),
    legend.background = element_rect(fill = alpha("white", 0.6), color = NA),
    legend.text = element_text(size = 18),
    legend.title = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18)
  )+
  scale_x_continuous(limits = c(-5, 15), breaks = seq(-5, 15, 1))



  