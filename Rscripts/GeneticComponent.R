# 10,000 (scroll for other pop size)
setwd("C:/Users/kushv/genArchSelect/QTL_finalModel_SURC/output_files/N_10k_L_10e7/")
files = list.files()

TenThousandGeneticComponent = data.frame(Burnin = c(1, 2, 3), Gen = c(0, 0, 0),
                                         Phe = c(0, 0, 0))

j = 1
for (i in c(1, 4, 7)){
  AlleleFile = files[i]
  file = files[i+1]
  
  data = read.table(AlleleFile, sep = ",") # loading in the data
  colnames(data) = c("Tick", "Mutation", "Freq") # naming column headers
  
  Ticks = unique(unlist(data[,1])) #vector of ticks (generations)
  
  # Pre allocates an empty list, with the length of the number of ticks
  SplitData = vector(mode = "list", length = length(Ticks))
  
  # for loop that essentially groups all mutations of a specific generation 
  # into one element of SplitData
  for (i in 1:length(Ticks)){
    SplitData[[i]] = data[data$Tick==Ticks[i],]  
  }
  
  names(SplitData) = Ticks
  
  # AlleleFreqs is the final dataset
  AlleleFreqs = SplitData
  
  # Calculate effect size and weighted contribution
  AlleleFreqs = lapply(AlleleFreqs, function(x){
    x$SelCoeff = as.numeric(gsub(".*:(-?[0-9.]+)>.*", "\\1", x$Mutation))
    x$WeightedContribution = x$SelCoeff * x$Freq
    return(x)
  })
  
  TenThousandGeneticComponent[j,2]=(sum(AlleleFreqs[[length(AlleleFreqs)]]$WeightedContribution))
  
  data = read.csv(file, header = F) # Loading in the data
  names(data) = c("Tick", "Heritability", "Mean phenotype") # Naming the columns
  
  TenThousandGeneticComponent[j,3]=data[nrow(data),3]
  j = j+1
  }


# 1000
setwd("C:/Users/kushv/genArchSelect/QTL_finalModel_SURC/output_files/N_1k_L_10e7/")
files = list.files()

ThousandGeneticComponent = data.frame(Burnin = c(1, 2, 3), Gen = c(0, 0, 0),
                                         Phe = c(0, 0, 0))

j = 1
for (i in c(1, 4, 7)){
  AlleleFile = files[i]
  file = files[i+1]
  
  data = read.table(AlleleFile, sep = ",") # loading in the data
  colnames(data) = c("Tick", "Mutation", "Freq") # naming column headers
  
  Ticks = unique(unlist(data[,1])) #vector of ticks (generations)
  
  # Pre allocates an empty list, with the length of the number of ticks
  SplitData = vector(mode = "list", length = length(Ticks))
  
  # for loop that essentially groups all mutations of a specific generation 
  # into one element of SplitData
  for (i in 1:length(Ticks)){
    SplitData[[i]] = data[data$Tick==Ticks[i],]  
  }
  
  names(SplitData) = Ticks
  
  # AlleleFreqs is the final dataset
  AlleleFreqs = SplitData
  
  # Calculate effect size and weighted contribution
  AlleleFreqs = lapply(AlleleFreqs, function(x){
    x$SelCoeff = as.numeric(gsub(".*:(-?[0-9.]+)>.*", "\\1", x$Mutation))
    x$WeightedContribution = x$SelCoeff * x$Freq
    return(x)
  })
  
  ThousandGeneticComponent[j,2]=(sum(AlleleFreqs[[length(AlleleFreqs)]]$WeightedContribution))
  
  data = read.csv(file, header = F) # Loading in the data
  names(data) = c("Tick", "Heritability", "Mean phenotype") # Naming the columns
  
  ThousandGeneticComponent[j,3]=data[nrow(data),3]
  j = j+1
}


# 100
setwd("C:/Users/kushv/genArchSelect/QTL_finalModel_SURC/output_files/N_100_L_10e7/")
files = list.files()

HundredGeneticComponent = data.frame(Burnin = c(1, 2, 3), Gen = c(0, 0, 0),
                                         Phe = c(0, 0, 0))

j = 1
for (i in c(1, 4, 7)){
  AlleleFile = files[i]
  file = files[i+1]
  
  data = read.table(AlleleFile, sep = ",") # loading in the data
  colnames(data) = c("Tick", "Mutation", "Freq") # naming column headers
  
  Ticks = unique(unlist(data[,1])) #vector of ticks (generations)
  
  # Pre allocates an empty list, with the length of the number of ticks
  SplitData = vector(mode = "list", length = length(Ticks))
  
  # for loop that essentially groups all mutations of a specific generation 
  # into one element of SplitData
  for (i in 1:length(Ticks)){
    SplitData[[i]] = data[data$Tick==Ticks[i],]  
  }
  
  names(SplitData) = Ticks
  
  # AlleleFreqs is the final dataset
  AlleleFreqs = SplitData
  
  # Calculate effect size and weighted contribution
  AlleleFreqs = lapply(AlleleFreqs, function(x){
    x$SelCoeff = as.numeric(gsub(".*:(-?[0-9.]+)>.*", "\\1", x$Mutation))
    x$WeightedContribution = x$SelCoeff * x$Freq
    return(x)
  })
  
  HundredGeneticComponent[j,2]=(sum(AlleleFreqs[[length(AlleleFreqs)]]$WeightedContribution))
  
  data = read.csv(file, header = F) # Loading in the data
  names(data) = c("Tick", "Heritability", "Mean phenotype") # Naming the columns
  
  HundredGeneticComponent[j,3]=data[nrow(data),3]
  j = j+1
}

