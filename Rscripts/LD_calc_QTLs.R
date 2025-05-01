library(dplyr)
library(tidyverse)
library(data.table)

# Functions -------

# Function to move nth row into the colName
colHeader <- function(df, n) {
  names(df) <- as.character(unlist(df[n, ]))
  df = df[-c(1:n), , drop = F]
}


# Function to select two loci in a haplotype matrix and return the frequency
# The frequency needs to be the last column

HaploCount <- function(HaploFreq){
  
  stopifnot(is.data.frame(HaploFreq))
  
  
  NumLoci = length(HaploFreq)
  FinalValues = c()
  
  for (i in 1:NumLoci){
    for (j in (i):(NumLoci)){
      Bool = HaploFreq[,i]&HaploFreq[,j]
      Freqs = HaploFreq[, length(HaploFreq)]
      FinalFreq = sum(Freqs[Bool])
      FinalValues = c(FinalValues, FinalFreq)
    }
  }
  return(FinalValues)
}

# Frequencies from HaploCount
HaploCompile <- function(HaploHisto, Frequencies){
  NumLoci = length(HaploHisto)-1
  FinalMatrix = matrix(nrow = NumLoci, ncol = NumLoci)
  vec = 1:NumLoci
  for (i in NumLoci:1){
    
    #Compliment value of i with resepect to the dimensions of the matrix
    i_compliment = NumLoci + 1 - i
    
    # Fills the column with values
    FinalMatrix[,i_compliment][i_compliment:NumLoci] = Frequencies[vec[1]:vec[i]]
    
    # Fills the row with values
    FinalMatrix[i_compliment, ][i_compliment:NumLoci] = Frequencies[vec[1]:vec[i]]
    
    # Increment someValue to next set of values
    vec = vec + length(vec)
    
    # Removes the last value as the length should decrease
    vec = vec[-i]
    
  }
  return(FinalMatrix)
  
}

# Reading in haplo information
LoadHaplo = function(HaploFile){
  lines = readLines(HaploFile, n = -1) # Read File
  data = as.data.frame(strsplit(lines, split = ","))
  
  # Splitting data into a list based on tick number, which is first row value
  Ticks = unique(unlist(data[1,]))
  
  SplitData = vector(mode = "list", length = length(Ticks))
  names(SplitData) = Ticks
  
  for (tick in 1:length(Ticks)){
    SplitData[[tick]] = data[,data[1,] == Ticks[tick], drop = F]
  }
  
  SplitData = lapply(SplitData, colHeader, 2)
  rownames(SplitData) = NULL
  
  # SplitData[[1]][,colnames(SplitData[[1]]) %in% positive$Mutation]
  
  
  SplitHaplotypes = vector(mode = "list", length = length(Ticks))
  for (i in 1:length(SplitData)){
    data_args = c(SplitData[[i]], sep ="") # Dataframe to List
    SplitHaplotypes[[i]] = as.data.frame(table(do.call(paste, data_args)))# Puts the haplotypes together into a frequency Table
  }
  
  #Split the haplotype's loci into seperarte cells and convert to Bool
  for (i in 1:length(SplitHaplotypes)){
    dataDim = dim(SplitHaplotypes[[i]])
    
    TempSplit = str_split_fixed(SplitHaplotypes[[i]]$Var1, "", nchar(as.character(SplitHaplotypes[[i]][1,1])))
    TempSplit = as.data.frame(ifelse(TempSplit == "T", T, F))
    SplitHaplotypes[[i]] = cbind(TempSplit, SplitHaplotypes[[i]]$Freq)
  }
  
  # ElimHaplo removes Haplotypes which have only one mutation present, that is the sum 
  # is equal to 1, as those haplotypes won't contribute to any frequencies as they have only one 
  # locus
  FilteredHaplotypes = vector(mode = "list", length = length(Ticks))
  for (i in 1:length(SplitHaplotypes)){
    Haplotype = SplitHaplotypes[[i]]
    ElimCheck = rowSums(Haplotype[,-ncol(Haplotype), drop = F])>=2
    FilteredHaplotypes[[i]] = drop_na(Haplotype[ElimCheck, ])
  }
  
  HaploFreqs =lapply(FilteredHaplotypes, function(df) HaploCompile(df, HaploCount(df)))
  names(HaploFreqs) = Ticks
  
  return(HaploFreqs)
}

# Reading in allele information
LoadAllele = function(AlleleFile){
  data = read.table(AlleleFile, sep = ",")
  colnames(data) = c("Tick", "Mutation", "Freq")
  
  Ticks = unique(unlist(data[,1]))
  SplitData = vector(mode = "list", length = length(Ticks))
  
  for (i in 1:length(Ticks)){
    SplitData[[i]] = data[data$Tick==Ticks[i],]  
  }
  
  names(SplitData) = Ticks
  
  AlleleFreqs = SplitData
  
  return(AlleleFreqs)
}

#Reading in Phenotype information 
LoadPheno = function(PhenoData){
  data = read.csv(PhenoData, header = F)
  names(data) = c("Tick", "Heritability", "Mean phenotype")
  return(data)
}

# Reading in the positions of the Mutations. This information was dervied from 
# a vcf file ----

MutID = read.csv("MID-POS-map_high.csv")
colnames(MutID) = c("MutID", "Pos", "EffectSize")

# Obtaining Mutation ID
AlleleFreqs = lapply(AlleleFreqs, function(x){
  x$MutID = as.numeric(gsub(".*<(\\d+):.*", "\\1", x$Mutation))
  return(x)
})

# Merging data using MutID as a dict
PosData = merge(AlleleFreqs[[1]], MutID, by = "MutID", all.x = T)
PosData = na.omit(PosData)
PosData = filter(PosData, PosData$Freq>0.1)


# Deriving chromosome number
PosData$Chr = (PosData$Pos %/% 1000000) + 1



MutsByChr <- PosData%>%
  group_by(Chr) %>%
  summarise(Mutations = paste(Mutation, collapse = ","))

Chromosomes = vector('list', 10)

# Chromosomes is a list containing mutations split by group
for (i in 1:10){
  a = (MutsByChr[i,2])
  Chromosomes[[i]] = (unlist(strsplit(as.character(a), split = ",")))
}

for (i in 1:10){
  print(dim(SplitData[[1]][Chromosomes[[i]]])[2])
}

# Same chromsosome heatmap ----

# Random sampling 
Haplo_smp_same = sample(SplitData[[1]][,Chromosomes[[2]]], 10)
Allele_smp_same = PosData[PosData$Mutation %in% colnames(Haplo_smp_same), ]

# Sorting to make sure it is the same order
Haplo_smp_same = Haplo_smp_same[,order(colnames(Haplo_smp_same))]
Allele_smp_same = Allele_smp_same[order(Allele_smp_same$Mutation), ]


# Calculating LD

SplitHaplotypes_smp_same = vector(mode = "list", length = length(Ticks))

data_args = c(Haplo_smp_same, sep ="") # Dataframe to List
Haplo_smp_same = as.data.frame(table(do.call(paste, data_args)))# Puts the haplotypes together into a frequency Table

#Split the haplotype's loci into seperarte cells and convert to Bool
  dataDim = dim(Haplo_smp_same)
  
  TempSplit = str_split_fixed(Haplo_smp_same$Var1, "", nchar(as.character(Haplo_smp_same[1,1])))
  TempSplit = as.data.frame(ifelse(TempSplit == "T", T, F))
  Haplo_smp_same = cbind(TempSplit, Haplo_smp_same$Freq)

FilteredHaplotypes_smp_same = vector(mode = "list", length = length(Ticks))

  Haplotype = Haplo_smp_same
  ElimCheck = rowSums(Haplotype[,-ncol(Haplotype), drop = F])>=2
  FilteredHaplotypes_smp_same = drop_na(Haplotype[ElimCheck, ])


HaploFreqs_smp_same =HaploCompile(FilteredHaplotypes_smp_same, HaploCount(FilteredHaplotypes_smp_same))


N = dim(SplitData[[1]])[1]
freqs = list()
pApB = list()
pAB = list()
Dmax = list()



  pAB = HaploFreqs_smp_same/(N)
  
  
  
  freqs = t(Allele_smp_same$Freq)
  colnames(freqs) = Allele_smp_same$Mutation
  
  MutNum = length(freqs)
  pApB = matrix(nrow = MutNum, ncol = MutNum)
  
  #[[i]][x] rep pA and [[i]][y] rep pB
  for (x in 1:MutNum){
    for (y in 1:MutNum){
      pApB[x, y] = freqs[x]*freqs[y]
      
    }
  }
  
  


D = pApB



D= pAB -pApB




Dmax = D


for(i in 1:length(Ticks)){
  for(x in 1:length(freqs[[i]])){
    for (y in 1:length(freqs[[i]])) {
      if (D[x, y]>= 0){
        Dmax[x, y] = min((freqs[[i]][x]*(1-freqs[[i]][y])), ((1 - freqs[[i]][x])*freqs[[i]][y]))
      } else {
        Dmax[x, y] = min((freqs[[i]][x]*freqs[[i]][y]), ((1-freqs[[i]][x])*(1-freqs[[i]][y])))
      }
    }
  }
}


Dprime = D

for(i in 1:length(Ticks)){
  Dprime[[i]] = D[[i]]/Dmax[[i]]
}



colnames(Dprime) = Allele_smp_same$Pos
rownames(Dprime) = Allele_smp_same$Pos

Dprime = as.data.frame(Dprime)
Dprime = rownames_to_column(Dprime, var = "POS_A")
DprimeLong = pivot_longer(cols = -POS_A, data = Dprime, names_to = "POS_B", values_to = "DPRIME")
DprimeLong = filter(DprimeLong, DprimeLong$DPRIME <=1)
DprimeLong = filter(DprimeLong, DprimeLong$DPRIME>= -1)
same = DprimeLong


a = ggplot(same, aes(POS_A, POS_B, fill = DPRIME))+
  geom_tile()+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-1, 1))+
  theme(axis.text.x = element_text(angle = 75, hjust = 1, vjust = 0.5))


# Diff Chromosomes heatmap ----

# Random sampling 
Haplo_smp_diff = sample(PosData$Mutation, 10)
Haplo_smp_diff = SplitData[[1]][,Haplo_smp_diff]
Allele_smp_diff = PosData[PosData$Mutation %in% colnames(Haplo_smp_diff), ]

# Sorting to make sure it is the same order
Haplo_smp_diff = Haplo_smp_diff[,order(colnames(Haplo_smp_diff))]
Allele_smp_diff = Allele_smp_diff[order(Allele_smp_diff$Mutation), ]


# Calculating LD

SplitHaplotypes_smp_diff = vector(mode = "list", length = length(Ticks))

data_args = c(Haplo_smp_diff, sep ="") # Dataframe to List
Haplo_smp_diff = as.data.frame(table(do.call(paste, data_args)))# Puts the haplotypes together into a frequency Table

#Split the haplotype's loci into seperarte cells and convert to Bool
dataDim = dim(Haplo_smp_diff)

TempSplit = str_split_fixed(Haplo_smp_diff$Var1, "", nchar(as.character(Haplo_smp_diff[1,1])))
TempSplit = as.data.frame(ifelse(TempSplit == "T", T, F))
Haplo_smp_diff = cbind(TempSplit, Haplo_smp_diff$Freq)

FilteredHaplotypes_smp_diff = vector(mode = "list", length = length(Ticks))

Haplotype = Haplo_smp_diff
ElimCheck = rowSums(Haplotype[,-ncol(Haplotype), drop = F])>=2
FilteredHaplotypes_smp_diff = drop_na(Haplotype[ElimCheck, ])


HaploFreqs_smp_diff =HaploCompile(FilteredHaplotypes_smp_diff, HaploCount(FilteredHaplotypes_smp_diff))


N = dim(SplitData[[1]])[1]
freqs = list()
pApB = list()
pAB = list()
Dmax = list()



pAB = HaploFreqs_smp_diff/(N)



freqs = t(Allele_smp_diff$Freq)
colnames(freqs) = Allele_smp_diff$Mutation

MutNum = length(freqs)
pApB = matrix(nrow = MutNum, ncol = MutNum)

#[[i]][x] rep pA and [[i]][y] rep pB
for (x in 1:MutNum){
  for (y in 1:MutNum){
    pApB[x, y] = freqs[x]*freqs[y]
    
  }
}




D = pApB



D= pAB -pApB




Dmax = D


for(i in 1:length(Ticks)){
  for(x in 1:length(freqs[[i]])){
    for (y in 1:length(freqs[[i]])) {
      if (D[x, y]>= 0){
        Dmax[x, y] = min((freqs[[i]][x]*(1-freqs[[i]][y])), ((1 - freqs[[i]][x])*freqs[[i]][y]))
      } else {
        Dmax[x, y] = min((freqs[[i]][x]*freqs[[i]][y]), ((1-freqs[[i]][x])*(1-freqs[[i]][y])))
      }
    }
  }
}


Dprime = D

for(i in 1:length(Ticks)){
  Dprime[[i]] = D[[i]]/Dmax[[i]]
}



colnames(Dprime) = Allele_smp_diff$Pos
rownames(Dprime) = Allele_smp_diff$Pos

Dprime = as.data.frame(Dprime)
Dprime = rownames_to_column(Dprime, var = "POS_A")
DprimeLong = pivot_longer(cols = -POS_A, data = Dprime, names_to = "POS_B", values_to = "DPRIME")
DprimeLong = filter(DprimeLong, DprimeLong$DPRIME <=1)
DprimeLong = filter(DprimeLong, DprimeLong$DPRIME>= -1)
diff = DprimeLong



b = ggplot(diff, aes(POS_A, POS_B, fill = DPRIME)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-1, 1)) +
  geom_tile(data = subset(diff, POS_A == POS_B), 
            aes(POS_A, POS_B), 
            fill = "black") +
  theme(axis.text.x = element_text(angle = 75, hjust = 0.5, vjust = 0.5))



save_plot('10k_10_low_same.png', plot = a, dpi = 900, base_height = 7, base_width = 7)
save_plot('10k_10_diff.png', plot = b, dpi = 900, base_height = 7, base_width = 7)

