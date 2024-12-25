setwd("C:\\Users\\kushv\\genArchSelect\\QTL_prototype\\output_files\\QTL_V2")
file = "QTL_V2_HaplotypeCount1241212155.csv"

library(dplyr)
library(tidyverse)
library(stringr)

# Functions ----

# Function to move nth row into the colName
colHeader <- function(df, n) {
  names(df) <- as.character(unlist(df[n, ]))
  df = df[-c(1:n), , drop = F]
}


# Function to select two loci in a haplotype matrix and return the frequency
# The frequency needs to be the last column
HaploCount <- function(x, i, j){
  temp = cbind(x[, i], x[, j])
  
  Bool = x[,i]&x[,j]
  Freqs = x[,length(colnames(x))]
  return(sum(Freqs[Bool]))
  
  
}


HaploCount <- function(HaploFreq){
  
  stopifnot(is.data.frame(HaploFreq))
  
  
  NumLoci = length(HaploFreq)-1
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

# Code -----

lines = readLines(file, n = -1) # Read File
data = as.data.frame(strsplit(lines, split = ","))

# Splitting data into a list based on tick number, which is first row value
Ticks = unique(unlist(data[1,]))

SplitData = vector(mode = "list", length = length(Ticks))
names(SplitData) = Ticks

for (tick in 1:length(Ticks)){
  SplitData[[tick]] = data[,data[1,] == Ticks[tick]]
}

SplitData = lapply(SplitData, colHeader, 2)
rownames(SplitData) = NULL


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
  ElimCheck = rowSums(Haplotype[,-ncol(Haplotype)])>=2
  FilteredHaplotypes[[i]] = drop_na(Haplotype[ElimCheck, ])
}

HaploFreqs =lapply(FilteredHaplotypes, function(df) HaploCompile(df, HaploCount(df)))
names(HaploFreqs) = Ticks

for (i in 1:length(HaploFreqs)){
  diag(HaploFreqs[[i]]) = NA
}







# Old Code ---------------------

# data = colHeader(data, 2)
# rownames(data) = NULL
# dataDim = dim(data)


# data_args = c(data, sep ="") # Dataframe to List
# Haplotypes= as.data.frame(table(do.call(paste, data_args)))# Puts the haplotypes together into a frequency Table

#Split the haplotype's loci into seperarte cells and convert to Bool
# TempSplit = str_split_fixed(Haplotypes$Var1, "", dataDim[2] )
# TempSplit = as.data.frame(ifelse(TempSplit == "T", T, F))
# Haplotypes = cbind(TempSplit, Haplotypes$Freq)

#Save Dimensions
#HaploDim = dim(Haplotypes)

# ElimHaplo removes Haplotypes which have only one mutation present, that is the sum 
# is equal to 1, as those haplotypes won't contribute to any frequencies as they have only one 
# locus
#ElimCheck = (rowSums(Haplotypes[, -HaploDim[2]]))>=2
#ElimHaplo = drop_na(Haplotype[ElimCheck, ])

#FinalHaplotypeFrequencies = HaploCompile(ElimHaplo, HaploCount(ElimHaplo))
#diag(FinalHaplotypeFrequencies) = NA



