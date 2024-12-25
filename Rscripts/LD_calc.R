library(dplyr)
library(tidyverse)

setwd("C:\\Users\\kushv\\genArchSelect\\QTL_prototype\\output_files\\QTL_V2\\N_1k\\2Mil")
HaploFile = ""



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





# Obtaining Haplo Frequency -------------------

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

# for (i in 1:length(HaploFreqs)){
#   diag(HaploFreqs[[i]]) = NA
# }






# Obtaining Allele Frequency ----
AlleleFile = "1092303715383681055_QTL_2_10k.csv"
data = read.table(AlleleFile, sep = ",")
colnames(data) = c("Tick", "Mutation", "Freq")

Ticks = unique(unlist(data[,1]))
SplitData = vector(mode = "list", length = length(Ticks))

for (i in 1:length(Ticks)){
  SplitData[[i]] = data[data$Tick==Ticks[i],]  
}

names(SplitData) = Ticks

AlleleFreqs = SplitData

AlleleFreqs = lapply(AlleleFreqs, function(x){
  x$SelCoeff = as.numeric(gsub(".*:(-?[0-9.]+)>.*", "\\1", x$Mutation))
  x$WeightedContribution = x$SelCoeff * x$Freq
  return(x)
})

Expected = unlist(lapply(AlleleFreqs, function(x)2*mean(x$WeightedContribution)))
ggplot(as.data.frame(Expected), aes(x = Ticks, y = Expected)) + geom_line()
# Calculating D' ----

N = sum(SplitHaplotypes[[1]]$`SplitHaplotypes[[i]]$Freq`)
freqs = list()
pApB = list()
pAB = list()
Dmax = list()
for (i in 1:length(Ticks)){
  pAB[[i]] = HaploFreqs[[i]]/2000 
  
  
  
  freqs[[i]] = t(AlleleFreqs[[i]]$Freq)
  colnames(freqs[[i]]) = AlleleFreqs[[i]]$Mutation
  
  MutNum = length(freqs[[i]])
  pApB[[i]] = matrix(nrow = MutNum, ncol = MutNum)
  
  #[[i]][x] rep pA and [[i]][y] rep pB
  for (x in 1:MutNum){
    for (y in 1:MutNum){
      pApB[[i]][x, y] = freqs[[i]][x]*freqs[[i]][y]
      
    }
  }
  
  
}

D = pApB
names(D) = Ticks

for(i in 1:length(Ticks)){
  D[[i]] = pAB[[i]]-pApB[[i]]
}



Dmax = D
names(Dmax) = Ticks

for(i in 1:length(Ticks)){
  for(x in 1:length(freqs[[i]])){
    for (y in 1:length(freqs[[i]])) {
      if (D[[i]][x, y]>= 0){
        Dmax[[i]][x, y] = min((freqs[[i]][x]*(1-freqs[[i]][y])), ((1 - freqs[[i]][x])*freqs[[i]][y]))
      } else {
        Dmax[[i]][x, y] = min((freqs[[i]][x]*freqs[[i]][y]), ((1-freqs[[i]][x])*(1-freqs[[i]][y])))
      }
    }
  }
}


Dprime = D

for(i in 1:length(Ticks)){
 Dprime[[i]] = D[[i]]/Dmax[[i]]
}
colnames(x) = AlleleFreqs[[]]$Mutation

for (i in 1:length(Ticks)){
  colnames(Dprime[[i]]) = AlleleFreqs[[i]]$Mutation
  rownames(Dprime[[i]]) = AlleleFreqs[[i]]$Mutation
}


# Final Data ----

# Haplotype Frequency
HaploFreqs

# Allele Frequency
AlleleFreqs

# Some Math Constants
pAB

pApB

D

Dmax

Dprime
