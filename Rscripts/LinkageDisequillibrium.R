HaplotFile = "Insert path to haplotype information file"
AlleleFile = "Insert path to allele information file"

# Obtaining the haplotype (pAB) frequencies


lines = readLines(HaploFile, n = -1) # Read File
data = as.data.frame(strsplit(lines, split = ",")) # convert to dataframe

# Splitting data into a list based on tick number, which is first row value
Ticks = unique(unlist(data[1,]))

SplitData = vector(mode = "list", length = length(Ticks)) # pre-allocating memory for a list
names(SplitData) = Ticks

# Putting the all data from each generation into the respective elemnt in SplitData
for (tick in 1:length(Ticks)){
  SplitData[[tick]] = data[,data[1,] == Ticks[tick], drop = F]
}

SplitData = lapply(SplitData, colHeader, 2)
rownames(SplitData) = NULL


SplitHaplotypes = vector(mode = "list", length = length(Ticks)) # pre-allocating memory

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

# ********************************************
# HaploFreqs is the final dataset that contains information about the frequency of each 
# haplotype present in the population
# The next section gets the allele frequencies
# ********************************************


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

AlleleFreqs = SplitData

# **********************************************
# AlleleFreqs is the final dataset of allele frequencies
# The next section does the calculations for Tajima's D and D'
# **********************************************

N = sum(SplitHaplotypes[[1]]$`SplitHaplotypes[[i]]$Freq`) # Obtaining N

# Pre-allocating
freqs = list()
pApB = list()
pAB = list()
Dmax = list()

# For each generation, we obtain pApB and pAB
for (i in 1:length(Ticks)){
  pAB[[i]] = HaploFreqs[[i]]/N 
  
  
  
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


# Calculating D
D = pApB
names(D) = Ticks

for(i in 1:length(Ticks)){
  D[[i]] = pAB[[i]]-pApB[[i]]
}


# Calculating Dmax and then Dprime
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


