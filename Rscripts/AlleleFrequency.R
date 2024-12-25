library(ggplot2)
library(dplyr)

AlleleFile = "Insert path to allele frequency information here"

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

# CombinedFrequencies is the final dataset but combined into one data.table
CombinedFrequencies = data.table::rbindlist(AlleleFreqs)


# Optionally extract the effect size
AlleleFreqs = lapply(AlleleFreqs, function(x){
  x$SelCoeff = as.numeric(gsub(".*:(-?[0-9.]+)>.*", "\\1", x$Mutation))
  x$WeightedContribution = x$SelCoeff * x$Freq
  return(x)
})

# Optionally calculate the weighted mean of effect size and frequency, to 
# estimate the additive genotype effect
Expected = unlist(lapply(AlleleFreqs, function(x)2*mean(x$WeightedContribution)))
ggplot(as.data.frame(Expected), aes(x = Ticks, y = Expected)) + geom_line()


# Plotting number of QTLs
ggplot(data.frame(Ticks = Ticks, Sites = unlist(lapply(AlleleFreqs, nrow))), aes(Ticks, Sites))+
  geom_line()+
  labs(title = FitFunc, y = "Number of QTLs") + 
  theme_bw()

# Plotting the number of fixed QTLs
FixedQTLs = data.frame(Ticks = Ticks, FixedQTLs = as.numeric(lapply(AlleleFreqs, function(x){
  return(sum(x$Freq == 1))
})))

ggplot(FixedQTLs, aes(x = Ticks, y = FixedQTLs)) +
  geom_line()+
  theme_bw()


# Plotting allele Frequencies

CombinedFrequencies %>%
  ggplot( aes(x = Tick,
              y = Freq,
              group = Mutation)) +
  geom_line() +
  theme_bw()


# Plotting changes in allele frequncy per generation


## Converting to wide type data
Wide_Alleledata = pivot_wider(CombinedFrequencies, names_from = Tick, values_from = Freq)

## Remove any fixed QTLs as they would not change in frequency and would thus
## confound the data
if(any(Wide_Alleledata[,2] == 1, na.rm = T)){
  Wide_Alleledata = Wide_Alleledata[-which(Wide_Alleledata[,2] == 1),]
}


## Convert to dataframe
deltaPvalues = data.frame(Mutation = Wide_Alleledata$Mutation)

## In case the population size is 1000, a larger generation width was used to make the
## graph more accessible
if(grepl("1k", FitFunc)){
  generationwidth=10
  col_index = 2
  
  for(i in seq(2+generationwidth, ncol(Wide_Alleledata), generationwidth)){
    deltaPvalues[,col_index] = Wide_Alleledata[i]-Wide_Alleledata[i-generationwidth]
    col_index = col_index + 1
  }
} else{ ## Otherwise, the generation width was 1 time point
  for(i in 3:ncol(Wide_Alleledata)){
    deltaPvalues[,i-1] = Wide_Alleledata[i]-Wide_Alleledata[i-1]
  }
}


## Convert to long type data
deltaPvalues_long = pivot_longer(deltaPvalues, !Mutation, names_to = "ticks", values_to = "values")


## Coercing ticks to numeric
deltaPvalues_long$ticks = as.numeric(deltaPvalues_long$ticks)

## Sorting by generations
deltaPvalues_long = deltaPvalues_long %>%
  arrange(ticks)

Ticks = unique(deltaPvalues_long$ticks)
BreakTicks = Ticks[seq(1, length(Ticks)-1, length.out = 8)]

## Plotting the data
## Note that scale_x_discrete may have to be tweaked to obtain an x-axis that is accessible
ggplot(deltaPvalues_long, aes(x = factor(ticks), 
                                  y = values)) +
  geom_boxplot(na.rm = T,
               outlier.size = 0.6)+
  labs(x = "Generations",
       y = paste("\u394", "Allele Frequency")) +
  scale_x_discrete(breaks = BreakTicks)+
  theme_bw() + 
  theme(axis.text.x  = element_text(angle = 90,
                                    hjust = 1),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
