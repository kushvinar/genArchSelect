# Loading in the data 
library(tidyverse)
library(ggplot2)

AlleleFile = files[3]

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



## Converting to wide type data
Wide_Alleledata = pivot_wider(CombinedFrequencies, names_from = Tick, values_from = Freq)

## Remove any fixed QTLs as they would not change in frequency and would thus
## confound the data
if(any(Wide_Alleledata[,2] == 1, na.rm = T)){
  Wide_Alleledata = Wide_Alleledata[-which(Wide_Alleledata[,2] == 1),]
}


## Convert to dataframe
deltaPvalues = data.frame(Mutation = Wide_Alleledata$Mutation)

for(i in 3:ncol(Wide_Alleledata)){
  deltaPvalues[,i-1] = Wide_Alleledata[i]-Wide_Alleledata[i-1]
}



## Convert to long type data
deltaPvalues_long = pivot_longer(deltaPvalues, !Mutation, names_to = "ticks", values_to = "values")


## Coercing ticks to numeric
deltaPvalues_long$ticks = as.numeric(deltaPvalues_long$ticks)

# ## Sorting by generations
# deltaPvalues_long = deltaPvalues_long %>%
#   arrange(ticks)



Ticks = unique(deltaPvalues_long$ticks)
BreakTicks = Ticks[seq(1, length(Ticks)-1, length.out = 8)]

# Note: if median is NA, that means the mutation existed for only one generation

finalData <- deltaPvalues_long %>%
  group_by(Mutation) %>%
  summarise(
    max_abs = values[which.max(abs(values))],
  ) %>%
  filter(!is.na(max_abs)) # Ensure no NA values in the result

# Extracting effect sizes
finalData$EffectSize = as.numeric(gsub(".*:(-?[0-9.]+)>.*", "\\1", finalData$Mutation))


# Final Plot
ggplot(data = finalData,
       aes(x = EffectSize,
           y = max_abs))+
  geom_point()+
  labs(x = "Effect size of mutation",
       y = "Largest AlleleFrequency Change",
       title = paste0(AlleleFile, "     ",
                      paste0("slope= ",coef(lm(finalData$max_abs~finalData$EffectSize))[2])))+
  geom_smooth(method = "lm", se = F)+
  theme_bw()
