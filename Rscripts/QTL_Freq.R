library(dplyr)
colHeader <- function(df) {
  names(df) <- as.character(unlist(df[1, ]))
  df = df[-1, , drop = F]
}

setwd("C:\\Users\\kushv\\genArchSelect\\QTL_prototype\\output_files\\QTL_V2")

file = "QTL_V2_1252395357.csv"

data = read.table(file, sep = ",")
colnames(data) = c("Tick", "Mutation", "Freq")

Ticks = unique(unlist(data[,1]))
SplitData = vector(mode = "list", length = length(Ticks))

for (i in 1:length(Ticks)){
  SplitData[[i]] = data[data[,1]==Ticks[i],]  
}

names(SplitData) = Ticks



# Old Code ------------

# lines = readLines(file, n = -1)
# test = (strsplit(lines, split = ","))
# test = lapply(test, '[', -1)
# test2 = lapply(test, strsplit, split = ": ", fixed = TRUE)
# test2 = lapply(test2, as.data.frame)
# test2 = lapply(test2, colHeader)
# 
# final2 = data.table::rbindlist(test2, fill = T)
# 
# final = do.call(bind_rows, test2)
# final <- Filter(function(x)!all(is.na(x)), final)

#data.table::fwrite(final, file = "test1.csv")


# while (TRUE) {
#   line = readLines(file, n = 1)
#   if (length(line) == 0) {
#     break
#   }
#   row = unlist(strsplit(line, split = ","))
#   test = strsplit(row[-1], split = ": ")
#   test = do.call(cbind, test)
#   test = data.frame(test)
#   test = colHeader(test)
#   data[[i]] = test
#   i = i + 1
# }
# data = do.call(bind_rows, data)
# data <- Filter(function(x)!all(is.na(x)), data)
# 
# close(file)
# 

