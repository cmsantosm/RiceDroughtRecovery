library(DESeq2)
library(biobroom)
library(tidyverse)

otu <- readRDS("../Data/otu_pers.RDS")
map <- readRDS("../Data/drought_map.RDS")

map <- filter(map, Compartment == "ES") 

map <- map %>% 
  mutate(Group = paste(Treatment2,Time, sep = "."))

otu <- otu[, colnames(otu) %in% map$SampleID]
otu <- otu[, match(map$SampleID, colnames(otu))]
otu <- otu[rowSums(otu) > 0, ]

dds<- DESeqDataSetFromMatrix(countData = otu,
                             colData = map,
                             design = ~ Group)

dds <- DESeq(dds)

contrasts <- vector(mode = "list")

for(i in 1:13) {
  contrasts[[paste("WC_TRN", i, sep = ".")]] <- c("Group", paste("WC", i, sep = "."), paste("WC_TRN", i, sep = "."))
  contrasts[[paste("D1", i, sep = ".")]] <- c("Group", paste("WC", i, sep = "."), paste("D1", i, sep = "."))
  contrasts[[paste("D2", i, sep = ".")]] <- c("Group", paste("WC", i, sep = "."), paste("D2", i, sep = "."))
  contrasts[[paste("D3", i, sep = ".")]] <- c("Group", paste("WC", i, sep = "."), paste("D3", i, sep = "."))
}

results <- vector(mode = "list")
shrinkFC <- vector(mode = "list")

for(i in seq_along(contrasts)) {
  results[[names(contrasts)[i]]] <- tidy(results(dds, contrast = contrasts[[i]])) %>% mutate(Day = i)
  shrinkFC[[names(contrasts)[i]]] <- tidy(lfcShrink(dds, contrast = contrasts[[i]])) %>% mutate(Day = i)
  print(contrasts[[i]])
}

results.df <- plyr::ldply(results, function(x) x)
names(results.df)[1] <- "Contrast"
names(results.df)[2] <- "OTU_ID"

shrinkFC.df <- plyr::ldply(shrinkFC, function(x) x)
names(shrinkFC.df)[1] <- "Contrast"
names(shrinkFC.df)[2] <- "OTU_ID"

saveRDS(results.df, "Data/es_dab_bal.RDS")
saveRDS(shrinkFC.df, "Data/es_shlfc.RDS")
saveRDS(dds, "Data/es_dds.RDS")

