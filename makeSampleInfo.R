library(argparse)
library(RMySQL)
library(plyr)
library(dplyr)
library(magrittr)
library(Biostrings)

setArguments <- function(){
  parser <- ArgumentParser(description = "Generate sampleInfo.tsv from databased information.")
  parser$add_argument("runDate", nargs = "?", default = "2015-09-29", help = "Run date in YYYY-MM-DD format.")
  parser$add_argument("-d", "--database", default = "hiv_specimen.database", help = "Group to use for specimen data.")
  parser$add_argument("-g", "--refGenome", default = "hg18", help = "Reference genome to use in intSiteCaller.")
  
  arguments <- parser$parse_args()
  arguments
}

arguments <- setArguments()

runDate <- arguments$runDate
refGenome <- arguments$refGenome
vectorSeq <- arguments$vectorSeq
specimenDatabase <- arguments$database

#Connect and query information from specimen database
junk <- sapply(dbListConnections(MySQL()), dbDisconnect) 
dbConn <- dbConnect(MySQL(), group = specimenDatabase)  
stopifnot(dbGetQuery(dbConn, "SELECT 1")==1)

query_selection_hivsam <- "SELECT parentAlias,childAlias,linkerNum,bcNum,uniqRegion FROM hivsam"
query_selection_contsam <- "SELECT parentAlias,childAlias,linkerNum,bcNum,uniqRegion FROM contsam"
query_condition_1 <- paste0("WHERE dateSeqd LIKE '", runDate, "'")
hivsam <- dbGetQuery(dbConn, paste(query_selection_hivsam, query_condition_1))
contsam <- dbGetQuery(dbConn, paste(query_selection_contsam, query_condition_1))

query_selection_hivsp <- "SELECT parentAlias,trial,gender FROM hivsp"
parentAliasString <- paste(unique(hivsam$parentAlias), collapse = "' OR parentAlias = '")
query_condition_2 <- paste0("WHERE parentAlias = '", parentAliasString, "'")
hivsp <- dbGetQuery(dbConn, paste(query_selection_hivsp, query_condition_2))

query_selection_trialseqinfo <- "SELECT trial,uniqRegion,primer,ltrBit,largeLTRFrag,vectorSeq FROM trialseqinfo"
query_condition_3 <- paste0("WHERE trial = '", unique(hivsp$trial), 
                           "' AND uniqRegion = '", unique(hivsam$uniqRegion), "'")
trialseqinfo <- dbGetQuery(dbConn, paste(query_selection_trialseqinfo, query_condition_3))

query_selection_linkerseqs <- "SELECT linkerNum,sequence FROM linkerseqs"
linkerNumString <- paste(unique(c(hivsam$linkerNum, contsam$linkerNum)), collapse = "' OR linkerNum = '")
query_condition_4 <- paste0("WHERE linkerNum = '", linkerNumString, "'")
linkerseqs <- dbGetQuery(dbConn, paste(query_selection_linkerseqs, query_condition_4))

query_selection_bcseqs <- "SELECT bcNum,revComp FROM bcseqs"
bcNumString <- paste(unique(c(hivsam$bcNum, contsam$bcNum)), collapse = "' OR bcNum = '")
query_condition_5 <- paste0("WHERE bcNum = '", bcNumString, "'")
bcseqs <- dbGetQuery(dbConn, paste(query_selection_bcseqs, query_condition_5))

dbDisconnect(dbConn)

#Join data together to generate the sampleInfo.tsv and vector.fasta file
contGender <- "m"
contTrial <- unique(hivsp$trial)
contsp <- data.frame(
  "parentAlias" = unique(contsam$parentAlias),
  "trial" = contTrial,
  "gender" = contGender)
specimens <- bind_rows(hivsp, contsp)

vector <- unlist(strsplit(trialseqinfo$vectorSeq, split = ":"))
writeXStringSet(DNAStringSet(vector[2]), filepath = paste0(vector[1], ".fasta"), format = "fasta")

sampleInfo <- bind_rows(hivsam, contsam) %>%
  left_join(., specimens, by = "parentAlias") %>%
  left_join(., linkerseqs, by = "linkerNum") %>%
  left_join(., bcseqs, by = "bcNum") %>%
  mutate(., primer = trialseqinfo$primer) %>%
  mutate(., ltrBit = trialseqinfo$ltrBit) %>%
  mutate(., largeLTRFrag = trialseqinfo$largeLTRFrag) %>%
  mutate(., vectorSeq = paste0(vector[1], ".fasta")) %>%
  select(., childAlias, sequence, revComp, gender, primer, ltrBit, uniqRegion, largeLTRFrag, vectorSeq)
names(sampleInfo) <- c("alias", "linkerSequence", "bcSeq", "gender", 
                    "primer", "ltrBit", "uniqRegion", "largeLTRFrag", "vectorSeq")
write.table(sampleInfo, file = "sampleInfo.tsv", sep = "\t", row.names = FALSE, quote = FALSE)


