library(xlsx)
library(dplyr)
library(plyr)

##Parsing functional annotations of protein sequences 
##by KEGG (Kyoto Encyclopedia of Genes and Genomes)

setwd("working_directory")

#read in the references and format
  ref <- read.csv("KO_Orthology_ko00001.txt",sep = "\t", header = FALSE)
  
  #remove kegg identifiers of pathways
  ref$V2 <- sub(".*? ", "", ref$V2)
  ref$V3 <- sub(".*? ", "", ref$V3)
  
  #make column names for refs
  ref_name <- c("Cat","Metabolism","Pathway","KO","Gene","Function")
  names(ref) <- ref_name
  ref <- ref[ref$Cat == "Metabolism", ] 
  
  ref <- ref[-1]

#for loop to run through each of the results
  bins <- read.csv("true_all.txt",sep = "\t", header = FALSE)
  new <- ref
  GB <- read.csv("GB_names.txt",sep = "\t", header = FALSE)
  
    for(i in 1:nrow(bins)){
      
      #import results
      kegg <- read.csv(paste0(bins[i,],".faa.tsv"),sep = "\t",header = FALSE,na.strings=c("","NA"),colClasses = c("character","character"))
      if(!( bins[i,] %in% GB$V1)){
        kegg <- kegg[2:3]
      }
      
      #remove proteins with no KO hits
      kegg <- na.omit(kegg)
      names(kegg) <- c(as.character(bins[i,]),"KO")
      
      new <- merge.data.frame(new,kegg, by = "KO", all=TRUE)
      
      #only keep unique values separated by the ;
      new <- new  %>%
        group_by(KO,Metabolism,Pathway) %>%
        summarise_all(funs(paste(unique(.), collapse=";")))
    }
  
  #make all NA values blanks
  new[is.na(new)] = ""
  
  #export data
  write.table(new, file='KEGG_results_all.tsv', quote=FALSE, sep='\t', row.names = FALSE)
