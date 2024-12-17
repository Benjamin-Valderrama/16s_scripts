library(dada2)
library(phyloseq)
library(ggplot2)
library(microbiome)
library(tidyr)


################################################
####   ARGUMENTS GIVEN BY THE MAIN SCRIPT   ####
study_folder <- commandArgs(trailingOnly = TRUE)[1]

################################################

pathF <- paste0(study_folder, "/00.rawdata/F")
pathR <- paste0(study_folder, "/00.rawdata/R")


filtpathF <- file.path(pathF, "filtered") # Filtered forward files will go into the pathF/filtered/ subdirectory
filtpathR <- file.path(pathR, "filtered") # ...

fastqFs <- sort(list.files(pathF, pattern="_1.fastq.gz"))
fastqRs <- sort(list.files(pathR, pattern="_2.fastq.gz"))

if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")




# FILTER AND TRIMMING (NEEDS CUSTOMIZATION)
print("FILTER AND TRIMMING ...")
out <- filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              trimLeft=30,
#	      trimRight=30,
              truncLen = 150,
	      maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)

filtFs <- list.files(filtpathF, pattern=".fastq.gz", full.names = TRUE) # HERE
filtRs <- list.files(filtpathR, pattern=".fastq.gz", full.names = TRUE) # HERE

sample.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
sample.namesR <- sapply(strsplit(basename(filtRs), "_"), `[`, 1)

if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")

names(filtFs) <- sample.names
names(filtRs) <- sample.names




# LEARN ERROR RATES 
print("LEARN ERROR RATES ...")
set.seed(01021997)

# forward 
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)
# reverse
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)



# DEREPLICATION
print("DEREPLICATION ...")

# sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

dadasF <- vector("list", length(sample.names))
names(dadasF) <- sample.names

dadasR <- vector("list", length(sample.names))
names(dadasR) <- sample.names


# next step it takes a lot of time
for(sam in sample.names){
  cat("Processing:", sam, "\n")

  derepF <- derepFastq(filtFs[[sam]])
  ddF <- dada(derepF, err=errF, multithread=TRUE, verbose=FALSE)
  dadasF[[sam]] <- ddF
  
  derepR <- derepFastq(filtRs[[sam]])
  ddR <- dada(derepR, err=errR, multithread=TRUE, verbose=FALSE)
  dadasR[[sam]] <- ddR
  
  merger <- mergePairs(ddF, derepF, ddR, derepR)
  mergers[[sam]] <- merger
}
rm(derepF); rm(derepR)



# CONSTRUCT SEQUENCE TABLE AND REMOVE CHIMERAS
seqtab <- makeSequenceTable(mergers)

print("REMOVE CHIMERAS ...")
# check the error here
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE)



# TRACK READS THROUGH THE PIPELINE
print("TRACK READS THROUGH THE PIPELINE ...")

getN <- function(x){sum(getUniques(x))}
track <- cbind(out, sapply(mergers, getN), sapply(dadasF, getN), sapply(dadasR, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR","merged", "nochim")
rownames(track) <- sample.names

write.table(x = track, file = paste0(study_folder, "/01.dada2_o/track_reads_through_pipeline.tsv"),
	    sep = "\t", 
            row.names = TRUE, col.names = TRUE)



# ASSIGN TAXONOMY
print("ASSIGN TAXONOMY ...")

tax <- assignTaxonomy(seqtab.nochim, "/data/databases/SILVA/silva_nr_v138_train_set.fa.gz", multithread=TRUE)



# WRITE DADA2 OUTPUTS
print("WRITE DADA2 OUTPUTS ...")

saveRDS(seqtab.nochim, paste0(study_folder, "/01.dada2_o/seqtab.rds"))
saveRDS(tax, paste0(study_folder, "/01.dada2_o/tax.rds"))

# I can't stop feeling like there are part of this code that are suboptimal... however, I decided to move forward and use it
# might be worth to double check and look for more efficient ways to achieve the same.


# get the names of the samples into a data frame for the phyloseq object
samples.out <- rownames(seqtab)
samdf <- data.frame(ID = samples.out)
row.names(samdf) = samples.out

# phyloseq object is created
ps <- phyloseq(otu_table(seqtab, taxa_are_rows=FALSE), sample_data(samdf), tax_table(tax))

# the count table is created. It has already the taxonomical annotation
count_table_tax  = t(rbind(ps@otu_table, t(ps@tax_table)))

# taxonomic annotation is removed (WHY?)
count_table_one_word = data.frame(count_table_tax[,1:(nrow(ps@otu_table))])


# taxonomic annotation is added to the count table again, this time in just one line
count_table_one_word$taxonomy = paste0(count_table_tax[,"Kingdom"], ";" , 
                                       count_table_tax[,"Phylum" ], ";" , 
                                       count_table_tax[,"Class"  ], ";" , 
                                       count_table_tax[,"Order"  ], ";" , 
                                       count_table_tax[,"Family" ], ";" , 
                                       count_table_tax[,"Genus"  ], ";")

write.csv(x = count_table_one_word, 
	  quote = F, 
	  file = paste0(study_folder, "/01.dada2_o/count_table.csv"))


# IT IS POSSIBLE TO AGGREGATE. I DON'T SEE THE POINT OF THIS CUZ GENUS IS ALREADY THE LOWEST TAXONOMIC LEVEL AVAILABLE.
# Genus_count_table

#agg <- aggregate_taxa(ps, level = "Genus",  verbose = TRUE) #Change me of you want a different level for your count table, Genus is reccomended.
#colap <- unite(as.data.frame(agg@tax_table@.Data), newCol, -unique) 
#genus_table <- agg@otu_table@.Data
#row.names(genus_table) <- colap$newCol

#write.csv(x = genus_table, quote = F, file = "../../data/rawdata/thomaz/south_africa/larissa/genus_table_from_dada2.csv")



# GENERATING THE OUTPUTS FOR PICRUSt2
print("WRITING PICRUSt2's INPUT FILES")

pip = t(ps@otu_table@.Data)

pip_otu_table = pip
row.names(pip_otu_table) = c(paste("otu", 1:nrow(pip_otu_table), sep = ""))

OTU_table = data.frame("OTU"=rownames(pip_otu_table),pip_otu_table)
OTU_table$taxonomy = paste(unlist(t(ps@tax_table@.Data)[5,]), unlist(t(ps@tax_table@.Data)[6,]), sep = "|")
OTU_table = OTU_table[apply(OTU_table[,-c(1, ncol(OTU_table))] == 0, 1, sum) <= ((ncol(OTU_table) -2)  * 0.9), ] #Remove features with prevalence < 10%


write.table(OTU_table, 
            file = paste0(study_folder, "/02.picrust2/input/otu_table_with_taxonomy.csv"), 
	    row.names=FALSE, quote = F, sep = ",")

write.table(OTU_table[,-ncol(OTU_table)], 
            file = paste0(study_folder, "/02.picrust2/input/otu_table.tsv"), 
	    row.names=FALSE, quote = F, sep = "\t")


repseqs <- c(paste(">otu", 1:nrow(pip_otu_table), "\n", row.names(pip), "\n", sep = ""))
write.table(unname(c(paste(">otu", 
                           1:nrow(pip_otu_table), 
                           "\n", row.names(pip), 
                           sep = ""))), 
            file = paste0(study_folder, "/02.picrust2/input/otu_sequences.fasta"), 
	    quote = F, row.names = F, col.names = F)



# SAVING THE IMAGE OF THE R ENVIRONMENT
save.image(file = paste0(study_folder, "/01.dada2_o/environment.RData"))

print("DADA2 FINISHED")
