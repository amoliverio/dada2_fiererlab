# -----------------------------------------------------------------------------------#
# dada2 tutorial with miseq dataset for Fierer Lab 
# -----------------------------------------------------------------------------------#

# This version runs the dada2 workflow for Big Data (paired-end) from Rstudio on our server.
    # More information here: https://benjjneb.github.io/dada2/bigdata_paired.html
    # I suggest opening the dada2 tutorial online to understand more about each step

# Note there is a slightly different pipeline for ITS and Non-"Big data" that takes a
    # closer look at each step here:     
    ## link to tutorial: https://benjjneb.github.io/dada2/tutorial.html

# -----------------------------------------------------------------------------------#
# Set up - Steps before starting pipeline
# -----------------------------------------------------------------------------------#

# If you are running it through fierer lab "microbe" server:
##Logging in: To use microbe, open a terminal window, and type and hit return:
    # ssh <your microbe user name>@microbe.colorado.edu

# Login to RStudio on the server (NOTE: you can also run from your own computer, slower!)
    # ssh -N -f -L localhost:####:localhost:8787 <yourUSERNAME>@microbe.colorado.edu

# -----------------------------------------------------------------------------------#
# Set up - You are logged in to Rstudio on server # (Or have it open on your computer)
# -----------------------------------------------------------------------------------#

# Install DADA2 & needed binaries
# WARNING: This may take a long time, so only do 
# this if these packages are not already installed!

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.8")
source("https://bioconductor.org/biocLite.R")
biocLite("ShortRead")

install.packages("dplyr")

# Load DADA2 and required packages
library(dada2); packageVersion("dada2")
library(ShortRead)
library(dplyr)




# pathway to idemp (demultiplexing tool)
idemp <- "/main/data/hollandh/bin/idemp" # CHANGE ME to the demultiplex idemp path;
system2(idemp) # Run shell commands from R

# Set path to shared data folder and contents
data.fp <- "/data/shared/2019_02_20_MicrMethods_tutorial"  ## CHANGE ME to the directory containing the fastq files.
list.files(data.fp) # Should list all files in shared folder
barcode.fp <- file.path(data.fp, "barcode_demultiplex_short.txt") # format is .txt file: barcode </t> sampleID
map.fp <- file.path(data.fp, "Molecular_Methods_18_515fBC_16S_Mapping_File_SHORT_vFinal_Fierer_10252018.txt")
I1.fp <- file.path(data.fp, "Undetermined_S0_L001_I1_001.fastq.gz") 
R1.fp <- file.path(data.fp, "Undetermined_S0_L001_R1_001.fastq.gz") 
R2.fp <- file.path(data.fp, "Undetermined_S0_L001_R2_001.fastq.gz") 

# Set up file paths in YOUR directory where you want data; you do not need to create these directories
project.fp <- "/data/hollandh/MicroMethods_dada2_tutorial"

# Set up names of sub directories to stay organized
preprocess.fp <- file.path(project.fp, "01_preprocess")
    demultiplex.fp <- file.path(preprocess.fp, "demultiplexed")
    filtN.fp <- file.path(preprocess.fp, "filtN")
    trimmed.fp <- file.path(preprocess.fp, "trimmed")
filter.fp <- file.path(project.fp, "02_filter") 
table.fp <- file.path(project.fp, "03_tabletax") 

# -----------------------------------------------------------------------------------#
# Pre-processing data for dada2 - demultiplex, cutadapt 
# -----------------------------------------------------------------------------------#

 # Call demultiplex script
flags <- paste("-b", barcode.fp, "-I1", I1.fp, "-R1", R1.fp, "-R2", R2.fp, "-o", demultiplex.fp) 
system2(idemp, args = flags)

# Look at output of demultiplex
list.files(demultiplex.fp)

# Move unassignable reads (so they are not included in downstream analyses)
unassigned_1 <- paste0("mv", " ", demultiplex.fp, "/Undetermined_S0_L001_R1_001.fastq.gz_unsigned.fastq.gz", " ", demultiplex.fp, "/Unassigned_reads1.fastq.gz")
unassigned_2 <- paste0("mv", " ", demultiplex.fp, "/Undetermined_S0_L001_R2_001.fastq.gz_unsigned.fastq.gz", " ", demultiplex.fp, "/Unassigned_reads2.fastq.gz")
system(unassigned_1)
system(unassigned_2)

# Rename files - gsub to get names in order!
R1_names <- gsub(paste0(demultiplex.fp, "/Undetermined_S0_L001_R1_001.fastq.gz_"), "", list.files(demultiplex.fp, pattern="R1", full.names = TRUE))
file.rename(list.files(demultiplex.fp, pattern="R1", full.names = TRUE), paste0(demultiplex.fp, "/R1_", R1_names))
R2_names <- gsub(paste0(demultiplex.fp, "/Undetermined_S0_L001_R2_001.fastq.gz_"), "", list.files(demultiplex.fp, pattern="R2", full.names = TRUE))
file.rename(list.files(demultiplex.fp, pattern="R2", full.names = TRUE), paste0(demultiplex.fp, "/R2_", R2_names))

# Get full paths for all files for cutadapt
# Forward and reverse fastq filenames have format: 
fnFs <- sort(list.files(demultiplex.fp, pattern="R1_", full.names = TRUE))
fnRs <- sort(list.files(demultiplex.fp, pattern="R2_", full.names = TRUE))

# Use cutadapt to trim reads
FWD <- "GTGYCAGCMGCCGCGGTAA"  ## CHANGE ME to your forward primer sequence # this is 515f # Not Rev-Comp
REV <- "GGACTACNVGGGTWTCTAAT"  ## CHANGE ME # this is 806Br # Not Rev-Comp

# Verify presence and orientation of primers
allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
                 RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# Pre-filter to remove Ns
fnFs.filtN <- file.path(preprocess.fp, "filtN", basename(fnFs)) # Put N-filtered files in filtN/ subdirectory
fnRs.filtN <- file.path(preprocess.fp, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

# Count how many time primers appear
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}

# Look at primer detection for first file
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

### Remove primers with cutadapt ###

# Load cutadapt
cutadapt <- "/usr/local/Python27/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Run shell commands from R

# Trim primers
if(!dir.exists(trimmed.fp)) dir.create(trimmed.fp)
fnFs.cut <- file.path(trimmed.fp, basename(fnFs))
fnRs.cut <- file.path(trimmed.fp, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 

# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# Run Cutadapt
for(i in seq_along(fnFs)) {
    system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                               "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                               fnFs.filtN[i], fnRs.filtN[i])) # input files
}

# As a sanity check, we will count the presence of primers in the first cutadapt-ed sample:
    ## should all be zero!
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# -----------------------------------------------------------------------------------#
###  Now start DADA2 pipeline ###
# -----------------------------------------------------------------------------------#

library(dada2); packageVersion("dada2")

# Put filtered reads into separate sub-directories for big data workflow
dir.create(filter.fp)
    subF.fp <- file.path(filter.fp, "preprocessed_F") 
    subR.fp <- file.path(filter.fp, "preprocessed_R") 
dir.create(subF.fp)
dir.create(subR.fp)

### copy R1 and R2 from trimmed to sub-directories
fnFs.Q <- file.path(subF.fp,  basename(fnFs)) 
fnRs.Q <- file.path(subR.fp,  basename(fnRs))
file.copy(from = fnFs.cut, to = fnFs.Q)
file.copy(from = fnRs.cut, to = fnRs.Q)

# File parsing
filtpathF <- file.path(subF.fp, "filtered") # Filtered forward files go into the pathF/filtered/ subdirectory
filtpathR <- file.path(subR.fp, "filtered") # ...
fastqFs <- sort(list.files(subF.fp, pattern="fastq.gz"))
fastqRs <- sort(list.files(subR.fp, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

# -----------------------------------------------------------------------------------#
# 1. FILTER AND TRIM FOR QUALITY
# -----------------------------------------------------------------------------------#

# Filtering: THESE PARAMETERS ARENT OPTIMAL FOR ALL DATASETS
filterAndTrim(fwd=file.path(subF.fp, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(subR.fp, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(150,150), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)

# -----------------------------------------------------------------------------------#
# 2. INFER sequence variants
# -----------------------------------------------------------------------------------#

# File parsing
filtFs <- list.files(filtpathF, pattern="fastq.gz", full.names = TRUE)
filtRs <- list.files(filtpathR, pattern="fastq.gz", full.names = TRUE)

# Sample names in order
sample.names <- substring(basename(filtFs), regexpr("_", basename(filtFs)) + 1) # doesn't drop fastq.gz
sample.names <- gsub(".fastq.gz", "", sample.names)
sample.namesR <- substring(basename(filtRs), regexpr("_", basename(filtRs)) + 1) # doesn't drop fastq.gz
sample.namesR <- gsub(".fastq.gz", "", sample.namesR)

# Double check
if(!identical(sample.names, sample.namesR)) stop("Forward and reverse files do not match.")
names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)

# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

# Sample inference and merger of paired-end reads
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

for(sam in sample.names) {
    cat("Processing:", sam, "\n")
    derepF <- derepFastq(filtFs[[sam]])
    ddF <- dada(derepF, err=errF, multithread=TRUE)
    derepR <- derepFastq(filtRs[[sam]])
    ddR <- dada(derepR, err=errR, multithread=TRUE)
    merger <- mergePairs(ddF, derepF, ddR, derepR)
    mergers[[sam]] <- merger
}

rm(derepF); rm(derepR)

# Construct sequence table
seqtab <- makeSequenceTable(mergers)

# save table
dir.create(table.fp)
saveRDS(seqtab, paste0(table.fp, "/seqtab.rds")) # CHANGE ME to where you want sequence table saved

# -----------------------------------------------------------------------------------#
# 3. Chimeras and Taxonomy Assignment
# -----------------------------------------------------------------------------------#

# Read in RDS 
st.all <- readRDS(paste0(table.fp, "/seqtab.rds"))

# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)

# Assign taxonomy
tax <- assignTaxonomy(seqtab, "/data/oliverioa/dada2_tutorials/silva_nr_v132_train_set.fa.gz",
                      multithread=TRUE)

# Write to disk
saveRDS(seqtab, paste0(table.fp, "/seqtab_final.rds")) # CHANGE ME to where you want sequence table saved
saveRDS(tax, paste0(table.fp, "/tax_final.rds")) # CHANGE ME ...

# -----------------------------------------------------------------------------------#
# 4. Optional - Format to obtain ESV IDs and repset, and Input for mctoolsr
# -----------------------------------------------------------------------------------#

# Flip table
seqtab.t <- as.data.frame(t(seqtab))

# Pull out ESV repset
rep_set_ESVs <- as.data.frame(rownames(seqtab.t))
rep_set_ESVs <- mutate(rep_set_ESVs, ESV_ID = 1:n())
rep_set_ESVs$ESV_ID <- sub("^", "ESV_", rep_set_ESVs$ESV_ID)
rep_set_ESVs$ESV <- rep_set_ESVs$`rownames(seqtab.t)` 
rep_set_ESVs$`rownames(seqtab.t)` <- NULL

# Add ESV numbers to table
rownames(seqtab.t) <- rep_set_ESVs$ESV_ID

# Add ESV numbers to taxonomy
taxonomy <- as.data.frame(tax)
taxonomy$ESV <- as.factor(rownames(taxonomy))
taxonomy <- merge(rep_set_ESVs, taxonomy, by = "ESV")
rownames(taxonomy) <- taxonomy$ESV_ID

# Mctoolsr format for downstream analyses

# Read in your map file w/sample IDs as rownames
map <- read.delim(map.fp, row.names = 1, as.is = TRUE) 

input <- list()
input$data_loaded <- seqtab.t
input$taxonomy_loaded <- taxonomy
input$map_loaded <- map

# Number of reads per sample
sort(colSums(input$data_loaded))

# Now the input should be compatible with MCTOOLSR, to save:
save(input, rep_set_ESVs, file = paste0(table.fp, "/input_tutorial.rda"))

# Also export files as .txt
write.table(seqtab.t, file = paste0(table.fp, "/seqtab_final.txt"))
write.table(tax, file = paste0(table.fp, "/tax_final.txt"))

# Note: you will need to think about the following in downstream applications:
    #1 Remove mitochondrial and chloroplast sequences, for example with mctoolsr:
            # input_filt <- filter_taxa_from_input(input, taxa_to_remove = "Chloroplast")
            # input_filt <- filter_taxa_from_input(input_filt, taxa_to_remove = "Mitochondria")
    #2 Normalize or rarefy your ESV table

# You can now transfer over the output files onto your local computer
    # .rda files are loaded into R like this: load(file = "filepath/input_tutorial.rda")
