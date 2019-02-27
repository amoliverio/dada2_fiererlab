#'# dada2 tutorial with MiSeq dataset for Fierer Lab 

#' This version runs the dada2 workflow for Big Data (paired-end) from Rstudio on the microbe server.
#' 
#' * More information here: [https://benjjneb.github.io/dada2/bigdata_paired.html]([https://benjjneb.github.io/dada2/bigdata_paired.html])
#'    
#' We suggest opening the dada2 tutorial online to understand more about each step
#'
#'--- 
#' **NOTE:** there is a slightly different pipeline for ITS and Non-"Big data" that takes a
#'      closer look at each step here: [https://benjjneb.github.io/dada2/tutorial.html](https://benjjneb.github.io/dada2/tutorial.html)
#'---
#' 
#' 
#' 
#' ## Set up (part 1) - Steps before starting pipeline ##
#' 
#' If you are running it through fierer lab "microbe" server:
#' #### Logging in: 
#' 
#' To use the microbe server, open a terminal window, and type and hit return:
#' 
#' ```bash    
#' ssh <your microbe user name>@microbe.colorado.edu 
#' ```
#'
#' #### Login to RStudio on the server 
#' 
#' (NOTE: you can also run from your own computer, but it will be slower!)
#' 
#'    1. Open a your web browser and start a new empty tab
#'    2. type `microbe.colorado.edu:8787` in the address bar
#'    3. use your server login credentials to log into rstudio server
#'
#' ## Set up (part 2) - You are logged in to Rstudio on server (or have it open on your computer) ##
#' 
#' Install DADA2 & other necessary packages

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("dada2", version = "3.8")

source("https://bioconductor.org/biocLite.R")
biocLite("ShortRead")

install.packages("dplyr")

#'---
#' **WARNING:** This installation may take a long time, so only do this if these 
#' packages are not already installed!
#'---
#'

#' Load DADA2 and required packages

library(dada2); packageVersion("dada2")
library(ShortRead)
library(dplyr)

#' Once the packages are installed, you can check to make sure the auxillary
#' software is working and set up some of the variables that you will need 
#' along the way.

# set up pathway to idemp (demultiplexing tool) and test
idemp <- "/main/data/hollandh/bin/idemp" # CHANGE ME to the demultiplex idemp path;
system2(idemp) # Check that idemp is in your path and you can run shell commands from R

# Set up pathway to cutadapt (primer trimming tool) and test
cutadapt <- "/usr/local/Python27/bin/cutadapt" # CHANGE ME to the cutadapt path on your machine
system2(cutadapt, args = "--version") # Check if it's working by running shell command from R


# Set path to shared data folder and contents
data.fp <- "/data/shared/2019_02_20_MicrMethods_tutorial"  ## CHANGE ME to the directory containing the fastq files.
list.files(data.fp) # Should list all files in shared folder
barcode.fp <- file.path(data.fp, "barcode_demultiplex_short.txt") # format is .txt file: barcode </t> sampleID
map.fp <- file.path(data.fp, "Molecular_Methods_18_515fBC_16S_Mapping_File_SHORT_vFinal_Fierer_10252018.txt") # CHANGE ME to the name of your mapping file
I1.fp <- file.path(data.fp, "Undetermined_S0_L001_I1_001.fastq.gz") 
R1.fp <- file.path(data.fp, "Undetermined_S0_L001_R1_001.fastq.gz") 
R2.fp <- file.path(data.fp, "Undetermined_S0_L001_R2_001.fastq.gz") 

#' Set up file paths in YOUR directory where you want data; 
#' you do not need to create the subdirectories but they are nice to have
#' for organizational purposes. 

project.fp <- "/data/hollandh/MicroMethods_dada2_tutorial" # CHANGE ME to the directory that you want your project to be in

# Set up names of sub directories to stay organized
preprocess.fp <- file.path(project.fp, "01_preprocess")
    demultiplex.fp <- file.path(preprocess.fp, "demultiplexed")
    filtN.fp <- file.path(preprocess.fp, "filtN")
    trimmed.fp <- file.path(preprocess.fp, "trimmed")
filter.fp <- file.path(project.fp, "02_filter") 
table.fp <- file.path(project.fp, "03_tabletax") 

#' ## Pre-processing data for dada2 - demultiplex, remove sequences with Ns, cutadapt 
#' 
#' #### Call the demultiplexing script
#' Demultiplexing splits your reads out into separate files based on the barcodes associated with each sample. 

flags <- paste("-b", barcode.fp, "-I1", I1.fp, "-R1", R1.fp, "-R2", R2.fp, "-o", demultiplex.fp) 
system2(idemp, args = flags) 

# Look at output of demultiplexing
list.files(demultiplex.fp)

#'---
#' **WARNING:** The demultiplexing step may take a while. If it takes too long
#' you can safely close RStudio on the server and the demultiplexing will run 
#' in the background. You should be able to resume the pipeline after demultiplexing
#' is complete by logging back into RStudio on the server.
#'---
#'

#' #### Clean up the output from idemp
#'
 
# Change names of unassignable reads so they are not included in downstream analyses (remove "R1/2" from names).
unassigned_1 <- paste0("mv", " ", demultiplex.fp, "/Undetermined_S0_L001_R1_001.fastq.gz_unsigned.fastq.gz", " ", demultiplex.fp, "/Unassigned_reads1.fastq.gz")
unassigned_2 <- paste0("mv", " ", demultiplex.fp, "/Undetermined_S0_L001_R2_001.fastq.gz_unsigned.fastq.gz", " ", demultiplex.fp, "/Unassigned_reads2.fastq.gz")
system(unassigned_1)
system(unassigned_2)

# Rename files - gsub to get names in order!
R1_names <- gsub(paste0(demultiplex.fp, "/Undetermined_S0_L001_R1_001.fastq.gz_"), "", list.files(demultiplex.fp, pattern="R1", full.names = TRUE))
file.rename(list.files(demultiplex.fp, pattern="R1", full.names = TRUE), paste0(demultiplex.fp, "/R1_", R1_names))
R2_names <- gsub(paste0(demultiplex.fp, "/Undetermined_S0_L001_R2_001.fastq.gz_"), "", list.files(demultiplex.fp, pattern="R2", full.names = TRUE))
file.rename(list.files(demultiplex.fp, pattern="R2", full.names = TRUE), paste0(demultiplex.fp, "/R2_", R2_names))

# Get full paths for all files and save them for downstream analyses
# Forward and reverse fastq filenames have format: 
fnFs <- sort(list.files(demultiplex.fp, pattern="R1_", full.names = TRUE))
fnRs <- sort(list.files(demultiplex.fp, pattern="R2_", full.names = TRUE))

#' #### Pre-filter to remove sequence reads with Ns
#' Ambiguous bases will make it hard for cutadapt to find short primer sequences in the reads.
#' To solve this problem, we will remove sequences with ambiguous bases (Ns)
#' 
fnFs.filtN <- file.path(preprocess.fp, "filtN", basename(fnFs)) # Name the N-filtered files to put them in filtN/ subdirectory
fnRs.filtN <- file.path(preprocess.fp, "filtN", basename(fnRs))

# filter Ns from reads and put them into the filtN directory
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # CHANGE multithread to FALSE on Windows (here and elsewhere in the program)


#' #### Prepare the primers sequences and custom functions for analyzing the results from cutadapt

# set up the primer sequences to pass along to cutadapt
FWD <- "GTGYCAGCMGCCGCGGTAA"  ## CHANGE ME to your forward primer sequence # this is 515f # Not Rev-Comp
REV <- "GGACTACNVGGGTWTCTAAT"  ## CHANGE ME # this is 806Br # Not Rev-Comp

# Write a function that creates a list of all orientations of the primers
allOrients <- function(primer) {
    # Create all orientations of the input sequence
    require(Biostrings)
    dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
    orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
                 RevComp = reverseComplement(dna))
    return(sapply(orients, toString))  # Convert back to character vector
}

# save the primer orientations to pass to cutadapt
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# Write a function that counts how many time primers appear in a sequence
primerHits <- function(primer, fn) {
    # Counts number of reads in which the primer is found
    nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
    return(sum(nhits > 0))
}

#' Before running cutadapt, we will look at primer detection
#' for the first sample, as a check.
#' There may be some primers here, we will remove them below using cutadapt
#' 
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#' #### Remove primers with cutadapt and assess the output using our custom functions

# create directory to hold the output from cutadapt
if (!dir.exists(trimmed.fp)) dir.create(trimmed.fp)
fnFs.cut <- file.path(trimmed.fp, basename(fnFs))
fnRs.cut <- file.path(trimmed.fp, basename(fnRs))

# Save the reverse complements of the primers to variables
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)

##  Create the cutadapt flags ##
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 

# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 

# Run Cutadapt
for (i in seq_along(fnFs)) {
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

#' # Now start DADA2 pipeline

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
filtpathF <- file.path(subF.fp, "filtered") # Filtered forward files go into the path preprocessed_F/filtered/ subdirectory
filtpathR <- file.path(subR.fp, "filtered") # ...
fastqFs <- sort(list.files(subF.fp, pattern="fastq.gz"))
fastqRs <- sort(list.files(subR.fp, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#' ### 1. FILTER AND TRIM FOR QUALITY
#' 
#' Before chosing sequence variants, we want to trim reads where their quality scores 
#' begin to drop (the `truncLen` and `truncQ` values) and remove any low-quality reads 
#' that are left over after we have finished trimming (the `maxEE` value).

#' |**WARNING:** THESE PARAMETERS ARE NOT OPTIMAL FOR ALL DATASETS. Make sure you determine
#' the trim and filtering parameters for your data. The following settings are 
#' generally appropriate for MiSeq runs.|
#' | --- |
#'

filterAndTrim(fwd=file.path(subF.fp, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(subR.fp, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(150,150), maxEE=2, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)

#' ### 2. INFER sequence variants
#' In this part of the pipeline dada2 will learn to distinguish error from biological 
#' differences using a subset of our data as a training set. After it understands the 
#' error rates, we will reduce the size of the dataset by combining all identical 
#' sequence reads into "unique sequences". Then, using the dereplicated data and 
#' error rates, dada2 will infer the sequence variants (OTUs) in our data. Finally, 
#' we will merge the coresponding forward and reverse reads to create a list of the 
#' fully denoised sequences and create a sequence table from the result.

#' #### Housekeeping - set up and verify the file names for the output
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

#' #### Learn the error rates
set.seed(100) # set seed to ensure that randomized steps are replicatable

# Learn forward error rates
errF <- learnErrors(filtFs, nbases=1e8, multithread=TRUE)

# Learn reverse error rates
errR <- learnErrors(filtRs, nbases=1e8, multithread=TRUE)

#' #### Dereplication, sequence inference, and merging of paired-end reads
# make a list to hold the loop output
mergers <- vector("list", length(sample.names))
names(mergers) <- sample.names

# for each sample, get a list of merged and denoised sequences
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

#' #### Construct sequence table
seqtab <- makeSequenceTable(mergers)

# save table as an r data object file
dir.create(table.fp)
saveRDS(seqtab, paste0(table.fp, "/seqtab.rds")) # CHANGE ME to where you want sequence table saved

#' ### 3. REMOVE Chimeras and ASSIGN Taxonomy
#' Although dada2 has searched for indel errors and subsitutions, there may still be chimeric
#' sequences in our dataset (sequences that are derived from forward and reverse sequences from 
#' two different organisms becoming fused together during PCR and/or sequencing). To identify 
#' chimeras, we will search for rare sequence variants that can be reconstructed by combining
#' left-hand and right-hand segments from two more abundant "parent" sequences.
#' 
#' After removing chimeras, we will use a taxonomy database to train a classifer-algorithm
#' to assign names to our sequence variants. 

# Read in RDS 
st.all <- readRDS(paste0(table.fp, "/seqtab.rds"))

# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)

# Assign taxonomy
tax <- assignTaxonomy(seqtab, "/data/oliverioa/dada2_tutorials/silva_nr_v132_train_set.fa.gz",
                      multithread=TRUE)

# Write results to disk
saveRDS(seqtab, paste0(table.fp, "/seqtab_final.rds"))
saveRDS(tax, paste0(table.fp, "/tax_final.rds"))

#' ### 4. Optional - FORMAT OUTPUT to obtain ESV IDs and repset, and input for mctoolsr
#' For convenience sake, we will now rename our ESVs with numbers, output our 
#' results as a traditional taxa table, and create a matrix with the representative
#' sequences for each ESV. 

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

## Create Mctoolsr format for downstream analyses

# Add ESV numbers to taxonomy
taxonomy <- as.data.frame(tax)
taxonomy$ESV <- as.factor(rownames(taxonomy))
taxonomy <- merge(rep_set_ESVs, taxonomy, by = "ESV")
rownames(taxonomy) <- taxonomy$ESV_ID

# Read in your map file w/sample IDs as rownames
map <- read.delim(map.fp, row.names = 1, as.is = TRUE) 

input <- list()
input$data_loaded <- seqtab.t
input$taxonomy_loaded <- taxonomy
input$map_loaded <- map

# Number of reads per sample
sort(colSums(input$data_loaded))

# Now the input should be compatible with MCTOOLSR, save files
save(input, rep_set_ESVs, file = paste0(table.fp, "/input_tutorial.rda"))

# Also export files as .txt
write.table(seqtab.t, file = paste0(table.fp, "/seqtab_final.txt"),
            sep = "\t", row.names = TRUE, col.names = NA)
write.table(tax, file = paste0(table.fp, "/tax_final.txt"), 
            sep = "\t", row.names = TRUE, col.names = NA)

#' ### Post-pipeline considerations
#' After following this pipline, you will need to think about the following in downstream applications:
#' 1. Remove mitochondrial and chloroplast sequences, for example with mctoolsr:
input_filt <- filter_taxa_from_input(input, taxa_to_remove = "Chloroplast")
input_filt <- filter_taxa_from_input(input_filt, taxa_to_remove = "Mitochondria")
#' 2. Normalize or rarefy your ESV table

#' You can also now transfer over the output files onto your local computer
#'
#'  .rda files are loaded into R like this: 

load(file = "filepath/input_tutorial.rda")
