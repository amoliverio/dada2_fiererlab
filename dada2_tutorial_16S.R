#'# dada2 tutorial with MiSeq dataset for Fierer Lab 

#+ setup, include=FALSE
knitr::opts_chunk$set(eval = FALSE, include = TRUE)

#' This version runs the dada2 workflow for Big Data (paired-end) from Rstudio on the microbe server.
#' 
#' * More information  on this pipeline can be found here: [https://benjjneb.github.io/dada2/bigdata_paired.html](https://benjjneb.github.io/dada2/bigdata_paired.html)
#'    
#' We suggest opening the dada2 tutorial online to understand more about each step
#'
#' || <span> ||
#' |-| :--- |-|
#' || _**NOTE:**_ there is a slightly different pipeline for ITS and Non-"Big data" that takes a closer look at each step here: [https://benjjneb.github.io/dada2/tutorial.html](https://benjjneb.github.io/dada2/tutorial.html) ||
#' || <span> ||
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
#' 
#' | _**NOTE:**_ if you are running on your local computer make sure you have idemp installed. Found here and it is a very quick install: [https://github.com/yhwu/idemp](https://github.com/yhwu/idemp) |
#' | :--- |
#' 

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
#' 
#' Note: If you are not working from microbe server, you will need to change the file paths for idemp and cutadapt to where they are stored on your computer/server.
#' 
#' For this tutorial we will be working with some samples that we obtained 16S amplicon data for, from a Illumina Miseq run. 
#' 
# Set up pathway to idemp (demultiplexing tool) and test
idemp <- "/usr/bin/idemp" # CHANGE ME if not on microbe
system2(idemp) # Check that idemp is in your path and you can run shell commands from R

# Set up pathway to cutadapt (primer trimming tool) and test
cutadapt <- "/usr/local/Python27/bin/cutadapt" # CHANGE ME if not on microbe
system2(cutadapt, args = "--version") # Check by running shell command from R

# Set path to shared data folder and contents
data.fp <- "/data/shared/2019_02_20_MicrMethods_tutorial"

# List all files in shared folder to check path
list.files(data.fp)

# Set file paths for barcodes file, map file, and fastqs
    # barcodes need to have 'N' on the end of each 12bp sequence for compatability
barcode.fp <- file.path(data.fp, "barcode_demultiplex_short.txt") # .txt file: barcode </t> sampleID
map.fp <- file.path(data.fp, "Molecular_Methods_18_515fBC_16S_Mapping_File_SHORT_vFinal_Fierer_10252018.txt")
I1.fp <- file.path(data.fp, "Undetermined_S0_L001_I1_001.fastq.gz") 
R1.fp <- file.path(data.fp, "Undetermined_S0_L001_R1_001.fastq.gz") 
R2.fp <- file.path(data.fp, "Undetermined_S0_L001_R2_001.fastq.gz") 

#' Set up file paths in YOUR directory where you want data; 
#' you do not need to create the subdirectories but they are nice to have
#' for organizational purposes. 

project.fp <- "/data/YOUR_USERNAME/MicroMethods_dada2_tutorial" # CHANGE ME to project directory

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
 
# Change names of unassignable reads so they are not included in downstream processing
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

# Name the N-filtered files to put them in filtN/ subdirectory
fnFs.filtN <- file.path(preprocess.fp, "filtN", basename(fnFs))
fnRs.filtN <- file.path(preprocess.fp, "filtN", basename(fnRs))

# filter Ns from reads and put them into the filtN directory
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE) # CHANGE multithread to FALSE on Windows (here and elsewhere in the program)


#' #### Prepare the primers sequences and custom functions for analyzing the results from cutadapt
#' Assign the primers you used to "FWD" and "REV" below. Note primers should be not be reverse complemented ahead of time. Our tutorial data uses 515f and 806br those are the primers below. Change if you sequenced with other primers.
#' 

# set up the primer sequences to pass along to cutadapt
FWD <- "GTGYCAGCMGCCGCGGTAA"  ## CHANGE ME # this is 515f
REV <- "GGACTACNVGGGTWTCTAAT"  ## CHANGE ME # this is 806Br

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
#' for the first sample, as a check. There may be some primers here, we will remove them below using cutadapt.
#' 
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#' #### Remove primers with cutadapt and assess the output

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

# As a sanity check, we will check for primers in the first cutadapt-ed sample:
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

# Move R1 and R2 from trimmed to separate forward/reverse sub-directories
fnFs.Q <- file.path(subF.fp,  basename(fnFs)) 
fnRs.Q <- file.path(subR.fp,  basename(fnRs))
file.rename(from = fnFs.cut, to = fnFs.Q)
file.rename(from = fnRs.cut, to = fnRs.Q)

# File parsing
filtpathF <- file.path(subF.fp, "filtered") # files go into preprocessed_F/filtered/
filtpathR <- file.path(subR.fp, "filtered") # ...
fastqFs <- sort(list.files(subF.fp, pattern="fastq.gz"))
fastqRs <- sort(list.files(subR.fp, pattern="fastq.gz"))
if(length(fastqFs) != length(fastqRs)) stop("Forward and reverse files do not match.")

#' ### 1. FILTER AND TRIM FOR QUALITY
#' 
#' Before chosing sequence variants, we want to trim reads where their quality scores 
#' begin to drop (the `truncLen` and `truncQ` values) and remove any low-quality reads 
#' that are left over after we have finished trimming (the `maxEE` value).

#' ---
#' **WARNING:** THESE PARAMETERS ARE NOT OPTIMAL FOR ALL DATASETS. Make sure you determine
#' the trim and filtering parameters for your data. The following settings are 
#' generally appropriate for MiSeq runs that are 2x150 bp.
#' 
#' You will want to change this depending on run chemistry and quality: for 2x250 bp runs you can try truncLen=c(240,160) as per the dada2 tutorial if your reverse reads drop off in quality and higher, for example, truncLen=c(240,240) if they do not.
#' **For ITS data:** Due to the expected variable read lengths in ITS data you should run this command without the trunclen parameter. See here for more information and appropriate parameters for ITS data:https://benjjneb.github.io/dada2/ITS_workflow.html.
#' 
#' **NOTE from dada2 tutorial:**
#' "If there is only one part of any amplicon bioinformatics workflow on which you spend time considering the parameters, it should be filtering! The parameters ... are not set in stone, and should be changed if they don’t work for your data. If too few reads are passing the filter, increase maxEE and/or reduce truncQ. If quality drops sharply at the end of your reads, reduce truncLen. If your reads are high quality and you want to reduce computation time in the sample inference step, reduce  maxEE."
#' ---
#'

filterAndTrim(fwd=file.path(subF.fp, fastqFs), filt=file.path(filtpathF, fastqFs),
              rev=file.path(subR.fp, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
              truncLen=c(150,140), maxEE=1, truncQ=11, maxN=0, rm.phix=TRUE,
              compress=TRUE, verbose=TRUE, multithread=TRUE)

#' ### 2. INFER sequence variants
#' In this part of the pipeline dada2 will learn to distinguish error from biological 
#' differences using a subset of our data as a training set. After it understands the 
#' error rates, we will reduce the size of the dataset by combining all identical 
#' sequence reads into "unique sequences". Then, using the dereplicated data and 
#' error rates, dada2 will infer the sequence variants (OTUs) in our data. Finally, 
#' we will merge the coresponding forward and reverse reads to create a list of the 
#' fully denoised sequences and create a sequence table from the result.

#' #### Housekeeping step - set up and verify the file names for the output:
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
saveRDS(seqtab, paste0(table.fp, "/seqtab.rds"))

#' ### 3. REMOVE Chimeras and ASSIGN Taxonomy
#' Although dada2 has searched for indel errors and subsitutions, there may still be chimeric
#' sequences in our dataset (sequences that are derived from forward and reverse sequences from 
#' two different organisms becoming fused together during PCR and/or sequencing). To identify 
#' chimeras, we will search for rare sequence variants that can be reconstructed by combining
#' left-hand and right-hand segments from two more abundant "parent" sequences. After removing chimeras, we will use a taxonomy database to train a classifer-algorithm
#' to assign names to our sequence variants.
#' 
#' For the tutorial 16S, we will assign taxonomy with Silva db v132, but you might want to use other databases for your data. Below are paths to some of the databases we use often. (If you are on your own computer you can download the database you need from this link https://benjjneb.github.io/dada2/training.html):
#' 
#' 16S bacteria and archaea (SILVA db): /db_files/dada2/silva_nr_v132_train_set.fa
#' ITS fungi (UNITE db): /db_files/dada2/unite_general_release_dynamic_02.02.2019.fasta
#' 18S protists (PR2 db): /db_files/dada2/pr2_version_4.11.1_dada2.fasta
#' 

# Read in RDS 
st.all <- readRDS(paste0(table.fp, "/seqtab.rds"))

# Remove chimeras
seqtab <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE)

# Assign taxonomy
tax <- assignTaxonomy(seqtab, "/db_files/dada2/silva_nr_v132_train_set.fa",
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

# Add ESV numbers to taxonomy
taxonomy <- as.data.frame(tax)
taxonomy$ESV <- as.factor(rownames(taxonomy))
taxonomy <- merge(rep_set_ESVs, taxonomy, by = "ESV")
rownames(taxonomy) <- taxonomy$ESV_ID
taxonomy_for_mctoolsr <- unite_(taxonomy, "taxonomy", c("Kingdom", "Phylum", "Class", "Order","Family", "Genus", "ESV_ID"), sep = ";")

# Merge taxonomy and table
seqtab_wTax <- merge(seqtab.t, taxonomy_for_mctoolsr, by = 0)
seqtab_wTax$ESV <- NULL

# Set name of table in mctoolsr format and save
out_fp <- paste0(table.fp, "/seqtab_wTax_mctoolsr.txt")
names(seqtab_wTax)[1] = "#ESV_ID"
write("#Exported for mctoolsr", out_fp)
suppressWarnings(write.table(seqtab_wTax, out_fp, sep = "\t", row.names = FALSE, append = TRUE))

# Also export files as .txt
write.table(seqtab.t, file = paste0(table.fp, "/seqtab_final.txt"),
            sep = "\t", row.names = TRUE, col.names = NA)
write.table(tax, file = paste0(table.fp, "/tax_final.txt"), 
            sep = "\t", row.names = TRUE, col.names = NA)

#' ### Post-pipeline considerations
#' After following this pipline, you will need to think about the following in downstream applications (example with 'mctoolsr' R package below):
#' 
#' 1. Remove mitochondrial and chloroplast sequences
#' 2. Remove reads assigned as eukaryotes
#' 3. Remove reads that are unassigned at domain level

input_filt <- filter_taxa_from_input(input, taxa_to_remove = c("Chloroplast","Mitochondria", "Eukaryota"))
input_filt <- filter_taxa_from_input(input_filt, at_spec_level = 1, taxa_to_remove = "NA")

#' 4. Normalize or rarefy your ESV table
#' You can also now transfer over the output files onto your local computer
#'
#'  .rda files are loaded into R like this: 

load(file = "filepath/input_tutorial.rda")
