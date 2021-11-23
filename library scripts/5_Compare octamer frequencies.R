## Set wd to location of the 2 bam files
setwd("D:/library")

#load the required R packages
library(Rsamtools)
library(parallel)
library(stringr)
library(data.table)
#install_github("cran/refGenome")
library(refGenome) 
library(VarCon)
library(ModCon)

## Increase memory limit to maximum
memory.limit(9990000)

## Define needed functions
Numextract <- function(string){
    unlist(regmatches(string, gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
}

## Get highest, second highest... forth highest number from vector
maxi2 <- function(x){ max(sort(x, decreasing = T)[-1])}
maxi3 <- function(x){ max(sort(x, decreasing = T)[c(-1,-2)])}
maxi4 <- function(x){ max(sort(x, decreasing = T)[c(-1,-2,-3)])}




## Get list of BAM files
all_bam_files <- list.files(path = ".", pattern = ".bam")
all_bam_files <- all_bam_files[!str_detect(all_bam_files, "bai")]

## Select first Bam file, of reads from RNA fragments upon transfection of the reporter plasmid
current_bam_file <- all_bam_files[1]

## Access header of Bam file
BAMheader <- scanBamHeader(current_bam_file, what=c("targets", "text"))
BAMheader <- BAMheader[[1]]
BAMheader <- sort(BAMheader$targets)
BAMheader <- BAMheader[!str_detect(names(BAMheader), "\\.")]

## Define "chromosome" of interest, which is called "plasmid" in this case
scafold <-  "plasmid"

## Define which information to gather from BAM file
what <- c("rname", "strand", "pos", "qwidth", "cigar","flag", "seq", "qname")

## define, which regions to exctract from BAM file (in this case again only reads from "plasmid")
if(scafold == "plasmid") which <-  IRangesList("plasmid" = IRanges(0, as.numeric(BAMheader[scafold])))

## Scan the BAM file for the information of reads, which were aligned to the region
param <- ScanBamParam(what = what, which = which)
bam <- scanBam(current_bam_file,index=current_bam_file, param=param) 
bam <- bam[[1]]

## Create data frame from BAM file
bam_df <- data.frame(seqid= bam$rname, strand= bam$strand, lstart = bam$pos, read_length = bam$qwidth, cigar= bam$cigar, seq= bam$seq, orientation= bam$flag, read_name= bam$qname )
bam_df$cigar <- as.character(bam_df$cigar)

## Delete unmatched reads
bam_df <- bam_df[!bam_df$orientation %in% c(73, 133, 89, 121,
                                            165, 181, 101, 117,
                                            153, 185, 69, 137,
                                            67, 131, 97, 145, 
                                            65, 129, 113, 177),]


## Detect reads which show gap in their alignment
bam_df$spliced <- str_detect(bam_df$cigar, "N")
bam_df_gapped <- bam_df[bam_df$spliced,]

## Report number of total reads before filtering
message(paste0("Number of reads before filtering for gapped reads (",scafold,"): ", as.character(nrow(bam_df))))

## Save data so far
save(bam_df_gapped, file="bam_df_gapped_BANDE_merged")
#load("bam_df_gapped_BANDE_merged")

## Remove reads with insertions or deletions, after checking the # of reads with this features
bam_df_gapped$deletions_and_others <- str_detect(bam_df_gapped$cigar, "[IDSHPX=]")
gapped <- bam_df_gapped[!bam_df_gapped$deletions_and_others,]
gapped$deletions_and_others <- NULL
gapped$spliced <- NULL

## Check whether there is always first matched nts than gap in alignment than matched nt,
## because there should not be multiple splicing events in one read
gapped$cigar <- as.character(gapped$cigar)

## Access only the cigar letters from the cigar string
gapped$cigar_letters <- as.character(lapply(gapped$cigar, function(x){   gsub("[[:digit:]]","",x) }))

## Select reads with a single gap in the alignment and meassure the loss of data
gapped2 <- gapped[gapped$cigar_letters == "MNM",]
gapped <- gapped2

## Save data so far
save(gapped, file="gapped")
#load("gapped")

## Rename uniquely gapped reads 
single_gapped <- gapped

## Calculate lend of the gap from lstart
single_gapped$first_matched_bases <- as.numeric(str_extract(single_gapped$cigar, "[0-9]+"))
single_gapped$lend  <- single_gapped$lstart+single_gapped$first_matched_bases-1
single_gapped$first_matched_bases <- NULL

## Access gap length from cigar string
single_gapped$gapped_matched_bases <- as.numeric(lapply(single_gapped$cigar, function(x){ Numextract(x)[2] }))
single_gapped$gaplen <- single_gapped$gapped_matched_bases

## Remove redundant columns
single_gapped$gapped_matched_bases <- NULL
single_gapped$cigar_letters<- NULL

## Get acceptor coordinates 
single_gapped$rstart <- single_gapped$lend + single_gapped$gaplen +1

## Calculate rend coordinate of the gap in the alignment from info of the cigar
single_gapped$gapped_matched_bases <- as.numeric(lapply(single_gapped$cigar, function(x){ Numextract(x)[3] }))
single_gapped$rend <- single_gapped$rstart + single_gapped$gapped_matched_bases -1
single_gapped$gapped_matched_bases <- NULL

single_gapped$gene_strand <- ""

load("library_type")

## Normalize alignment strand for paired and unpaired depending on library type
if(library_type == "F"){
    
    ## Check if this is correct for paired
    single_gapped$gene_strand[single_gapped$orientation == 0]  <- "+"
    single_gapped$gene_strand[single_gapped$orientation == 16] <- "-"
    
    single_gapped$gene_strand[single_gapped$orientation == 99] <- "+"
    single_gapped$gene_strand[single_gapped$orientation == 147] <- "+"
    
    single_gapped$gene_strand[single_gapped$orientation == 83] <- "-"
    single_gapped$gene_strand[single_gapped$orientation == 163] <- "-"
    
}

if(library_type == "R"){
    
    ## Check if this is correct for paired
    single_gapped$gene_strand[single_gapped$orientation == 0]  <- "-"
    single_gapped$gene_strand[single_gapped$orientation == 16] <- "+"
    
    single_gapped$gene_strand[single_gapped$orientation == 99] <- "-"
    single_gapped$gene_strand[single_gapped$orientation == 147] <- "-"
    
    single_gapped$gene_strand[single_gapped$orientation == 83] <- "+"
    single_gapped$gene_strand[single_gapped$orientation == 163] <- "+"
}


## Normalize gap coordinates
single_gapped$norm_lstart <- single_gapped$lstart
single_gapped$norm_lstart[single_gapped$gene_strand == "-"] <- single_gapped$rend[single_gapped$gene_strand == "-"]

single_gapped$norm_lend <- single_gapped$lend
single_gapped$norm_lend[single_gapped$gene_strand == "-"] <- single_gapped$rstart[single_gapped$gene_strand == "-"]

single_gapped$norm_rstart <- single_gapped$rstart
single_gapped$norm_rstart[single_gapped$gene_strand == "-"] <- single_gapped$lend[single_gapped$gene_strand == "-"]

single_gapped$norm_rend <- single_gapped$rend
single_gapped$norm_rend[single_gapped$gene_strand == "-"] <- single_gapped$lstart[single_gapped$gene_strand == "-"]

single_gapped$lstart <- NULL
single_gapped$lend <- NULL
single_gapped$rstart <- NULL
single_gapped$rend <- NULL

## Check abundance of different gaps
table(single_gapped$norm_lend)

## Mark reads with "correct" junction, which indicates downstream SD or upstream SD usage
single_gapped$downSDusage <- 0
single_gapped$downSDusage[single_gapped$norm_lend == 1157 | single_gapped$norm_rstart == 1157] <- 1

single_gapped$upSDusage <- 0
single_gapped$upSDusage[single_gapped$norm_lend == 1094 | single_gapped$norm_rstart == 1094] <- 1

## Count reads with downstream SD usage
table(single_gapped$downSDusage)
# 0       1 
# 5743358 1005604 

## Count reads with upstream SD usage
table(single_gapped$upSDusage)
# 0       1 
# 1011296 5737666 

## Only select reads with downstream SD usage, to access sequence in octamer library
downsd <- single_gapped[single_gapped$downSDusage == 1,]

## Use regular expression with 5nt anchor sequences to find position of octamer library in each read
downsd$octa_pos <- gregexpr("ATTGG[ATGC][ATGC][ATGC][ATGC][ATGC][ATGC][ATGC][ATGC]CCTAT", downsd$seq)

## Access the sequence of the octamer library in each read
downsd$octa_seq <- as.character(lapply(c(1:nrow(downsd)), function(x){
    
    oct <- downsd$seq[x]
    posi <- as.numeric(downsd$octa_pos[x][[1]])
    octa <- ""
    
    if(posi[1] != -1){
        
        posi <- max(posi)
        
        octa <- substr(oct, posi+5, posi+12)
         
    }
    
    octa

}))

## Save the data so far
save(downsd, file="downsd")

## Count abundance of every octamer found
downsd$c <- 1
sds <- tapply(downsd$c, downsd$octa_seq, sum)
sds <- data.frame(sds)
sds$octamer <- row.names(sds)
row.names(sds) <- NULL

## Count and remove reads with no octamer found
sds <- sds[sds$octamer != "",]

## Save data for later
write.csv2(sds, file="octamers_band.csv", row.names = F)







##############
### Counting octamers in plasmid sequencing data
##############

## Select BAM file with the plasmid sequencing data
current_bam_file <- all_bam_files[2]

## Access header of Bam file
BAMheader <- scanBamHeader(current_bam_file, what=c("targets", "text"))
BAMheader <- BAMheader[[1]]
BAMheader <- sort(BAMheader$targets)
BAMheader <- BAMheader[!str_detect(names(BAMheader), "\\.")]

## Define "chromosome" of interest, which is called "plasmid" in this case
scafold <-  "plasmid"

## Define which information to gather from BAM file
what <- c("rname", "strand", "pos", "qwidth", "cigar","flag", "seq", "qname")

## Define, which regions to exctract from BAM file (in this case again only reads from "plasmid")
if(scafold == "plasmid") which <-  IRangesList("plasmid" = IRanges(0, as.numeric(BAMheader[scafold])))

## Scan the BAM file for the information of reads, which were aligned to the region
param <- ScanBamParam(what = what, which = which)
bam <- scanBam(current_bam_file,index=current_bam_file, param=param) 
bam <- bam[[1]]

## Create data frame from BAM file
bam_df <- data.frame(seqid= bam$rname, strand= bam$strand, lstart = bam$pos, read_length = bam$qwidth, cigar= bam$cigar, seq= bam$seq, orientation= bam$flag, read_name= bam$qname )
bam_df$cigar <- as.character(bam_df$cigar)

## Delete unmatched reads
bam_df <- bam_df[!bam_df$orientation %in% c(73, 133, 89, 121,
                                            165, 181, 101, 117,
                                            153, 185, 69, 137,
                                            67, 131, 97, 145, 
                                            65, 129, 113, 177),]


## Detect reads which show gap in their alignment
bam_df$spliced <- str_detect(bam_df$cigar, "N")
bam_df_gapped <- bam_df[!bam_df$spliced,]

## Report number of total reads before filtering
message(paste0("Number of reads before filtering (",scafold,"): ", as.character(nrow(bam_df))))

## Save data so far
save(bam_df_gapped, file="bam_df_gapped_PLASMID_merged")
#load("bam_df_gapped_BANDE")

## Remove reads with insertions or deletions, after checking the # of reads with this features
bam_df_gapped$deletions_and_others <- str_detect(bam_df_gapped$cigar, "[IDSHPX=]")
gapped <- bam_df_gapped[!bam_df_gapped$deletions_and_others,]
gapped$deletions_and_others <- NULL
gapped$spliced <- NULL
gapped$cigar <- as.character(gapped$cigar)

## Rename plasmid reads
plasmid <- gapped

## Use regular expression with 5nt anchor sequences to find position of octamer library in each read
plasmid$octa_pos <- gregexpr("ATTGG[ATGC][ATGC][ATGC][ATGC][ATGC][ATGC][ATGC][ATGC]CCTAT", plasmid$seq)

## Save the data so far
save(plasmid, file="plasmid_plasmid")
#load("plasmid_plasmid")

## Access the sequence of the octamer library in each read
plasmid$octa_seq <- as.character(lapply(c(1:nrow(plasmid)), function(x){
    
    oct <- plasmid$seq[x]
    posi <- as.numeric(plasmid$octa_pos[x][[1]])
    octa <- ""
    
    if(posi[1] != -1){
        
        posi <- max(posi)
        
        octa <- substr(oct, posi+5, posi+12)
        
    }
    
    octa
    
}))

## Save the data so far
save(plasmid, file="plasmid_plasmid2")

## Count abundance of every octamer found
plasmid$c <- 1
sds <- tapply(plasmid$c, plasmid$octa_seq, sum)
sds <- data.frame(sds)
sds$octamer <- row.names(sds)
row.names(sds) <- NULL

## Count and remove reads with no octamer found
sds <- sds[sds$octamer != "",]

## Save data for later
write.csv2(sds, file="octamers_plasmid.csv", row.names = F)






##############
### Compare band and plasmid octamer frequencies
##############

## Load band and plasmid octamer frequencies
ba <- read.csv2("D:/library/octamers_band.csv")
pl <- read.csv2("D:/library/octamers_plasmid.csv")

## Create combined list of all found unique octamers
octameres <- data.frame(octamer=c(unique(c(ba$octamer, pl$octamer))))

## Enter frequency of each octamer in band or plasmid sequencing data
octameres$freq_plasmid <- pl$sds[match(octameres$octamer, pl$octamer)]
octameres$freq_band <- ba$sds[match(octameres$octamer, ba$octamer)]
octameres[is.na(octameres)] <- 0

## Calculate the HZEI of a completely neutral insertion between the two donor sites
octaneutral <- "CCTATTGGCCTATTGGCCTATTGGCCTATTGGCCTATTGG"
octameres_merged$exonHZEI <- calculateHZEIint(octaneutral)  # = -4.526667

## Create new sequence with 2 neutral octamers upstream and downstream of the determined octamer seq
octameres_merged$exon_with_okta <- paste0("CCTATTGGCCTATTGG",octameres_merged$octamer,"CCTATTGGCCTATTGG")

## Calculate total HZEI for the exon sequence 
octameres_merged$exon_with_oktaHZEI  <-  as.numeric(lapply(octameres_merged$exon_with_okta, calculateHZEIint))

## Calcualte HZEI difference between completely neutral exon and the exon with the respective octamer in the middel
octameres_merged$HZ_EI_diff_exon <- octameres_merged$exon_with_oktaHZEI - octameres_merged$exonHZEI

## Calculate enrichment score
norm_total_band <- octameres_merged$freq_band[octameres_merged$freq_band >0 & freq_band$freq_plasmid >0]
norm_total_plasmid <- octameres_merged$freq_plasmid[octameres_merged$freq_band >0 & freq_band$freq_plasmid >0]
octameres_merged$enrichment_score <- (octameres_merged$freq_band / octameres_merged$freq_plasmid) / ( norm_total_band / norm_total_plasmid)

## Save the data
write.csv2(octameres_merged, file="octameres_merged_exonHZEI.csv", row.names=F)


