## Setting pipline required values 
mainDirectory <- "~/Desktop/HomeOffice/SD paper/same for endothelial cells"
Projectshortcut <- "endo"

## Process BAM files to gap table
library(spliceSites)
library(refGenome)

## Define working directory and annotation directory
prdir <-"/path/to/Bam/files/"
setwd(prdir)

## Load bam files and bam indices
bam <- list.files(path = ".", pattern = ".bam")
bam <- bam[!str_detect(bam, "bai")]

## Create own bam index
bam_idx <- character()

for(i in 1:length(bam)){
  
  bam_idx[i] <- paste0(bam[i],".bai")
  reader <- bamReader(bam[i])
  createIndex(reader,bam_idx[i])
  loadIndex(reader,bam_idx[i])
  
  reader <- bamReader(bam[i], idx=TRUE)
  count <- bamCountAll(reader)
  print(bam[i])
  print(count)
}

## Create empty list for gaptables per BAM file
end_list <- list()

## For every BAM file create Gaptable
for(i in 1:length(bam)){
  grs <- readPooledBamGaps(bam[i], idxInfiles = bam_idx[i])
  end_list[[bam[i]]] <- as.data.frame(grs)
}


## Save the list of gap tables
save(end_list,file=file.path(prdir, "end_list"))


## Create new directory for the occasion
setwd(mainDirectory)

## Load packages
library(spliceSites)
library(refGenome) 

## Load RNA-seq alignment data
load("end_list")

## Load reference FASTA GRCh38. To generate this file, please download the human reference fasta file
## from Ensembl, and load it into the working directory e.g. with the function of the package VarCon
## "prepareReferenceFasta(filepath)" and rename the chromosomes according to the nomenclature in "end_list"
load("chromosome_sequences_hg38")

## Load annotation GTF file for GRCh38 from Ensembl version 100, can be downloaded from Ensembl archive
ens_gtf <- "Homo_sapiens.GRCh38.100.gtf"
beg <- ensemblGenome()
basedir(beg) <- "C:/Users/Jo/Desktop/HomeOffice/reference"
read.gtf(beg, ens_gtf) 
juc <- getSpliceTable(beg)
rm(ens_gtf)

## Create new list with processed end_list elements
end_list_processed <- list()

## Only choose wildtype endothelial cells
end_list <- end_list[c(1,5,9,13)]

## Calculate rpmgs
for(i in 1:length(end_list)){
  cur <-  end_list[[i]]
  cur$name <- names(end_list)[i]
  cur$rpmg <- round((cur$nAligns/sum(cur$nAligns))*1000000,digits=3)
  end_list_processed[[i]] <- cur
  
}

## Save the new end_list with rpmgs
save(end_list_processed, file=paste0(Projectshortcut,("_end_list_processed")))

## Create new object for annoatated single BAM file information
anno <- list()

## For every BAM file
for(grs_element_number in 1:length(end_list_processed)){ 
  
  ## Create gapsites object
  ref_table <- end_list_processed[[grs_element_number]]
  file_name <- substr((ref_table[1,14]),43,47)
  
  ref_table$nAligns <- as.numeric(ref_table$nAligns)
  ref_table$seqid <- as.character(ref_table$seqid)
  
  ## Create a gapSites Object from the data frame for annotation etc
  rtbg <- gapSites( seqid=ref_table$seqid,lstart=ref_table$lstart,lend=ref_table$lend,
                    rstart=ref_table$rstart,rend=ref_table$rend,gaplen=ref_table$gaplen,
                    nr_aligns=ref_table$nAligns)
  
  
  ## Add maxent scores
  mes <- load.maxEnt()
  rtbgx <- addMaxEnt(rtbg, s, mes)
  rtbgx <- addHbond(rtbgx, s)
  
  ## Add annotation
  annotation(rtbgx) <- annotate(rtbgx, juc)
  strand(rtbgx) <- getAnnStrand(rtbgx)
  rtbgx <- addGeneAligns(rtbgx)
  
  ## Combine with annotation
  rtbg <- annotation(rtbgx)
  
  ## Convert gapSItes-Data to data frame
  end_table <- as.data.frame.gapSites(rtbg)
  end_table2 <- as.data.frame.gapSites(rtbgx)
  
  
  ref_table$ID <- paste(ref_table$seqid,ref_table$lend,ref_table$rstart,ref_table$gaplen,ref_table$nAligns)
  end_table$ID <-  paste(end_table$seqid,end_table$lend,end_table$rstart,end_table$gaplen,end_table$nAligns)
  end_table2$ID <-  paste(end_table2$seqid,end_table2$lend,end_table2$rstart,end_table2$gaplen,end_table2$nAligns)
  
  end_table$nlstart <- ref_table$nlstart[match(end_table$ID,ref_table$ID)]
  end_table$qsm <- ref_table$qsm[match(end_table$ID,ref_table$ID)]
  end_table$gqs <- ref_table$gqs[match(end_table$ID,ref_table$ID)]
  end_table$nProbes <- 1
  end_table$name <- ref_table$name[match(end_table$ID,ref_table$ID)]
  end_table$strand <- end_table2$strand[match(end_table$ID,end_table2$ID)]
  end_table$strand.x <- NULL
  end_table$strand.y <- NULL
  end_table$meStrand <- NULL
  
  ab <- na.omit(end_table)
  ab <- unique(ab)
  
  abs <- end_table[!end_table$id %in% ab$id,]
  
  ## Add Sequence of Donor site
  ab$number <- 1:nrow(ab)
  tbl <- ab
  
  ## Access donor site sequence
  ltbl <- tbl
  ltbl$begin <- tbl$lend-2
  ltbl$end <- tbl$lend+8
  
  ltbl$begin[ltbl$strand =="-"] <- tbl$rstart[ltbl$strand =="-"]-8
  ltbl$end[ltbl$strand =="-"] <- tbl$rstart[ltbl$strand =="-"]+2
  
  # Sort table
  ltbl <- ltbl[order(ltbl$seqid, ltbl$begin), ]
  seqnames <- sort(unique(ltbl$seqid))
  n <- length(seqnames)
  dnames <- names(s)
  s$MT <- NULL
  res_list <- list()
  
  library( seqinr)
  
  ## Get the sequence of a specific position of a specific chromosome and gene
  j <- 1
  for(i in 1:n)
  {
    lg <- match(seqnames[i], dnames)
    if(!is.na(lg))
    {
      dt <- ltbl[ltbl$seqid== seqnames[i], ]
      test2 <- Views(s[[lg]], start=dt$begin, end=dt$end)
      test <- as.character(test2)
      test2 <- test2@ranges
      test2@NAMES <- as.character(test)
      chr_res <- as.data.frame(test2)
      chr_res$chr <- seqnames[i]
      chr_res$ID <- paste(chr_res$chr, chr_res$start, chr_res$end)
      
      res_list[[j]] <- chr_res
      j <- j + 1
      
    }else{
      cat("Seqname mismatch: ", i, " ", seqnames[i])
    }
  }
  
  library(data.table)
  
  res <- rbindlist(res_list)
  ltbl$ID <- paste(ltbl$seqid , ltbl$begin, ltbl$end)
  
  ltbl$lseq <- res$names[match(ltbl$ID, res$ID)]
  
  ## Minus strand sequecne manipulation
  ltbl$lseq[ltbl$strand=="-"] <- gsub("A","t",ltbl$lseq[ltbl$strand=="-"])
  ltbl$lseq[ltbl$strand=="-"] <- gsub("C","g",ltbl$lseq[ltbl$strand=="-"])
  ltbl$lseq[ltbl$strand=="-"] <- gsub("T","a",ltbl$lseq[ltbl$strand=="-"])
  ltbl$lseq[ltbl$strand=="-"] <- gsub("G","c",ltbl$lseq[ltbl$strand=="-"])
  
  ltbl$lseq <- toupper(ltbl$lseq)
  ltbl$lseq[ltbl$strand=="-"] <- reverse(ltbl$lseq[ltbl$strand=="-"])
  
  ab$donor_seq <-  ltbl$lseq[match(ab$number, ltbl$number)]
  ab$donor_site <- substr(ab$donor_seq,4,5)
  
  ## Test for Donor site frequency on minus strand as controll
  table(ab$donor_site[ab$strand=="-"])  
  
  
  
  ## Add Sequence of Acceptor site
  tbl <- ab
  ltbl <- tbl
  ltbl$begin <- tbl$rstart-20
  ltbl$end <- tbl$rstart+2
  
  ltbl$begin[ltbl$strand =="-"] <- tbl$lend[ltbl$strand =="-"]-2
  ltbl$end[ltbl$strand =="-"] <- tbl$lend[ltbl$strand =="-"]+20
  
  ## Sort table
  ltbl <- ltbl[order(ltbl$seqid, ltbl$begin), ]
  seqnames <- sort(unique(ltbl$seqid))
  n <- length(seqnames)
  dnames <- names(s)
  s$MT <- NULL
  res_list <- list()
  
  library( seqinr)
  
  ## Get the sequence of a specific position of a specific chromosome and gene
  j <- 1
  for(i in 1:n)
  {
    lg <- match(seqnames[i], dnames)
    if(!is.na(lg))
    {
      dt <- ltbl[ltbl$seqid== seqnames[i], ]
      test2 <- Views(s[[lg]], start=dt$begin, end=dt$end)
      test <- as.character(test2)
      test2 <- test2@ranges
      test2@NAMES <- as.character(test)
      chr_res <- as.data.frame(test2)
      chr_res$chr <- seqnames[i]
      chr_res$ID <- paste(chr_res$chr, chr_res$start, chr_res$end)
      
      res_list[[j]] <- chr_res
      j <- j + 1
      
    }else{
      cat("Seqname mismatch: ", i, " ", seqnames[i])
    }
  }
  
  library(data.table)
  
  res <- rbindlist(res_list)
  ltbl$ID <- paste(ltbl$seqid , ltbl$begin, ltbl$end)
  
  ltbl$lseq <- res$names[match(ltbl$ID, res$ID)]
  
  ## Minus strand sequecne manipulation
  ltbl$lseq[ltbl$strand=="-"] <- gsub("A","t",ltbl$lseq[ltbl$strand=="-"])
  ltbl$lseq[ltbl$strand=="-"] <- gsub("C","g",ltbl$lseq[ltbl$strand=="-"])
  ltbl$lseq[ltbl$strand=="-"] <- gsub("T","a",ltbl$lseq[ltbl$strand=="-"])
  ltbl$lseq[ltbl$strand=="-"] <- gsub("G","c",ltbl$lseq[ltbl$strand=="-"])
  
  ltbl$lseq <- toupper(ltbl$lseq)
  ltbl$lseq[ltbl$strand=="-"] <- reverse(ltbl$lseq[ltbl$strand=="-"])
  
  ab$acceptor_seq <-  ltbl$lseq[match(ab$number, ltbl$number)]
  ab$acceptor_site <- substr(ab$acceptor_seq,19,20)
  
  ab$number <- NULL
  
  
  ## Normalize gap coordinates
  ab$norm_lstart <- ab$lstart 
  ab$norm_lstart[ab$strand == "-"] <- ab$rend[ab$strand == "-"]
  ab$norm_lend <- ab$lend
  ab$norm_lend[ab$strand == "-"] <- ab$rstart[ab$strand == "-"]
  
  ## right
  ab$norm_rstart <- ab$rstart 
  ab$norm_rstart[ab$strand == "-"] <- ab$lend[ab$strand == "-"]
  ab$norm_rend <- ab$rend
  ab$norm_rend[ab$strand == "-"] <- ab$lstart[ab$strand == "-"]
  
  ## Delete unnormalized coordinates
  ab$lend <- NULL
  ab$rend <- NULL
  ab$lstart <- NULL
  ab$rstart  <- NULL
  
  ## Normalize the ldiff and rdiff entries of the gaps by the respective strand
  #ldiff
  ab$norm_ldiff <- ab$ldiff
  ab$norm_ldiff[ab$strand == "-"] <- ab$rdiff[ab$strand == "-"]
  ab$ldiff <- ab$norm_ldiff
  #rdiff
  ab$norm_rdiff <- ab$rdiff
  ab$norm_rdiff[ab$strand == "-"] <- ab$ldiff[ab$strand == "-"]
  ab$rdiff <- ab$norm_rdiff
  
  ## Delete transformation variables
  ab$norm_ldiff <- NULL
  ab$norm_rdiff<- NULL
  
  ## The ldiff on the minus strand is multiplied with -1, to assign negative relative distances to an annotated SD with positons on the exon.
  ab$ldiff[ab$strand == "-"]  <- -(ab$ldiff[ab$strand == "-"] )
  ab$rdiff[ab$strand == "-"]  <- -(ab$rdiff[ab$strand == "-"] )
  
  anno[[file_name]] <- ab
  
  
}

save(anno, file=paste0(Projectshortcut,("_annotated_grs_file_sample_seperated")))

## Increase memory
memory.limit(size = 9999999999)
Projectshortcut <- "endo"
load(paste0(Projectshortcut,("_annotated_grs_file_sample_seperated")))

## Tidy up
rm(ab, abs, beg, dt, end_table, end_table2, enjc80, end_list, end_list_processed, uc, ltbl, mes, nn , proj, ref_table, rtbg, rtbgx, rvl, tbl )

library("data.table")

## Generate junction table "fib"
fib <- rbindlist(anno)
fib$id <- NULL
fib$nAligns <- NULL
fib$gptm <- NULL
fib$rpmg <- NULL
fib$s5strand <- NULL
fib$s3strand <- NULL
fib$gene_aligns <- NULL
fib$refid <- NULL
fib$nref <- NULL
fib$first_refid <- NULL
fib$last_refid <- NULL
fib$nadv <- NULL
fib$nrev <- NULL
fib$ID <- NULL
fib$name <- NULL
fib$nlstart <- NULL
fib$qsm <- NULL
fib$gqs <- NULL
fib$norm_lstart <- NULL
fib$norm_rend <- NULL
fib$transcript_id <- NULL
fib$nProbes <- NULL
fib <- unique(fib)

## Create Junction ID
fib$ID <- paste(fib$seqid, fib$gene_id, fib$norm_lend, fib$norm_rstart)

## Save the precursor tables
save(anno, fib, file = paste0(Projectshortcut,("Summarize_objects")))


anno2 <- list()

## For every sample
for(i in 1:length(anno)){
  
  cur_tab <- anno[[i]]
  
  ## Generate new variable, to mark events, which have to be deleted later
  cur_tab$auswahlfaktor <- 1
  
  ## Using the variable auswahlfaktor, gaps with a length of more than 26914 nucleotides are deleted 
  cur_tab$auswahlfaktor[cur_tab$gaplen > 26914] <- NA
  
  ## Gaps which could not be assigned to a gene were marked to be deleted later on
  cur_tab$auswahlfaktor[is.na(cur_tab$gene_name)] <- NA
  
  ## Delet the marked gaps
  cur_tab<- na.omit(cur_tab)
  
  ## Create a variable for gqs transformation
  cur_tab$gqsfaktor <- 1
  
  ## Gqsfactor for different naligns and nlstart combinations is added
  cur_tab$gqsfaktor[(cur_tab$nAligns == 1) ] <- 32
  
  cur_tab$gqsfaktor[(cur_tab$nAligns == 2)&(cur_tab$nlstart == 1) ] <- 16
  cur_tab$gqsfaktor[(cur_tab$nAligns == 2)&(cur_tab$nlstart == 2) ] <- 8
  
  cur_tab$gqsfaktor[(cur_tab$nAligns == 3)&(cur_tab$nlstart == 1) ] <- (32/3)
  cur_tab$gqsfaktor[(cur_tab$nAligns == 3)&(cur_tab$nlstart == 2) ] <- (16/3)
  cur_tab$gqsfaktor[(cur_tab$nAligns == 3)&(cur_tab$nlstart == 3) ] <- (32/9)
  
  cur_tab$gqsfaktor[(cur_tab$nAligns >= 4)&(cur_tab$nlstart == 1) ] <- 8
  cur_tab$gqsfaktor[(cur_tab$nAligns >= 4)&(cur_tab$nlstart == 2) ] <- 4
  cur_tab$gqsfaktor[(cur_tab$nAligns >= 4)&(cur_tab$nlstart == 3) ] <- (8/3)
  cur_tab$gqsfaktor[(cur_tab$nAligns >= 4)&(cur_tab$nlstart == 4) ] <- 2
  cur_tab$gqsfaktor[(cur_tab$nAligns > 4)&(cur_tab$nlstart == 5) ] <- 1.6
  cur_tab$gqsfaktor[(cur_tab$nAligns > 4)&(cur_tab$nlstart == 6) ] <- (4/3)
  cur_tab$gqsfaktor[(cur_tab$nAligns > 4)&(cur_tab$nlstart == 7) ] <- (8/7)
  
  
  ## Controll the outcome
  head(cur_tab)
  
  ## The new gqs is calculated
  cur_tab$gqs <- cur_tab$gqs*cur_tab$gqsfaktor
  
  ## Delete the gqsfactor
  cur_tab$gqsfaktor <- NULL
  
  ## Delete all gaps with a gqs, less than 400, are marked to be deleted
  cur_tab$auswahlfaktor[cur_tab$gqs < 400] <- NA
  
  ## Delet the marked gaps
  cur_tab <- na.omit(cur_tab)
  
  ## Select reads on main scaffolds
  table(cur_tab$donor_site)
  cur_tab <- cur_tab[cur_tab$seqid %in% c(1:24,"X","Y"), ]
  cur_tab <- na.omit(cur_tab)
  cur_tab$id <- 1:nrow(cur_tab)
  cur_tab$seqid <- as.character(cur_tab$seqid)
  cur_tab$strand <- as.character(cur_tab$strand)
  
  
  ## Generate SD aggreagted table
  cur_tab$IDS <- paste(cur_tab$gene_id, cur_tab$norm_lend)
  
  aln <- tapply(cur_tab$rpmg ,cur_tab$IDS, sum)
  adf <- data.frame(aln)
  adf$name <- names(aln)
  
  gnr <- tapply(cur_tab$nAligns ,cur_tab$IDS, sum)
  gnrdf <- data.frame(gnr)
  gnrdf$name <- names(gnr)
  
  ## Creating SD table with summed up positional aligns and GNR
  sd <- data.frame(IDS= cur_tab$IDS,gene_name="", norm_lend=0, pos_rpmg=0, pos_aln=0)
  sd$gene_name <- cur_tab$gene_name[match(sd$IDS,cur_tab$IDS)]
  sd$gene_id <- cur_tab$gene_id[match(sd$IDS,cur_tab$IDS)]
  sd$norm_lend <- cur_tab$norm_lend[match(sd$IDS,cur_tab$IDS)]
  sd$pos_rpmg <- adf$aln[match(sd$IDS,adf$name)]
  sd$pos_aln <- gnrdf[match(sd$IDS,gnrdf$name),1]
  
  ## Add aln
  aln <- tapply(sd$pos_aln ,cur_tab$gene_id, max)
  adf <- data.frame(aln)
  adf$name <- names(aln)
  sd$MRIG_SD <- adf$aln[match(sd$gene_id,adf$name)]
  sd$GNR_SD <- sd$pos_aln/sd$MRIG_SD
  
  ## Make df unique
  sd <- unique(sd)
  
  ## Insert other informations from cur_tab table
  sd$seqid <- cur_tab$seqid[match(sd$IDS,cur_tab$IDS)]
  sd$strand <- cur_tab$strand[match(sd$IDS,cur_tab$IDS)]
  sd$donor_site <- cur_tab$donor_site[match(sd$IDS,cur_tab$IDS)]
  sd$donor_sequence <- cur_tab$donor_seq[match(sd$IDS,cur_tab$IDS)]
  sd$maxent <- 0
  sd$norm_lhbond <- 0
  file_name <- names(anno)[[i]]
  
  ## Create sequence columns for HEX ints
  sd$Donor_HEXseq     <- ""
  sd$Acceptor_HEXseq  <- ""
  sd$Donor_HZEI     <- 0
  sd$Acceptor_HZEI  <- 0
  
  ## Get sequence for upstream and downstream 50nt hexplorer integral
  sd$SD_50nt_up <- ""
  sd$SD_50nt_down <- ""
  sd$HEX_50up_SD <- 0
  sd$HEX_50down_SD <- 0
  
  sd$nProbe <- cur_tab$nProbes[match(sd$IDS,cur_tab$IDS)]
  sd$name <- file_name
  
  anno2[[file_name]] <- sd
}

## Save file
save(anno2, file="anno2")

anno <- list()

## Test distribution of MRIG and aligns in the single BAM files for common treshold
for(i in 1:length(anno2)){ 
  
  ## Create gapsites object
  cur_tab <- anno2[[i]]
  file_name <- names(anno2)[[i]]
  
  cur_tab$gene_id <- as.character(cur_tab$gene_id)
  
  max_transc_expr <- data.frame(tapply(cur_tab$pos_rpmg ,cur_tab$gene_id, max))
  max_transc_expr[2] <- row.names(max_transc_expr)
  names(max_transc_expr) <- c("max_rpmg","gene_id")
  
  cur_tab$MRIG <- max_transc_expr$max_rpmg[match(cur_tab$gene_id, max_transc_expr$gene_id)]
  cur_tab$MRIG <- as.numeric(cur_tab$MRIG)
  
  max_transc_expr <- data.frame(tapply(cur_tab$pos_aln,cur_tab$gene_id, max))
  max_transc_expr[2] <- row.names(max_transc_expr)
  names(max_transc_expr) <- c("max_nAligns","gene_id")
  
  cur_tab$MRIG_aln <- max_transc_expr$max_nAligns[match(cur_tab$gene_id, max_transc_expr$gene_id)]
  cur_tab$MRIG_aln <- as.numeric(cur_tab$MRIG_aln)
  
  
  #Sum aln pro gen
  max_transc_expr <- data.frame(tapply(cur_tab$pos_aln,cur_tab$gene_id, sum))
  max_transc_expr[2] <- row.names(max_transc_expr)
  names(max_transc_expr) <- c("max_nAligns","gene_id")
  
  cur_tab$sum_aln <- max_transc_expr$max_nAligns[match(cur_tab$gene_id, max_transc_expr$gene_id)]
  cur_tab$sum_aln <- as.numeric(cur_tab$sum_aln)
  
  
  cur_tab$file_name_ID <- i
  anno[[file_name]] <- cur_tab
  
}

dt <- rbindlist(anno)

dt_test <- data.frame(file_ID= dt$file_name_ID, MRIG= dt$MRIG, MRIG_aln = dt$MRIG_aln, sum_aln= dt$sum_aln ,gene_name=dt$gene_name)
dt_test <- unique(dt_test)

## Get number of genes with MRIG > 1
dt_test_MRIG1 <- dt_test[dt_test$MRIG > 1, ]
length(unique(dt_test_MRIG1$gene_name))

write.csv2(dt_test, file="Gene_per_sample_MRIG_MRIGaln.csv", row.names=F)
getwd()


## Create SD table "fib"
dt <- rbindlist(anno)
fib <- dt

dt$ID <- paste(dt$seqid, dt$gene_id, dt$norm_lend)
fib$ID <- paste(fib$seqid, fib$gene_id, fib$norm_lend)
fib$pos_rpmg <- NULL
fib$pos_aln <- NULL
fib$MRIG_SD <- NULL
fib$GNR_SD <- NULL
fib$name <- NULL
fib$MRIG <- NULL
fib$MRIG_aln <- NULL
fib$file_name_ID <- NULL
fib$sum_aln <- NULL
fib <- unique(fib)

## Insert info in fib table

## nProbes the SD is found
value <- data.frame(tapply(dt$nProbe,dt$ID,sum))
value$name <- row.names(value)
fib$nProbes <- value[match(fib$ID, value$name),1]

## GNR and MRIG per sample
value <- data.frame(tapply(dt$GNR_SD,dt$ID,mean))
value$name <- row.names(value)
fib$mean_GNR <- value[match(fib$ID, value$name),1]

mrig <- data.frame(gene_id=dt$gene_id, MRIG=dt$MRIG, sample=dt$name)
value <- data.frame(tapply(mrig$MRIG,mrig$gene_id,mean))
value$name <- row.names(value)
fib$mean_MRIG  <- value[match(fib$gene_id, value$name),1]

fib[is.na(fib)] <- 0


## Corrected GNR per group
fib$mean_GNR_corrected <- 0
value <- data.frame(tapply(dt$GNR_SD[dt$MRIG >=1],dt$ID[dt$MRIG >=1],mean))
value$name <- row.names(value)
fib$mean_GNR_corrected <- value[match(fib$ID, value$name),1]
fib$GNR_diff <- fib$mean_GNR_corrected - fib$mean_GNR

## Analyse difference between GNR and corrected GNR
plot(density(fib$mean_GNR_corrected, na.rm=T), xlim=c(0,1))
plot(density(fib$mean_GNR, na.rm=T), xlim=c(0,1))
plot(density(fib$GNR_diff, na.rm=T), xlim=c(-.5,.5))

## Save file
Projectshortcut <- "endo"
save(fib, file=paste(Projectshortcut,("_origin_file")))

load(paste(Projectshortcut,("_origin_file")))



## Add Sequence of SD surrounding
fib$id  <- 1:nrow(fib)
tbl <- fib

ltbl <- tbl
ltbl$begin <- tbl$norm_lend-90
ltbl$end <- tbl$norm_lend+90

ltbl$begin[ltbl$strand =="-"] <- tbl$norm_lend[ltbl$strand =="-"]-90
ltbl$end[ltbl$strand =="-"] <- tbl$norm_lend[ltbl$strand =="-"]+90

## Sort table
ltbl <- ltbl[order(ltbl$seqid, ltbl$begin), ]
seqnames <- sort(unique(ltbl$seqid))
n <- length(seqnames)
dnames <- names(s)
s$MT <- NULL
res_list <- list()

library( seqinr)

## get the sequence of a specific position of a specific chromosome and gene
j <- 1
for(i in 1:n)
{
  lg <- match(seqnames[i], dnames)
  if(!is.na(lg))
  {
    dt <- ltbl[ltbl$seqid== seqnames[i], ]
    test2 <- Views(s[[lg]], start=dt$begin, end=dt$end)
    test <- as.character(test2)
    test2 <- test2@ranges
    test2@NAMES <- as.character(test)
    chr_res <- as.data.frame(test2)
    chr_res$chr <- seqnames[i]
    chr_res$ID <- paste(chr_res$chr, chr_res$start, chr_res$end)
    
    res_list[[j]] <- chr_res
    j <- j + 1
    
  }else{
    cat("Seqname mismatch: ", i, " ", seqnames[i])
  }
}

library(data.table)

res <- rbindlist(res_list)
ltbl$ID <- paste(ltbl$seqid , ltbl$begin, ltbl$end)

ltbl$lseq <- res$names[match(ltbl$ID, res$ID)]

## Minus strand sequecne manipulation
ltbl$lseq[ltbl$strand=="-"] <- gsub("A","t",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("C","g",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("T","a",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("G","c",ltbl$lseq[ltbl$strand=="-"])

ltbl$lseq <- toupper(ltbl$lseq)
ltbl$lseq[ltbl$strand=="-"] <- reverse(ltbl$lseq[ltbl$strand=="-"])

fib$SD_surrounding <-  ltbl$lseq[match(fib$id, ltbl$id)]

table(fib$donor_site)

save(fib, file="endo_processed1")


## Get intrinsic ss strength

## Create a variable for the MaxEnt score
fib$maxent <- 0

## Generate the requiered donor sequence for the MaxEnt score
fib$seqmaxent <- substr(fib$donor_seq, 1,9)

## Load the table with the maxent scores and the associated SD sequences
library(ModCon)

## Calculate SD maxent score
fib$maxent<- calculateMaxEntScanScore(fib$seqmaxent, 5)

## Remove the generated special SD sequence, which was used for MaxEnt generation
fib$seqmaxent <- NULL

## Round the MaxEnt score of the gaps 
fib$maxent <- round(fib$maxent, digits = 1)

## Get the HBS
fib$norm_lhbond <- 0
mtc <- match(fib$donor_seq,hbg$seq)
fib$norm_lhbond <- hbg$hbs[mtc]
fib$norm_lhbond[is.na(fib$norm_lhbond)== T] <- 0

#tidy up the workspace
rm(hbg,m2)


#Save
save(fib, file="SD_Table")



table(fib$donor_site)

## Create sequence columns for HEX ints
fib$Donor_HEXseq     <- as.character(lapply(fib$SD_surrounding, function(x){ substr(x, 84, 104)}))
fib$Donor_HZEI       <- 0

## Create columns with sequence and HZEI information
## Get sequence for upstream and downstream 50nt hexplorer integral
fib$SD_50nt_up <- substr(fib$SD_surrounding, 39,88)
fib$SD_50nt_down <- substr(fib$SD_surrounding, 100,149)

fib$HEX_50up_SD <- 0
fib$HEX_50down_SD <- 0

## Filter SDs with no GT at SD position +1/+2
fib <- fib[fib$donor_site == "GT",]


save(fib, file="SDs with position expression")


## Select SDs which can be found in TSL1 transcripts
SDs <- fib

## Load table with transcript coordinates
load("trans_cord_GRCh38.rda")
trans_cord <- trans_cord_GRCh38

SDs$seqid <- as.character(SDs$seqid)
SDs$strand <- as.character(SDs$strand)
SDs$gene_name <- as.character(SDs$gene_name)

## Check which transcript shows TSL1 Transcript support level 1 from Ensembl version 100
tsl <- read.csv("tsl.txt", stringsAsFactors=FALSE)
library(stringr)
tsl$tsl1 <- str_detect(tsl$Transcript.support.level..TSL., "tsl1")
tsl1 <- tsl$Transcript.stable.ID[tsl$tsl1]

## Add TSL1 info to data frame
trans_cord$tsl1 <- 0
trans_cord$tsl1[trans_cord$Transcript.stable.ID %in% tsl1] <- 1
trans_cord <- trans_cord[trans_cord$tsl1 == 1,]

## Calculate SD coordinate in transcript table
trans_cord$SD_cord <- trans_cord$Exon.region.end..bp.
trans_cord$SD_cord[trans_cord$Strand == -1] <- trans_cord$Exon.region.start..bp.[trans_cord$Strand == -1]

## Get maximal exon rank per transcript
tans_max_ranks <- data.frame(tapply(trans_cord$Exon.rank.in.transcript,trans_cord$Transcript.stable.ID,  max))
tans_max_ranks$name <- row.names(tans_max_ranks)
trans_cord$tmr <- tans_max_ranks$tapply.trans_cord.Exon.rank.in.transcript..trans_cord.Transcript.stable.ID..[match(trans_cord$Transcript.stable.ID,tans_max_ranks$name )]

tsl1_sds <- trans_cord[trans_cord$tmr != trans_cord$Exon.rank.in.transcript,]
tsl1_sds$SD_cord <- paste(tsl1_sds$SD_cord, tsl1_sds$Chromosome.scaffold.name, tsl1_sds$Gene.stable.ID )
tsl1_sds <- unique(tsl1_sds$SD_cord)

## Only select endothelial SDs from TSL1 tanscripts
SDs$SDID <- paste(SDs$norm_lend, SDs$seqid, SDs$gene_id)
SDs$tsl1 <- 0
SDs$tsl1[SDs$SDID %in% tsl1_sds] <- 1
SDs <- SDs[SDs$tsl1 == 1,]
SDs$SDID <- NULL
SDs$tsl1 <- NULL

## Calculate potential nearest SA coordinate from annotated ss of transcript table
i <- 1
SDs$potential_SA_coordinate <- 0
options(warn=2)

## Per SD
for(i in 1:nrow(SDs)){
  
  SD_coord <- SDs$norm_lend[i]
  chr <- SDs$seqid[i]
  strand <- SDs$strand[i]
  geneName<- SDs$gene_id[i]
  
  ## Get list of annotated exons of the same gene
  gene_exons <- trans_cord[trans_cord$Gene.stable.ID == geneName,]
  
  geneSD <- gene_exons$Exon.region.end..bp.
  if(strand == "-")   geneSD <- gene_exons$Exon.region.start..bp.
  names(geneSD) <- rep("SD", length(geneSD))
  
  geneSA <- gene_exons$Exon.region.start..bp.
  if(strand == "-")   geneSA <- gene_exons$Exon.region.end..bp.
  names(geneSA) <- rep("SA", length(geneSA))
  
  SD_coordnamed <- SD_coord
  names(SD_coordnamed) <- "SD_coord"
  
  if(!SD_coord %in% geneSD){
    
    genesplicesites <- c(geneSD,geneSA,SD_coordnamed)
    genesplicesites <- sort(genesplicesites)
    
    SDposinvector <- which(genesplicesites== SD_coord)
    downvec <- 1
    upvec <- 1
    
    if(SDposinvector != 1) while(((SDposinvector-downvec)>=2)&(names(genesplicesites)[SDposinvector-downvec]!="SA")) downvec <- downvec+1
    while((names(genesplicesites)[SDposinvector+upvec]!="SA")&((SDposinvector+upvec)<=length(genesplicesites))) upvec <- upvec+1
    
    selected <- genesplicesites[(SDposinvector-downvec):(SDposinvector+upvec)]
    selected <- selected[!is.na(selected)]
    selected_dist <- abs(selected-SD_coord)
    selected_dist <- selected_dist[selected_dist != 0]
    
    if(length(selected_dist) != 0){ 
      best_select <- which(selected_dist == min(selected_dist))
      selected <- selected[best_select]
      selected <- unique(selected)
      if(strand == "+") gene_exons <- gene_exons[gene_exons$Exon.region.end..bp. == selected,]
      if(strand == "-") gene_exons <- gene_exons[gene_exons$Exon.region.start..bp. == selected,]
    }
    
  }else{
    
    if(strand == "+") gene_exons <- gene_exons[gene_exons$Exon.region.end..bp. == SD_coord,]
    if(strand == "-") gene_exons <- gene_exons[gene_exons$Exon.region.start..bp. == SD_coord,]
    
  }
  
  gene_exons <- gene_exons[gene_exons$Exon.rank.in.transcript != 1,]
  
  if( nrow(gene_exons)>0 ){
    
    relevant_SAs <- gene_exons$Exon.region.start..bp.
    if(strand == "-") relevant_SAs <- gene_exons$Exon.region.end..bp.
    relevant_SAs <- unique(relevant_SAs)
    
    sa_df_dist <- abs(SD_coord - relevant_SAs)
    
    SDs$potential_SA_coordinate[i] <- relevant_SAs[which(sa_df_dist == min(sa_df_dist))]
    
  }
  
  i <- i+1
}

options(warn=1)

save(SDs, file="SDs")
getwd()

## Get SA sequence
SDs$SA_seq <- ""

SDs <- SDs[SDs$potential_SA_coordinate != 0,]
SDs$id <- 1:nrow(SDs)
tbl <- SDs

ltbl <- tbl
ltbl$begin <- tbl$potential_SA_coordinate - 20
ltbl$end <- tbl$potential_SA_coordinate + 2 

ltbl$begin[ltbl$strand =="-"] <- tbl$potential_SA_coordinate[ltbl$strand =="-"]-2
ltbl$end[ltbl$strand =="-"] <- tbl$potential_SA_coordinate[ltbl$strand =="-"]+20

# Sort table
ltbl <- ltbl[order(ltbl$seqid, ltbl$begin), ]
seqnames <- sort(unique(ltbl$seqid))
n <- length(seqnames)
dnames <- names(s)
#s$MT <- NULL
res_list <- list()

library( seqinr)

## Get the sequence of a specific position of a specific chromosome and gene
j <- 1
for(i in 1:n)
{
  lg <- match(seqnames[i], dnames)
  if(!is.na(lg))
  {
    dt <- ltbl[ltbl$seqid== seqnames[i], ]
    test2 <- Views(s[[lg]], start=dt$begin, end=dt$end)
    test <- as.character(test2)
    test2 <- test2@ranges
    test2@NAMES <- as.character(test)
    chr_res <- as.data.frame(test2)
    chr_res$chr <- seqnames[i]
    chr_res$ID <- paste(chr_res$chr, chr_res$start, chr_res$end)
    
    res_list[[j]] <- chr_res
    j <- j + 1
    
  }else{
    cat("Seqname mismatch: ", i, " ", seqnames[i])
  }
}

library(data.table)

res <- rbindlist(res_list)
ltbl$ID <- paste(ltbl$seqid , ltbl$begin, ltbl$end)

ltbl$lseq <- res$names[match(ltbl$ID, res$ID)]

## Minus strand sequence manipulation
ltbl$lseq[ltbl$strand=="-"] <- gsub("A","t",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("C","g",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("T","a",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("G","c",ltbl$lseq[ltbl$strand=="-"])

ltbl$lseq <- toupper(ltbl$lseq)
ltbl$lseq[ltbl$strand=="-"] <- reverse(ltbl$lseq[ltbl$strand=="-"])

SDs$SA_seq <-  ltbl$lseq[match(SDs$id, ltbl$id)]

curwd <- getwd()
getwd()

## Calcualte MaxEntScan score per SA
library(ModCon)
SDs$SA_maxent <- 0

SDs_maxent <- unique(SDs$SA_seq)
SDs_maxent <- SDs_maxent[SDs_maxent != ""]

SMS <- data.frame(seq=SDs_maxent, mxe= 0)
SMS$seq <- as.character(SMS$seq)
SMS$mxe <- calculateMaxEntScanScore(SMS$seq, 3)

SDs$SA_maxent <- SMS$mxe[match(SDs$SA_seq, SMS$seq)]

setwd(curwd)

save(SDs, file="SDs2")



#data for SD pair search
SDs <- na.omit(SDs)

#Calculate exon length
SDs$exon_length <- 0

SDs$exon_length <- lapply(1:nrow(SDs), function(x){
  
  if(SDs$potential_SA_coordinate[x] != 0){ 
    groß <-  max(c(SDs$norm_lend[x], SDs$potential_SA_coordinate[x]))
    klein <- min(c(SDs$norm_lend[x], SDs$potential_SA_coordinate[x]))
    
    SDs$exon_length[x] <- groß-klein+1
  }else SDs$exon_length[x] <- 0
})


SDs$exon_length <- as.numeric(SDs$exon_length)

SDs$SA_ID <- paste(SDs$potential_SA_coordinate, SDs$seqid)

sas <- data.frame(SDs$SA_ID, SDs$exon_length)
sas <- unique(sas)

plot(density(sas$SDs.exon_length), xlim=c(0,500))
mean(sas$SDs.exon_length)
median(sas$SDs.exon_length)

save(SDs, file="SDs3")


SDs$`SD coordinate` <- SDs$norm_lend

## Searching for next best GT site within exon
SDs$id <- 1:nrow(SDs)
SDs$seqid <- as.character(SDs$seqid)
SDs <- na.omit(SDs)

SDs$begin <- SDs$potential_SA_coordinate
SDs$end <- SDs$`SD coordinate`

SDs$begin[SDs$strand =="-"] <- SDs$`SD coordinate`[SDs$strand =="-"]
SDs$end[SDs$strand =="-"] <- SDs$potential_SA_coordinate[SDs$strand =="-"]


# Sort table
SDs <- SDs[order(SDs$seqid, SDs$begin), ]
seqnames <- sort(unique(SDs$seqid))
n <- length(seqnames)
dnames <- names(s)
#s$MT <- NULL
res_list <- list()

library( seqinr)

## Get the sequence of a specific position of a specific chromosome and gene
j <- 1
for(i in 1:n)
{
  lg <- match(seqnames[i], dnames)
  if(!is.na(lg))
  {
    dt <- SDs[SDs$seqid== seqnames[i], ]
    test2 <- Views(s[[lg]], start=dt$begin, end=dt$end)
    test <- as.character(test2)
    test2 <- test2@ranges
    test2@NAMES <- as.character(test)
    chr_res <- as.data.frame(test2)
    chr_res$chr <- seqnames[i]
    chr_res$ID <- paste(chr_res$chr, chr_res$start, chr_res$end)
    
    res_list[[j]] <- chr_res
    j <- j + 1
    
  }else{
    cat("Seqname mismatch: ", i, " ", seqnames[i])
  }
}

library(data.table)

res <- rbindlist(res_list)
SDs$ID <- paste(SDs$seqid , SDs$begin, SDs$end)

SDs$lseq <- res$names[match(SDs$ID, res$ID)]

## Minus strand sequecne manipulation
SDs$lseq[SDs$strand=="-"] <- gsub("A","t",SDs$lseq[SDs$strand=="-"])
SDs$lseq[SDs$strand=="-"] <- gsub("C","g",SDs$lseq[SDs$strand=="-"])
SDs$lseq[SDs$strand=="-"] <- gsub("T","a",SDs$lseq[SDs$strand=="-"])
SDs$lseq[SDs$strand=="-"] <- gsub("G","c",SDs$lseq[SDs$strand=="-"])

SDs$lseq <- toupper(SDs$lseq)
SDs$lseq[SDs$strand=="-"] <- reverse(SDs$lseq[SDs$strand=="-"])

SDs$exon_seq <- SDs$lseq
SDs$lseq <- NULL
SDs$begin <- NULL
SDs$end <- NULL
SDs$ID <- NULL




## Get intronic nucleotide sequence
SDs$begin <- SDs$`SD coordinate`
SDs$end <- SDs$`SD coordinate`+100

SDs$begin[SDs$strand =="-"] <- SDs$`SD coordinate`[SDs$strand =="-"]-100
SDs$end[SDs$strand =="-"] <- SDs$`SD coordinate`[SDs$strand =="-"] 


SDs <- SDs[order(SDs$seqid, SDs$begin), ]
seqnames <- sort(unique(SDs$seqid))
n <- length(seqnames)
dnames <- names(s)
#s$MT <- NULL
res_list <- list()

library( seqinr)

## Get the sequence of a specific position of a specific chromosome and gene
j <- 1
for(i in 1:n)
{
  lg <- match(seqnames[i], dnames)
  if(!is.na(lg))
  {
    dt <- SDs[SDs$seqid== seqnames[i], ]
    test2 <- Views(s[[lg]], start=dt$begin, end=dt$end)
    test <- as.character(test2)
    test2 <- test2@ranges
    test2@NAMES <- as.character(test)
    chr_res <- as.data.frame(test2)
    chr_res$chr <- seqnames[i]
    chr_res$ID <- paste(chr_res$chr, chr_res$start, chr_res$end)
    
    res_list[[j]] <- chr_res
    j <- j + 1
    
  }else{
    cat("Seqname mismatch: ", i, " ", seqnames[i])
  }
}

library(data.table)

res <- rbindlist(res_list)
SDs$ID <- paste(SDs$seqid , SDs$begin, SDs$end)

SDs$lseq <- res$names[match(SDs$ID, res$ID)]

## Minus strand sequence manipulation
SDs$lseq[SDs$strand=="-"] <- gsub("A","t",SDs$lseq[SDs$strand=="-"])
SDs$lseq[SDs$strand=="-"] <- gsub("C","g",SDs$lseq[SDs$strand=="-"])
SDs$lseq[SDs$strand=="-"] <- gsub("T","a",SDs$lseq[SDs$strand=="-"])
SDs$lseq[SDs$strand=="-"] <- gsub("G","c",SDs$lseq[SDs$strand=="-"])

SDs$lseq <- toupper(SDs$lseq)
SDs$lseq[SDs$strand=="-"] <- reverse(SDs$lseq[SDs$strand=="-"])

SDs$intro_seq <- SDs$lseq
SDs$lseq <- NULL
SDs$begin <- NULL
SDs$end <- NULL
SDs$ID <- NULL


## Search for exonic GTs
SDs$alt_gt_pos <- 0
SDs$alt_gt_hbs <- 0
SDs$alt_gt_seq <- ""

res<-data.frame()

save(SDs, file="SDs3")


res <- list()

library(parallel)
no_cores <- detectCores()-1
clust <- makeCluster(no_cores)
clusterExport(clust, list("SDs"), envir=environment())

## Search GT sites and their position within exon sequence
res <- parLapply(clust, SDs$exon_seq, function(x){
  
  after <- x
  ag_pos_after <- gregexpr(pattern ='GT',after)[[1]]
  ag_pos_after <- ag_pos_after[(ag_pos_after > 12) & (ag_pos_after <= nchar(after)-16)]
  
  #If there are GT sites left
  if(length(ag_pos_after)!= 0){
    
    gt_seq <- lapply(ag_pos_after, function(x){substr(after,x-3,x+7)})
    gt_seq <- as.character(gt_seq)
    
    gt_seq
    
  }else{ "" }
  
})

res_pos <- list()

## Get GT site position
res_pos <- parLapply(clust, SDs$exon_seq, function(x){
  
  after <- x
  ag_pos_after <- gregexpr(pattern ='GT',after)[[1]]
  ag_pos_after <- ag_pos_after[(ag_pos_after > 12) & (ag_pos_after <= nchar(after)-16)]
  
  #If there are GT sites left
  if(length(ag_pos_after)!= 0){
    
    gt_seq <- lapply(ag_pos_after, function(x){substr(after,x-3,x+7)})
    gt_seq <- as.character(gt_seq)
    
    ag_pos_after
    
  }else{ 0 }
  
})


save(SDs,res,res_pos, file="SDs3 and additional info")


res_list <- list()

## Generat list of 5'ss pairs
res_list <- lapply(1:nrow(SDs), function(x){
  
  number_of_GTs <- length(res[[x]])
  output_df <- SDs[rep(x, number_of_GTs),]
  
  output_df$alt_gt_pos[1:nrow(output_df)] <- res_pos[[x]]
  output_df$alt_gt_seq[1:nrow(output_df)] <- res[[x]]
  
  output_df
})



library(data.table)
SDs <- rbindlist(res_list)


SDs$alt_gt_hbs <- hbg$hbs[match(SDs$alt_gt_seq, hbg$seq)]
save(SDs, file="SDs4")


## Calcualte metric differences between both SDs
SDs$HBS_diff <- SDs$alt_gt_hbs-SDs$norm_lhbond
SDs$dist_to_SD <- SDs$alt_gt_pos - SDs$exon_length
SDs <- na.omit(SDs)

save(SDs, file="SDs4")

load("SDs4")

## Quickly plot SD pair metrics
plot(density(SDs$exon_length), xlab="Exon length", main ="", xlim=c(0,500))
plot(density(SDs$dist_to_SD), xlab="Distance to SD", main ="", xlim=c(-300,0))
plot(density(SDs$HBS_diff), xlab="HBS diff", main ="", xlim=c(-20,20))
plot(density(SDs$mean_GNR), xlab="GNR per SD", main ="", xlim=c(0,1))


table(SDs$mean_GNR > 0.2)
702261/(nrow(SDs))





## Load SD usage information, saved in line 681
load("SDs with position expression")
sd <- fib
rm(fib)
library(ModCon)
load("SDs4")


## enter information about GT site usage
SDs$pos_sum_mean_GNR_cryptic_GT <- 0
SDs$cryptic_GT_coordinate <- SDs$`SD coordinate` + SDs$dist_to_SD -1
SDs$cryptic_GT_coordinate[SDs$strand == "-"] <- SDs$`SD coordinate` [SDs$strand == "-"] - SDs$dist_to_SD[SDs$strand == "-"] +1
SDs$cryptic_GT_ID  <- paste(SDs$gene_id , SDs$cryptic_GT_coordinate  )

SDs$pos_sum_mean_GNR_cryptic_GT <- sd$mean_GNR[match(SDs$cryptic_GT_ID, sd$IDS)]
SDs[is.na(SDs)] <- 0
SDs$GNR_diff <- SDs$pos_sum_mean_GNR_cryptic_GT / SDs$mean_GNR

## check GT site usage
table(SDs$pos_sum_mean_GNR_cryptic_GT > 0)

## Save SD pair table
save(SDs, file="SDs4")

## Rename columns for clarity
SDs$`SD coordinate`  <-  SDs$norm_lend
SDs$norm_lend <-NULL
SDs$test <- NULL
SDs$pos_sum_mean_GNR <-  SDs$mean_GNR
SDs$mean_GNR <- NULL
SDs$mean_GNR_corrected <- NULL
SDs$id <- NULL
SDs$SA_ID <- NULL
SDs$cryptic_GT_ID <- NULL
SDs$SD_surrounding <- NULL
SDs$nProbe <- NULL
sd <- SDs
save(sd, file="rearranged SD pair table")

## Load reference genome fasta
load("chromosome_sequences_hg38")


## Access nucleotide sequence segments a,b,c upstream of cryptic GT, inbetween both Gt sites and downstream of the used GT
library(spliceSites)


## Add sequence segment a, 60nt (upstream of GT site)
sd$seq_a <- ""
sd$id <- 1:nrow(sd)
tbl <- sd
ltbl <- tbl
ltbl$begin <- tbl$cryptic_GT_coordinate - 57
ltbl$end <- tbl$cryptic_GT_coordinate+2

ltbl$begin[ltbl$strand =="-"] <- tbl$cryptic_GT_coordinate[ltbl$strand =="-"]-2
ltbl$end[ltbl$strand =="-"] <- tbl$cryptic_GT_coordinate[ltbl$strand =="-"]+57

## Sort table
ltbl <- ltbl[order(ltbl$seqid, ltbl$begin), ]
seqnames <- sort(unique(ltbl$seqid))
n <- length(seqnames)
dnames <- names(s)
res_list <- list()

library( seqinr)

## get the sequence of a specific position of a specific chromosome and gene
j <- 1
for(i in 1:n)
{
  lg <- match(seqnames[i], dnames)
  if(!is.na(lg))
  {
    dt <- ltbl[ltbl$seqid== seqnames[i], ]
    test2 <- Views(s[[lg]], start=dt$begin, end=dt$end)
    test <- as.character(test2)
    test2 <- test2@ranges
    test2@NAMES <- as.character(test)
    chr_res <- as.data.frame(test2)
    chr_res$chr <- seqnames[i]
    chr_res$ID <- paste(chr_res$chr, chr_res$start, chr_res$end)
    
    res_list[[j]] <- chr_res
    j <- j + 1
    
  }else{
    cat("Seqname mismatch: ", i, " ", seqnames[i])
  }
}

library(data.table)

res <- rbindlist(res_list)
ltbl$ID <- paste(ltbl$seqid , ltbl$begin, ltbl$end)

ltbl$lseq <- res$names[match(ltbl$ID, res$ID)]

## Minus strand sequence manipulation
ltbl$lseq[ltbl$strand=="-"] <- gsub("A","t",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("C","g",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("T","a",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("G","c",ltbl$lseq[ltbl$strand=="-"])

ltbl$lseq <- toupper(ltbl$lseq)
ltbl$lseq[ltbl$strand=="-"] <- reverse(ltbl$lseq[ltbl$strand=="-"])
sd$seq_a <-  ltbl$lseq[match(sd$id, ltbl$id)]
sd$seq_a_50nt <- sd$seq_a
sd$seq_a <- ""
sd$id <- 1:nrow(sd)


## Add sequence segment a, 50nt (upstream of GT site)
tbl <- sd
ltbl <- tbl
ltbl$begin <- tbl$cryptic_GT_coordinate - 47
ltbl$end <- tbl$cryptic_GT_coordinate+2

ltbl$begin[ltbl$strand =="-"] <- tbl$cryptic_GT_coordinate[ltbl$strand =="-"]-2
ltbl$end[ltbl$strand =="-"] <- tbl$cryptic_GT_coordinate[ltbl$strand =="-"]+47

## Sort table
ltbl <- ltbl[order(ltbl$seqid, ltbl$begin), ]
seqnames <- sort(unique(ltbl$seqid))
n <- length(seqnames)
dnames <- names(s)
res_list <- list()

j <- 1
for(i in 1:n)
{
  lg <- match(seqnames[i], dnames)
  if(!is.na(lg))
  {
    dt <- ltbl[ltbl$seqid== seqnames[i], ]
    test2 <- Views(s[[lg]], start=dt$begin, end=dt$end)
    test <- as.character(test2)
    test2 <- test2@ranges
    test2@NAMES <- as.character(test)
    chr_res <- as.data.frame(test2)
    chr_res$chr <- seqnames[i]
    chr_res$ID <- paste(chr_res$chr, chr_res$start, chr_res$end)
    
    res_list[[j]] <- chr_res
    j <- j + 1
    
  }else{
    cat("Seqname mismatch: ", i, " ", seqnames[i])
  }
}

library(data.table)

res <- rbindlist(res_list)
ltbl$ID <- paste(ltbl$seqid , ltbl$begin, ltbl$end)

ltbl$lseq <- res$names[match(ltbl$ID, res$ID)]

## Minus strand sequence manipulation
ltbl$lseq[ltbl$strand=="-"] <- gsub("A","t",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("C","g",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("T","a",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("G","c",ltbl$lseq[ltbl$strand=="-"])

ltbl$lseq <- toupper(ltbl$lseq)
ltbl$lseq[ltbl$strand=="-"] <- reverse(ltbl$lseq[ltbl$strand=="-"])
sd$seq_a <-  ltbl$lseq[match(sd$id, ltbl$id)]
sd$seq_a_40nt <- sd$seq_a
sd$seq_a <- NULL


## Add sequence segment b (between GT site and annotatged 5'ss)

sd$seq_b <- ""
sd$id <- 1:nrow(sd)
tbl <- sd

#Sequence to the left, generate a new dataset with following informations of the given table
ltbl <- tbl
ltbl$begin <- tbl$cryptic_GT_coordinate+4
ltbl$end <- tbl$`SD coordinate`+2

ltbl$begin[ltbl$strand =="-"] <- tbl$`SD coordinate`[ltbl$strand =="-"]-2
ltbl$end[ltbl$strand =="-"] <- tbl$cryptic_GT_coordinate[ltbl$strand =="-"]-4

# Sort table
ltbl <- ltbl[order(ltbl$seqid, ltbl$begin), ]
seqnames <- sort(unique(ltbl$seqid))
n <- length(seqnames)
dnames <- names(s)
res_list <- list()

j <- 1
for(i in 1:n)
{
  lg <- match(seqnames[i], dnames)
  if(!is.na(lg))
  {
    dt <- ltbl[ltbl$seqid== seqnames[i], ]
    test2 <- Views(s[[lg]], start=dt$begin, end=dt$end)
    test <- as.character(test2)
    test2 <- test2@ranges
    test2@NAMES <- as.character(test)
    chr_res <- as.data.frame(test2)
    chr_res$chr <- seqnames[i]
    chr_res$ID <- paste(chr_res$chr, chr_res$start, chr_res$end)
    
    res_list[[j]] <- chr_res
    j <- j + 1
    
  }else{
    cat("Seqname mismatch: ", i, " ", seqnames[i])
  }
}

library(data.table)

res <- rbindlist(res_list)
ltbl$ID <- paste(ltbl$seqid , ltbl$begin, ltbl$end)

ltbl$lseq <- res$names[match(ltbl$ID, res$ID)]

## Minus strand sequence manipulation
ltbl$lseq[ltbl$strand=="-"] <- gsub("A","t",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("C","g",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("T","a",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("G","c",ltbl$lseq[ltbl$strand=="-"])

ltbl$lseq <- toupper(ltbl$lseq)
ltbl$lseq[ltbl$strand=="-"] <- reverse(ltbl$lseq[ltbl$strand=="-"])

sd$seq_b <-  ltbl$lseq[match(sd$id, ltbl$id)]
sd$seq_b_len <- nchar(sd$seq_b)-10

## Plot sequence b length
plot(density(sd$seq_b_len))


## Get intron sequence
sd$intro_seq <- ""
sd$id <- 1:nrow(sd)
tbl <- sd

ltbl <- tbl
ltbl$begin <- tbl$`SD coordinate`+1
ltbl$end <- tbl$`SD coordinate`+100

ltbl$begin[ltbl$strand =="-"] <- tbl$`SD coordinate`[ltbl$strand =="-"]-100
ltbl$end[ltbl$strand =="-"] <- tbl$`SD coordinate`[ltbl$strand =="-"]-1
ltbl <- ltbl[order(ltbl$seqid, ltbl$begin), ]
seqnames <- sort(unique(ltbl$seqid))
n <- length(seqnames)
dnames <- names(s)
res_list <- list()

## get the sequence of a specific position of a specific chromosome and gene
j <- 1
for(i in 1:n)
{
  lg <- match(seqnames[i], dnames)
  if(!is.na(lg))
  {
    dt <- ltbl[ltbl$seqid== seqnames[i], ]
    test2 <- Views(s[[lg]], start=dt$begin, end=dt$end)
    test <- as.character(test2)
    test2 <- test2@ranges
    test2@NAMES <- as.character(test)
    chr_res <- as.data.frame(test2)
    chr_res$chr <- seqnames[i]
    chr_res$ID <- paste(chr_res$chr, chr_res$start, chr_res$end)
    
    res_list[[j]] <- chr_res
    j <- j + 1
    
  }else{
    cat("Seqname mismatch: ", i, " ", seqnames[i])
  }
}

res <- rbindlist(res_list)
ltbl$ID <- paste(ltbl$seqid , ltbl$begin, ltbl$end)

ltbl$lseq <- res$names[match(ltbl$ID, res$ID)]

#Minus strand sequence manipulation
ltbl$lseq[ltbl$strand=="-"] <- gsub("A","t",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("C","g",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("T","a",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("G","c",ltbl$lseq[ltbl$strand=="-"])

ltbl$lseq <- toupper(ltbl$lseq)
ltbl$lseq[ltbl$strand=="-"] <- reverse(ltbl$lseq[ltbl$strand=="-"])
sd$intro_seq <-  ltbl$lseq[match(sd$id, ltbl$id)]



## Definition of B1 and B2 and B_50nt and B_40nt

## B >= 100 Nt
sd$seq_b1_50nt <- ""
sd$seq_b2_50nt <- ""
sd$seq_b_50nt <- ""

## If length of sequence segment b > 100nt, split b into b1 (for SSHW of GT site) and b2 (for SSHW of 5'ss) which are 60nt long
sd$seq_b1_50nt[sd$seq_b_len >= 100 ] <- substr(sd$seq_b[sd$seq_b_len >= 100 ], 1, 60)
sd$seq_b2_50nt[sd$seq_b_len >= 100 ] <- substr(sd$seq_b[sd$seq_b_len >= 100 ],(sd$seq_b_len[sd$seq_b_len >= 100 ]+10-59) , (sd$seq_b_len[sd$seq_b_len >= 100 ]+10))


## 50 Nt. <= B < 100 Nt
sd$seq_b1_50nt[sd$seq_b_len >= 50] <- substr(sd$seq_b[sd$seq_b_len >= 50 ], 1, 60)
sd$seq_b2_50nt[sd$seq_b_len >= 50 ] <- substr(sd$seq_b[sd$seq_b_len >= 50 ],(sd$seq_b_len[sd$seq_b_len >= 50 ]+10-59) , (sd$seq_b_len[sd$seq_b_len >= 50 ]+10))


## If length of sequence segment b > 6nt and b < 50, use the same b sequence for SSHW of GT and 5'ss
sd$seq_b_50nt[sd$seq_b_len >= 6 & sd$seq_b_len < 50] <- sd$seq_b[sd$seq_b_len >= 6 & sd$seq_b_len < 50]

## Save sd pair table
save(sd, file="rearranged SD pair table")


## Delete pairs with distance smaller 6nt
sd <- sd[sd$seq_b_len > 5, ]

## Remove any cases where sequence was not correctly determind
sd <- sd[nchar(sd$seq_a_40nt)>2 & nchar(sd$seq_a_50nt)>2 & nchar(sd$seq_b)>2 & (nchar(sd$seq_b1_50nt)>2 | sd$seq_b1_50nt == "") & (nchar(sd$seq_b2_50nt)>2 | sd$seq_b2_50nt == ""),]
sd <- na.omit(sd)



## Do the same for 40nt long b subsequences

# B >= 80 Nt
sd$seq_b1_40nt <- ""
sd$seq_b2_40nt <- ""
sd$seq_b_40nt <- ""

sd$seq_b1_40nt[sd$seq_b_len >= 80 ] <- substr(sd$seq_b[sd$seq_b_len >= 80 ], 1, 50)
sd$seq_b2_40nt[sd$seq_b_len >= 80 ] <- substr(sd$seq_b[sd$seq_b_len >= 80 ],(sd$seq_b_len[sd$seq_b_len >= 80 ]+10-49) , (sd$seq_b_len[sd$seq_b_len >= 80 ]+10))


## 50 Nt. <= B < 100 Nt
sd$seq_b1_40nt[sd$seq_b_len >= 40] <- substr(sd$seq_b[sd$seq_b_len >= 40 ], 1, 50)
sd$seq_b2_40nt[sd$seq_b_len >= 40 ] <- substr(sd$seq_b[sd$seq_b_len >= 40 ],(sd$seq_b_len[sd$seq_b_len >= 40 ]+10-49) , (sd$seq_b_len[sd$seq_b_len >= 40 ]+10))


## 6 < B < 50 Nt.: 
sd$seq_b_40nt[sd$seq_b_len >= 6 & sd$seq_b_len < 40] <- sd$seq_b[sd$seq_b_len >= 6 & sd$seq_b_len < 40]



## Delete pairs with distance smaller 6nt
sd <- sd[sd$seq_b_len > 5, ]

sd <- sd[nchar(sd$seq_a_40nt)>2 &  (nchar(sd$seq_b1_40nt)>2 | sd$seq_b1_40nt == "") &  (nchar(sd$seq_b2_40nt)>2 | sd$seq_b2_40nt == ""),]
sd <- na.omit(sd)




## Add sequence segment c
sd$seq_c <- ""
sd$id <- 1:nrow(sd)
tbl <- sd


ltbl <- tbl
ltbl$begin <- tbl$`SD coordinate`+4
ltbl$end <- tbl$`SD coordinate`+63

ltbl$begin[ltbl$strand =="-"] <- tbl$`SD coordinate`[ltbl$strand =="-"]-63
ltbl$end[ltbl$strand =="-"] <- tbl$`SD coordinate`[ltbl$strand =="-"]-4

## Sort table
ltbl <- ltbl[order(ltbl$seqid, ltbl$begin), ]
seqnames <- sort(unique(ltbl$seqid))
n <- length(seqnames)
dnames <- names(s)
res_list <- list()

## Get the sequence of a specific position of a specific chromosome and gene
j <- 1
for(i in 1:n)
{
  lg <- match(seqnames[i], dnames)
  if(!is.na(lg))
  {
    dt <- ltbl[ltbl$seqid== seqnames[i], ]
    test2 <- Views(s[[lg]], start=dt$begin, end=dt$end)
    test <- as.character(test2)
    test2 <- test2@ranges
    test2@NAMES <- as.character(test)
    chr_res <- as.data.frame(test2)
    chr_res$chr <- seqnames[i]
    chr_res$ID <- paste(chr_res$chr, chr_res$start, chr_res$end)
    
    res_list[[j]] <- chr_res
    j <- j + 1
    
  }else{
    cat("Seqname mismatch: ", i, " ", seqnames[i])
  }
}

library(data.table)

res <- rbindlist(res_list)
ltbl$ID <- paste(ltbl$seqid , ltbl$begin, ltbl$end)

ltbl$lseq <- res$names[match(ltbl$ID, res$ID)]

#Minus strand sequence manipulation
ltbl$lseq[ltbl$strand=="-"] <- gsub("A","t",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("C","g",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("T","a",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("G","c",ltbl$lseq[ltbl$strand=="-"])

ltbl$lseq <- toupper(ltbl$lseq)
ltbl$lseq[ltbl$strand=="-"] <- reverse(ltbl$lseq[ltbl$strand=="-"])

sd$seq_c <-  ltbl$lseq[match(sd$id, ltbl$id)]
sd$seq_c_50nt <- sd$seq_c
sd$seq_c_40nt <- substr(sd$seq_c_50nt, 1, 50)

sd$seq_c <- NULL
sd <- sd[nchar(sd$seq_c_50nt)>2 & nchar(sd$seq_c_40nt)>2 ,]
save(sd, file="rearranged SD pair table")

## Check creation of sequence segment b
sum(nchar(sd$seq_b)< 11)
sum((nchar(sd$seq_b_40nt)< 11)&((sd$seq_b_40nt != "")))
sum((nchar(sd$seq_b_50nt)< 11)&((sd$seq_b_50nt != "")))
sum((nchar(sd$seq_b1_40nt)< 11)&((sd$seq_b1_40nt != "")))
sum((nchar(sd$seq_b1_50nt)< 11)&((sd$seq_b1_50nt != "")))
sum((nchar(sd$seq_b2_40nt)< 11)&((sd$seq_b2_40nt != "")))
sum((nchar(sd$seq_b2_50nt)< 11)&((sd$seq_b2_50nt != "")))




## All SD pair sequence segments ready for HZEI calculation

## HZEI/nt
## Overlap Vector Function
splitWithOverlap <- function(vec, seg.length, overlap){
  starts = seq(1, length(vec), by=seg.length-overlap)
  ends   = starts + seg.length - 1
  ends[ends > length(vec)] = length(vec)
  
  lapply(1:length(starts), function(i) vec[starts[i]:ends[i]])
}

## HEXplorer plotting Function
calculate_HEXint <- function(seq){
  
  out <- splitWithOverlap(strsplit(toupper(seq), "")[[1]],6,5)
  out <- lapply(out, function(x) paste(x,collapse=""))
  out <- out[nchar(unlist(out))==6]
  out <- as.character(hex$value[match(out,hex$seq)])
  out9 <- splitWithOverlap(out,6,5)
  out9 <- out9[lapply(out9,length)==6]
  out9 <- lapply(out9, function(x){mean(as.numeric(x))})
  out9 <- sum(unlist(out9))
  
}

library(ModCon)

## make it parallel
## Setup cluster
no_cores <- detectCores()-1
clust <- makeCluster(no_cores)

## Export needed data to clusters
clusterExport(clust, list("sd", "calculate_HEXint", "splitWithOverlap","hex"), envir=environment())

## Calcualte total HZEI for every sequence segment
sd$seq_a_50nt_HZEIint    <-  as.numeric(parLapply(clust, sd$seq_a_50nt, calculate_HEXint))
sd$seq_a_40nt_HZEIint    <-  as.numeric(parLapply(clust, sd$seq_a_40nt, calculate_HEXint))

sd$seq_b_50nt_HZEIint[sd$seq_b_50nt != ""]    <-  as.numeric(parLapply(clust, sd$seq_b_50nt[sd$seq_b_50nt != ""], calculate_HEXint))
sd$seq_b_40nt_HZEIint[sd$seq_b_40nt != ""]    <-  as.numeric(parLapply(clust, sd$seq_b_40nt[sd$seq_b_40nt != ""], calculate_HEXint))

sd$seq_b1_50nt_HZEIint[sd$seq_b1_50nt != ""]    <-  as.numeric(parLapply(clust, sd$seq_b1_50nt[sd$seq_b1_50nt != ""], calculate_HEXint))
sd$seq_b2_50nt_HZEIint[sd$seq_b2_50nt != ""]    <-  as.numeric(parLapply(clust, sd$seq_b2_50nt[sd$seq_b2_50nt != ""], calculate_HEXint))

sd$seq_b1_40nt_HZEIint[sd$seq_b1_40nt != ""]    <-  as.numeric(parLapply(clust, sd$seq_b1_40nt[sd$seq_b1_40nt != ""], calculate_HEXint))
sd$seq_b2_40nt_HZEIint[sd$seq_b2_40nt != ""]    <-  as.numeric(parLapply(clust, sd$seq_b2_40nt[sd$seq_b2_40nt != ""], calculate_HEXint))

sd$seq_c_50nt_HZEIint    <-  as.numeric(parLapply(clust, sd$seq_c_50nt, calculate_HEXint))
sd$seq_c_40nt_HZEIint    <-  as.numeric(parLapply(clust, sd$seq_c_40nt, calculate_HEXint))


## Stop the Cluster
stopCluster(clust) 

## Save sd pairs
save(sd, file="SDs_out")

## Create GT id
sd$cryptGT_ID <- paste(sd$IDS, sd$alt_gt_pos)
test$cryptGT_ID <- paste(test$IDS, test$alt_gt_pos)
sd$alt_gt_seq <- test$alt_gt_seq[match(sd$cryptGT_ID, test$cryptGT_ID)]



## calculate SSHW (called SRE support in the following)

sd$crypt_GT_SRE_support_50nt <- 0
sd$SD_SRE_support_50nt <- 0

# 50 Nt. <=  B
sd$crypt_GT_SRE_support_50nt[sd$seq_b_len >= 50] <- sd$seq_a_50nt_HZEIint[sd$seq_b_len >= 50] - sd$seq_b1_50nt_HZEIint[sd$seq_b_len >= 50]
sd$SD_SRE_support_50nt[sd$seq_b_len >= 50] <-  sd$seq_b2_50nt_HZEIint[sd$seq_b_len >= 50] - sd$seq_c_50nt_HZEIint[sd$seq_b_len >= 50]

# 6 < B < 50 Nt.: 
sd$crypt_GT_SRE_support_50nt[sd$seq_b_len < 50] <- sd$seq_a_50nt_HZEIint[sd$seq_b_len < 50] - sd$seq_b_50nt_HZEIint[sd$seq_b_len < 50]
sd$SD_SRE_support_50nt[sd$seq_b_len < 50] <- sd$seq_b_50nt_HZEIint[sd$seq_b_len < 50] - sd$seq_c_50nt_HZEIint[sd$seq_b_len < 50]





sd$crypt_GT_SRE_support_40nt <- 0
sd$SD_SRE_support_40nt <- 0

# 50 Nt. <=  B
sd$crypt_GT_SRE_support_40nt[sd$seq_b_len >= 40] <- sd$seq_a_40nt_HZEIint[sd$seq_b_len >= 40] - sd$seq_b1_40nt_HZEIint[sd$seq_b_len >= 40]
sd$SD_SRE_support_40nt[sd$seq_b_len >= 40] <-  sd$seq_b2_40nt_HZEIint[sd$seq_b_len >= 40] - sd$seq_c_40nt_HZEIint[sd$seq_b_len >= 40]

# 6 < B < 50 Nt.: 
sd$crypt_GT_SRE_support_40nt[sd$seq_b_len < 40] <- sd$seq_a_40nt_HZEIint[sd$seq_b_len < 40] - sd$seq_b_40nt_HZEIint[sd$seq_b_len < 40]
sd$SD_SRE_support_40nt[sd$seq_b_len < 40] <- sd$seq_b_40nt_HZEIint[sd$seq_b_len < 40] - sd$seq_c_40nt_HZEIint[sd$seq_b_len < 40]

sd$donor_seq <- sd$donor_sequence
sd$donor_sequence <- NULL


## Focusing on GT sites
sd <- sd[sd$donor_site == "GT",]

## Check if GT or 5'ss contains U12 donor site or potential TT woble at SD pos +1/+2
library(stringr)
sd$motif_tatcct_sd <- 0
sd$motif_tatcct_sd <- str_detect(sd$donor_seq, "TATCCT")

sd$crypt_GT_seq <- sd$alt_gt_seq
sd$motif_tatcct_cypt <- 0
sd$motif_tatcct_cypt <- str_detect(sd$crypt_GT_seq, "TATCCT")

sd$motif_tatcct <- sd$motif_tatcct_sd + sd$motif_tatcct_cypt
sd$motif_tatcct_sd   <- NULL
sd$motif_tatcct_cypt <- NULL

sd$gtt_search <- 0
sd$gtt_search[substr(sd$donor_seq, 4,6) == "GTT"] <- 1

sd$gtt_search_crypt <- 0
sd$gtt_search_crypt[substr(sd$crypt_GT_seq, 4,6) == "GTT"] <- 1

sd$gtt_motif <- sd$gtt_search + sd$gtt_search_crypt
sd$gtt_search <- NULL
sd$gtt_search_crypt <- NULL
sd <- sd[sd$motif_tatcct == 0 & sd$gtt_motif == 0,]
plot(density(sd$HBS_diff))
save(sd, file="SD pairs")



## Load full table of potential u12 donor sequences, taken from Alioto, T.S. (2007) U12DB: a database of orthologous U12-type spliceosomal introns. Nucleic Acids Res
u12 <- read.table("U12 donors.txt", quote="\"", comment.char="", stringsAsFactors=FALSE)

u12$donor_seq <- substr(u12$V1, 11, 18)
u12 <- unique(u12$donor_seq)
sd$motif_tatcct <- NULL
sd$U12_donor <- 0

## Check occurence of every motif in both SDs
for(i in 1:nrow(sd)){
  
  sd_seq <- sd$donor_seq[i]
  gt_seq <- sd$crypt_GT_seq[i]
  
  if(any(sapply(u12, grepl, c(gt_seq,sd_seq)))) sd$U12_donor[i] <- 1
  
}


sd$test_seq <- paste0(sd$exon_seq, sd$intro_seq)
sd$up_seq <- substr(sd$exon_seq, nchar(sd$exon_seq)-3,nchar(sd$exon_seq))
sd$down_seq <- substr(sd$intro_seq, 3,5)

library(ModCon)

sd$alt_anno_hbs <- 0
sd$alt_anno_seq <- ""

save(sd, file="zwischenspeicher")


## Get alternative donor sites, which overlap the annotated 5'ss and which would distort the SD usage calculations
position <- 4
sd$up_seq1 <- substr(sd$test_seq,( sd$exon_length+2-position-4), (sd$exon_length+2-position+6))
position <- 3
sd$up_seq2 <- substr(sd$test_seq,( sd$exon_length+2-position-4), (sd$exon_length+2-position+6))
position <- 2
sd$up_seq3 <- substr(sd$test_seq,( sd$exon_length+2-position-4), (sd$exon_length+2-position+6))
position <- 1
sd$up_seq4 <- substr(sd$test_seq,( sd$exon_length+2-position-4), (sd$exon_length+2-position+6))


position <- 1
sd$down_seq1 <- substr(sd$test_seq,( sd$exon_length+2+position-3), (sd$exon_length+2+position+7))
position <- 2
sd$down_seq2 <- substr(sd$test_seq,( sd$exon_length+2+position-3), (sd$exon_length+2+position+7))
position <- 3
sd$down_seq3 <- substr(sd$test_seq,( sd$exon_length+2+position-3), (sd$exon_length+2+position+7))
position <- 4
sd$down_seq4 <- substr(sd$test_seq,( sd$exon_length+2+position-3), (sd$exon_length+2+position+7))


## Get hbs for near overlapping 5'ss 
sd$up_seq1_hbs <- hbg$hbs[match(sd$up_seq1, hbg$seq)]
sd$up_seq2_hbs <- hbg$hbs[match(sd$up_seq2, hbg$seq)]
sd$up_seq3_hbs <- hbg$hbs[match(sd$up_seq3, hbg$seq)]
sd$up_seq4_hbs <- hbg$hbs[match(sd$up_seq4, hbg$seq)]

sd$down_seq1_hbs <- hbg$hbs[match(sd$down_seq1, hbg$seq)]
sd$down_seq2_hbs <- hbg$hbs[match(sd$down_seq2, hbg$seq)]
sd$down_seq3_hbs <- hbg$hbs[match(sd$down_seq3, hbg$seq)]
sd$down_seq4_hbs <- hbg$hbs[match(sd$down_seq4, hbg$seq)]

sd[is.na(sd)] <- 0

head(sd)


## make it parallel
## Setup cluster
no_cores <- detectCores()-1
clust <- makeCluster(no_cores)

clusterExport(clust, list("sd"), envir=environment())

## Access the  maximal HBS of overlapping SD sites
sd$alt_anno_hbs <- parLapply(clust, c(1:nrow(sd)), function(x){
  
  max(sd$up_seq1_hbs[x], sd$up_seq2_hbs[x], sd$up_seq3_hbs[x], sd$up_seq4_hbs[x], sd$down_seq1_hbs[x],
                         sd$down_seq2_hbs[x], sd$down_seq3_hbs[x], sd$down_seq4_hbs[x])
  
})

sd$alt_anno_hbs <- as.numeric(sd$alt_anno_hbs)

## save sd pair table
save(sd, file="zwischenspeicher")

sd$norm_lhbond <- as.numeric(sd$norm_lhbond)
sd$hbs_diff_annoGTs <- sd$alt_anno_hbs-sd$norm_lhbond
table(sd$hbs_diff_annoGTs < 0)

## Delete SD pairs, where 5'ss shows overlapping SD of higher HBS, where U12 ss occur or GTT SD sites
sd_save <- sd
sd <- sd[sd$hbs_diff_annoGTs < 0,]
sd <- sd[, -c(67:88)]
sd <- sd[sd$U12_donor == 0 & sd$gtt_motif == 0,]

sd$gtt_motif <- NULL
sd$U12_donor <- NULL
sd$id <- NULL



## Before saving reverse point in hbs column
sd$norm_lhbond <- gsub("\\.",",",sd$norm_lhbond)
save(sd, file="zwischenspeicher2")

## Delete not needed columns
sd$Donor_HEXseq <- NULL
sd$Donor_HZEI <- NULL
sd$Acceptor_HEXseq <- NULL
sd$Acceptor_HZEI <- NULL
sd$HEX_50down_SD  <- NULL
sd$HEX_50up_SD <- NULL

sd$intro_seq_downstream <- sd$intro_seq
sd$intro_seq <- NULL



## Add upstream intronic sequence
load("chromosome_sequences_hg38")
sd$intro_seq_upstream <- ""
sd$id <- 1:nrow(sd)
tbl <- sd

ltbl <- tbl
ltbl$begin <- tbl$potential_SA_coordinate - 70
ltbl$end <- tbl$potential_SA_coordinate

ltbl$begin[ltbl$strand =="-"] <- tbl$potential_SA_coordinate[ltbl$strand =="-"]
ltbl$end[ltbl$strand =="-"] <- tbl$potential_SA_coordinate[ltbl$strand =="-"]+70

# Sort table
ltbl <- ltbl[order(ltbl$seqid, ltbl$begin), ]
seqnames <- sort(unique(ltbl$seqid))
n <- length(seqnames)
dnames <- names(s)
res_list <- list()


j <- 1
for(i in 1:n)
{
  lg <- match(seqnames[i], dnames)
  if(!is.na(lg))
  {
    dt <- ltbl[ltbl$seqid== seqnames[i], ]
    test2 <- Views(s[[lg]], start=dt$begin, end=dt$end)
    test <- as.character(test2)
    test2 <- test2@ranges
    test2@NAMES <- as.character(test)
    chr_res <- as.data.frame(test2)
    chr_res$chr <- seqnames[i]
    chr_res$ID <- paste(chr_res$chr, chr_res$start, chr_res$end)
    
    res_list[[j]] <- chr_res
    j <- j + 1
    
  }else{
    cat("Seqname mismatch: ", i, " ", seqnames[i])
  }
}

library(data.table)

res <- rbindlist(res_list)
ltbl$ID <- paste(ltbl$seqid , ltbl$begin, ltbl$end)

ltbl$lseq <- res$names[match(ltbl$ID, res$ID)]

## Minus strand sequence manipulation
ltbl$lseq[ltbl$strand=="-"] <- gsub("A","t",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("C","g",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("T","a",ltbl$lseq[ltbl$strand=="-"])
ltbl$lseq[ltbl$strand=="-"] <- gsub("G","c",ltbl$lseq[ltbl$strand=="-"])

ltbl$lseq <- toupper(ltbl$lseq)
ltbl$lseq[ltbl$strand=="-"] <- reverse(ltbl$lseq[ltbl$strand=="-"])

sd$intro_seq_upstream <-  ltbl$lseq[match(sd$id, ltbl$id)]

sd$id <- NULL
sd$nProbe <- NULL
sd$SD_surrounding <- NULL
sd$SA_maxent <- gsub("\\.",",",sd$SA_maxent)

save(sd, file="SDs before reconstruction")

Projectshortcut <- 
load("SDs before reconstruction")
Projectshortcut <- "endo"
load(paste(Projectshortcut,("_origin_file")))

SDs <- sd

SDs$IDS <- paste(SDs$gene_id, SDs$`SD coordinate`)
fib$IDS <- paste(fib$gene_id, fib$norm_lend)

fib$mean_GNR <-  fib$mean_GNR_corrected
SDs$pos_sum_mean_GNR <- fib$mean_GNR[match(SDs$IDS, fib$IDS)]

SDs$crypt_GT_seq <- NULL
SDs$GNR_diff <- NULL


SDs$pos_sum_mean_GNR_cryptic_GT <- 0
SDs$cryptic_GT_coordinate <- SDs$`SD coordinate` + SDs$dist_to_SD -1
SDs$cryptic_GT_coordinate[SDs$strand == "-"] <- SDs$`SD coordinate` [SDs$strand == "-"] - SDs$dist_to_SD[SDs$strand == "-"] +1

SDs$cryptic_GT_ID  <- paste(SDs$gene_id , SDs$cryptic_GT_coordinate  )

SDs$pos_sum_mean_GNR_cryptic_GT <- fib$mean_GNR[match(SDs$cryptic_GT_ID, fib$IDS)]
SDs[is.na(SDs)] <- 0

SDs$cryptic_GT_ID <- NULL
SDs$GNR_ratio <- SDs$pos_sum_mean_GNR_cryptic_GT/SDs$pos_sum_mean_GNR

SDs2 <- SDs[,c(1:24,57,25:56)]

SDs2$seqid[SDs2$seqid == "X"] <- 23
SDs2$seqid[SDs2$seqid == "Y"] <- 24

SDs <- SDs2

SDs <- na.omit(SDs)

sum(is.na(SDs$pos_sum_mean_GNR) )
sum(SDs$pos_sum_mean_GNR_cryptic_GT > 0)

SDs <- SDs[SDs$pos_sum_mean_GNR > 0,]

save(SDs, file="SDs final 2")
write.csv2(SDs, file="SD pairs 9.1.csv", row.names = F)
getwd()



names(SDs)

SDs$seq_a_50nt_HZEI_nt <-  SDs$seq_a_50nt_HZEIint/(nchar(SDs$seq_a_50nt)-10 )     
SDs$seq_a_40nt_HZEI_nt <-  SDs$seq_a_40nt_HZEIint/(nchar(SDs$seq_a_40nt)-10 ) 
SDs$seq_b_50nt_HZEI_nt <-  SDs$seq_b_50nt_HZEIint/(nchar(SDs$seq_b_50nt)-10 ) 
SDs$seq_b_40nt_HZEI_nt <-  SDs$seq_b_40nt_HZEIint /(nchar(SDs$seq_b_40nt)-10 )   
SDs$seq_b1_50nt_HZEI_nt <- SDs$seq_b1_50nt_HZEIint /(nchar(SDs$seq_b1_50nt)-10 ) 
SDs$seq_b2_50nt_HZEI_nt <- SDs$seq_b2_50nt_HZEIint/(nchar(SDs$seq_b2_50nt)-10 ) 
SDs$seq_b1_40nt_HZEI_nt <- SDs$seq_b1_40nt_HZEIint/(nchar(SDs$seq_b1_40nt)-10 ) 
SDs$seq_b2_40nt_HZEI_nt <- SDs$seq_b2_40nt_HZEIint/(nchar(SDs$seq_b2_40nt)-10 ) 
SDs$seq_c_50nt_HZEI_nt <-  SDs$seq_c_50nt_HZEIint /(nchar(SDs$seq_c_50nt)-10 ) 
SDs$seq_c_40nt_HZEI_nt <-  SDs$seq_c_40nt_HZEIint/(nchar(SDs$seq_c_40nt)-10 ) 


sd <- SDs

## calculate SRE suppors

sd$crypt_GT_SRE_support_50nt <- 0
sd$SD_SRE_support_50nt <- 0

# 50 Nt. <=  B
sd$crypt_GT_SRE_support_50nt[sd$seq_b_len >= 50] <- sd$seq_a_50nt_HZEI_nt[sd$seq_b_len >= 50] - sd$seq_b1_50nt_HZEI_nt[sd$seq_b_len >= 50]
sd$SD_SRE_support_50nt[sd$seq_b_len >= 50] <-  sd$seq_b2_50nt_HZEI_nt[sd$seq_b_len >= 50] - sd$seq_c_50nt_HZEI_nt[sd$seq_b_len >= 50]

# 6 < B < 50 Nt.: 
sd$crypt_GT_SRE_support_50nt[sd$seq_b_len < 50] <- sd$seq_a_50nt_HZEI_nt[sd$seq_b_len < 50] - sd$seq_b_50nt_HZEI_nt[sd$seq_b_len < 50]
sd$SD_SRE_support_50nt[sd$seq_b_len < 50] <- sd$seq_b_50nt_HZEI_nt[sd$seq_b_len < 50] - sd$seq_c_50nt_HZEI_nt[sd$seq_b_len < 50]





sd$crypt_GT_SRE_support_40nt <- 0
sd$SD_SRE_support_40nt <- 0

# 50 Nt. <=  B
sd$crypt_GT_SRE_support_40nt[sd$seq_b_len >= 40] <- sd$seq_a_40nt_HZEI_nt[sd$seq_b_len >= 40] - sd$seq_b1_40nt_HZEI_nt[sd$seq_b_len >= 40]
sd$SD_SRE_support_40nt[sd$seq_b_len >= 40] <-  sd$seq_b2_40nt_HZEI_nt[sd$seq_b_len >= 40] - sd$seq_c_40nt_HZEI_nt[sd$seq_b_len >= 40]

# 6 < B < 50 Nt.: 
sd$crypt_GT_SRE_support_40nt[sd$seq_b_len < 40] <- sd$seq_a_40nt_HZEI_nt[sd$seq_b_len < 40] - sd$seq_b_40nt_HZEI_nt[sd$seq_b_len < 40]
sd$SD_SRE_support_40nt[sd$seq_b_len < 40] <- sd$seq_b_40nt_HZEI_nt[sd$seq_b_len < 40] - sd$seq_c_40nt_HZEI_nt[sd$seq_b_len < 40]





sd$norm_lhbond <- gsub(",","\\.",sd$norm_lhbond)
sd$norm_lhbond <- as.numeric(sd$norm_lhbond)
sd$maxent <- gsub(",","\\.",sd$maxent)
sd$maxent <- as.numeric(sd$maxent)

## Select sd pairs for classification
sd <- sd[sd$GNR_ratio < 0.01,]
sd <- sd[sd$alt_gt_hbs >= 10,]


## Retrieve the ldiff values from file from above
load(paste0("endo",("_annotated_grs_file_sample_seperated")))
library("data.table")
fib <- rbindlist(anno)
fib$IDS <- paste(fib$gene_id, fib$norm_lend)

sd$ldiff <- fib$ldiff[match(sd$IDS, fib$IDS)]
sd <- sd[sd$ldiff == 0,]

## Create table of all GT sites and 5'ss with their respective HBS and SSHW
res <- data.frame(SRE_support= c(sd$crypt_GT_SRE_support_50nt, sd$SD_SRE_support_50nt),
                  donor=c(rep("crypt GT",nrow(sd)),rep("SD",nrow(sd) )),
                  HBS=c(sd$alt_gt_hbs, sd$norm_lhbond))


res$id <- 1:nrow(res)

## Subset data in train "randomly" and test data sets in 3:1 ratio
train <- res[sample(1:nrow(res),nrow(res)*0.75),]
test  <- res[!res$id %in% train$id,]

## save test and train data
save(test,train, file="endo_ROC_data")

library("ROCit")

test$donor_0 <- 0
test$donor_0[test$donor == "SD"] <- 1


## Model using only SSHW
## First, fit a logistic model
logistic.model_SRE_support <- glm(as.factor(donor)~SRE_support,
                                  data = train,
                                  family = "binomial")

summary(logistic.model_SRE_support)

## Create the score and get its associated class
class <- test$donor_0
test$score <- predict(logistic.model_SRE_support, test, type="response")

## score = log odds
score <- qlogis(test$score)

## Create rocit object
rocit_emp_SRE_support <- rocit(score = score, 
                               class = class, 
                               method = "emp")

## Model using SSHW and HBS
## First, fit a logistic model
logistic.model_SRE_support_HBS <- glm(as.factor(donor)~SRE_support*HBS,
                                      data = train,
                                      family = "binomial")

summary(logistic.model_SRE_support_HBS)

## Create the score and get its associated class
class <- test$donor_0
test$score <- predict(logistic.model_SRE_support_HBS, test, type="response")

# score = log odds
score <- qlogis(test$score)

## Create rocit object
rocit_emp_SRE_support_HBS <- rocit(score = score, 
                                   class = class, 
                                   method = "emp")




## Model using only HBS
logistic.model_HBS <- glm(as.factor(donor)~HBS,
                          data = train,
                          family = "binomial")

summary(logistic.model_HBS)

## Create the score and get its associated class
class <- test$donor_0
test$score <- predict(logistic.model_HBS, test, type="response")

# score = log odds
score <- qlogis(test$score)


## Create rocit object
rocit_emp_HBS <- rocit(score = score, 
                       class = class, 
                       method = "emp")


## Asses the quality of the models
summary(rocit_emp_HBS)
# AUC= 0.9253  
hbs <- data.frame(FPR=rocit_emp_HBS$FPR, TPR=rocit_emp_HBS$TPR)

summary(rocit_emp_SRE_support)
# AUC= 0.8911   
sre <- data.frame(FPR=rocit_emp_SRE_support$FPR, TPR=rocit_emp_SRE_support$TPR)

summary(rocit_emp_SRE_support_HBS)
# AUC=0.9779 
sreHBS <- data.frame(FPR=rocit_emp_SRE_support_HBS$FPR, TPR=rocit_emp_SRE_support_HBS$TPR)

## Generate ROC curve
library(ggplot2)
hbs$parameter <- "HBS"
sre$parameter <- "SSHW"
sreHBS$parameter <- "HBS and SSHW"

vis <- rbind(hbs, sre, sreHBS)

vis$parameter <- factor(vis$parameter, levels=  c("HBS and SSHW", "HBS", "SSHW"))

largesvg <- ggplot(vis, aes(y=TPR, x=FPR, color=parameter))+geom_line(size=1.6)  +
  theme(legend.position="bottom", legend.title = element_blank(), text = element_text(size = 19)  )

smallsvg <- ggplot(vis, aes(y=TPR, x=FPR, color=parameter))+geom_line(size=1.6) + coord_cartesian(xlim=c(0, 0.25), ylim=c(0.75,1)) +
  theme(legend.position="bottom", legend.title = element_blank(), text = element_text(size = 19)  )

largesvg
smallsvg

ggsave(file="ENDO ROC.svg", plot=largesvg, width=5, height=5)
ggsave(file="ENDO ROC small.svg", plot=smallsvg, width=5, height=5)

write.csv2(SDs, file= "Endo SD pairs 9.1.csv", row.names=F)

test$type <- "test"
train$type <- "train"
test$class <- NULL
test$score <- NULL
test$donor_0 <- NULL
out <- rbind(test, train)
write.csv2(out, file= "ROC test and train Endo.csv", row.names=F)
