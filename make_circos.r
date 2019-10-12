#!/usr/bin/env Rscript
library(RCircos)

# usage: Rscript make_circos.r <sv table> <sample name> <gene label table> <cnv data> <out>

# parse args
args = commandArgs(trailingOnly=TRUE)
sv.file <- args[1]
sample.name <- args[2]
gene.label.file <- args[3]
cnv.file <- args[4]
out.file <- args[5]
# TMP <- Sys.getenv("TMP_DIR") 
# tmp.bed = paste0(TMP ,"/" , sample.name, "_bkpts.bed")
tmp.bed = paste0(sample.name, "_bkpts.bed")

# load prereq data
data(UCSC.HG19.Human.CytoBandIdeogram)

# set core parameters
chr.exclude <- NULL;
cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;
tracks.inside <- 10;
tracks.outside <- 5;
RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside);
rcircos.params <- RCircos.Get.Plot.Parameters();
rcircos.params$text.size <- 1
RCircos.Reset.Plot.Parameters(rcircos.params)
rcircos.cyto <- RCircos.Get.Plot.Ideogram();
rcircos.position <- RCircos.Get.Plot.Positions();
RCircos.List.Plot.Parameters()

link.data <- tryCatch(read.table(sv.file, sep = ',', stringsAsFactors = F, header = T), error=function(e) data.frame())
                    
if (nrow(link.data) != 0) {
  
  link.data <- transform(link.data,
                         chromStart = as.numeric(chromStart),
                         chromEnd = as.numeric(chromEnd),
                         chromStart.1 = as.numeric(chromStart.1),
                         chromEnd.1 = as.numeric(chromEnd.1))
  
  # write a bed file of all breakpoints to intersect with gene label table
  bkpts.1 <- link.data[c("Chromosome", "chromStart", "chromEnd")]
  bkpts.2 <- link.data[c("Chromosome.1", "chromStart.1", "chromEnd.1")]
  colnames(bkpts.2) <- colnames(bkpts.1)
  write.table(rbind(bkpts.1, bkpts.2), tmp.bed, sep = '\t', quote = F, col.names = F, row.names = F)
  
  # only keep labels that fall within an event
  print(paste0("bedtools intersect -wb -a ", tmp.bed, " -b ", gene.label.file))
  gene.labels <- system(paste0("bedtools intersect -wb -a ", tmp.bed, " -b ", gene.label.file), intern = T)
  gene.labels <- data.frame(do.call('rbind', strsplit(gene.labels, '\t', fixed=TRUE)), stringsAsFactors = F)
  if (nrow(gene.labels) > 0) {
    gene.labels <- gene.labels[,4:ncol(gene.labels)]
    
    # deduplicate labels
    gene.labels <- gene.labels[!duplicated(gene.labels),]
    colnames(gene.labels) <- c("Chromosome", "chromStart", "chromEnd", "gene")
    gene.labels <- transform(gene.labels,
                             chromStart = as.numeric(chromStart),
                             chromEnd = as.numeric(chromEnd))
  }
  
  # make the plot
  png(file=out.file, height=3000, width=3000, res = 500)
  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot()
  track.num <- 2
  RCircos.Link.Plot(link.data, track.num, TRUE)
  title(sample.name, line=-1)
  
  # label the genes
  if (nrow(gene.labels) > 0) {
    name.col <- 4
    side <- "out"
    track.num <- 1
    RCircos.Gene.Connector.Plot(gene.labels, track.num, side);
    track.num <- 2
    RCircos.Gene.Name.Plot(gene.labels, name.col, track.num, side);
  }
    
  # remove intermediate file
  system(paste0("rm -f ", tmp.bed))
  
} else {
  # make empty plot
  png(file=out.file, height=3000, width=3000, res = 500)
  RCircos.Set.Plot.Area()
  RCircos.Chromosome.Ideogram.Plot()
  title(sample.name, line=-1)
}
                    
# parse cnv data
cnv = tryCatch(read.csv(cnv.file, stringsAsFactors = F), error=function(e) data.frame())
               
if (nrow(cnv) != 0) {
    colnames(cnv) <- c("Chromosome", "chromStart", "chromEnd", "cnv")
    cnv$Chromosome <- paste0('chr', cnv$Chromosome)
    cnv$GeneName <- "gene"
    cnv <- cnv[, c("Chromosome", "chromStart", "chromEnd", "GeneName", "cnv")]
}
                    
# add CNV heatmap track
if (nrow(cnv) != 0) {
  RCircos.Heatmap.Plot(cnv, data.col = 5, track.num = 1, side = "in")
}
                   
dev.off()


