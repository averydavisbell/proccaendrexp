#! /usr/bin/env/ Rscript
# Get CeNDR RNA-seq (data from Zhang et al 2022) per-strain per-gene quantification estimates - after workflow run of strain-specific salmon
# By Avery Davis Bell, begun 2023.12.11
require(argparser, quietly = T)
require(data.table, quietly = T)
require(tximport, quietly = T)
require(DESeq2, quietly = T)
require(ggplot2, quietly = T)

#### Functions ####
tximportdata<-function(sampinfo, exampquantsf, tx2genef, tmpdir = paste0("tmp",sample.int(10e03, 1))){
  # ** from 
  # Gets in salmon RNA quantification data via tximport
  # In: sampinfo, data.table of sample information for samples to read in. Must include column SampleID
  #     exampquantsf, example filepath to salmon quant.sf (or quant.sf.gz) RNA quantifiaction file for one sample. Transcripts in name-sorted order. Where each Sample ID goes, needs to have _sampid_ (e.g. path/to/file/_sampid__genecounts.txt.gz).
  #               For any other differences in filepath, include * for interpolation.
  #     tx2genef, Path to file mapping transcripts to genes. Two columns (transcript ID, gene ID), no header.
  #     tmpdir, directory that will be created for temporary salmon files. Created and deleted by function.
  # Out: list of - 
  #   $missingts, record of transcripts that weren't in 1 or more input files and had to be added (non-reference strains). Columns:
  #      transcript_id, fromgene (gene this comes from), nSampMissingIn (# samples that didn't have this - multiple of strains that didn't have it),
  #      ntranscriptsthisgene, # transcripts this gene has, in service of:
  #      genemissing, T or F - T if this gene only has this one transcript
  #      MissingIn, sample IDs that had transcript missing
  #   $txi, tximport output for this data. To be passed directly to DESeqDataSetFromTximport
  #   $tx2gene, transcript ID to gene ID mapping data.table. Two columns: transcript_id, gene_id
  
  # Subfunctions [from previous script!]
  getalltrs<-function(sfile, tx2gene, ofile){
    # Reads salmon file, adds trancsripts that weren't included, save which ones these were
    # Writes OUT to ofile
    
    s<-fread(sfile)
    setkey(s, Name)
    s<-s[tx2gene$tr,]
    
    missingts<-tx2gene[tr%in%s[is.na(Length), Name],]
    if(nrow(missingts)>0){
      s[which(Name%in%missingts$tr), `:=`(Length = 10, EffectiveLength = 10, TPM = 0, NumReads = 0)]
      
      write.table(s, ofile, sep = "\t", row.names = F, quote = F)
      return(list(newfile = ofile, missingts = missingts, tx2gene = tx2gene))
    }else{
      return(NULL)
    }
  }
  
  # Get all files
  myfiles<-sapply(sampinfo$SampleID, function(x){
    Sys.glob(gsub("_sampid_", x, exampquantsf))
  })
  
  # Make sure all files have all transcripts (they don't necessarily if strain-specific/different strains!)
  ## Get in transcripts
  tx2gene<-fread(tx2genef, header = F)
  setnames(tx2gene, c("tr","g"))
  setkey(tx2gene, "tr")
  ## New files with all transcripts are temporarily generated: Write TEMPORARY salmon files that have ALL transcripts, zeroing out those that don't exist in alt strains
  if(!dir.exists(tmpdir)){dir.create(tmpdir, recursive = T)}
  newfs<-lapply(names(myfiles), function(x) getalltrs(myfiles[[x]], tx2gene = tx2gene, ofile = paste0(tmpdir, "/", x, "_quant.sf")))
  ## Update to read these files where needed
  files.use<-sapply(1:length(myfiles), function(x) ifelse(is.null(newfs[[x]]), myfiles[[x]], newfs[[x]]$newfile))
  ## Format missing transcripts to save/return
  missingts<-rbindlist(lapply(1:length(myfiles), function(x){
    if(!is.null(newfs[[x]])){
      out<-newfs[[x]]$missingts
      out[, Samp:=names(myfiles)[x]]
      setnames(out, c("tr", "g"), c("transcript_id", "gene_id"))
    }else{
      return(NULL)
    }
  }))
  setkey(missingts, "transcript_id")
  missingts<-missingts[,.(gene_id = gene_id[1], nSampsMissingIn = .N, MissingIn = paste(Samp, collapse = ",")), by = "transcript_id"]
  ### Record if this is only transcript from gene
  setkey(tx2gene, g)
  tperg<-tx2gene[, .N, by = g]
  setkey(missingts, gene_id)
  missingts<-tperg[missingts]
  setnames(missingts, c("g", "N"), c("fromgene", "ntranscriptsthisgene"))
  missingts[,genemissing:=(ntranscriptsthisgene==1)]
  setcolorder(missingts, c("transcript_id", "fromgene", "nSampsMissingIn", "ntranscriptsthisgene", "genemissing"))
  
  # Read with tximport 
  setkey(tx2gene, "tr")
  names(files.use)<-names(myfiles)
  txi.salmon <- tximport(files.use, type = "salmon", tx2gene = tx2gene)
  
  # Delete TEMPORARY salmon files
  unlink(tmpdir, recursive = T)
  
  # Return
  setnames(tx2gene, c("transcript_id", "gene_id"))
  return(list(missingts = missingts, txi = txi.salmon, t2g = tx2gene))
}

plotPCA_givePCs<-function (object, intgroup = "condition", ntop = 500, returnData = FALSE, xpc = 1, ypc = 2, mytitle = "") {
  # This is DESeq2's plotPCA function (accessed via ESeq2:::plotPCA.DESeqTransform) with added functionality to plot other than first 1 or 2 PCs
  # In: see ?plotPCA except for:
  #   xpc: integer. Which PC to plot on X axis.
  #   ypc: integer. Which PC to plot on Y axis
  #   mytitle: character. Title for plot.
  
  xpcname<-paste0("PC", xpc)
  ypcname<-paste0("PC", ypc)
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = ":"))
  }
  else {
    colData(object)[[intgroup]]
  }
  d <- data.frame(xpcname = pca$x[, xpc], ypcname = pca$x[, ypc], group = group, 
                  intgroup.df, name = colnames(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[xpc:ypc]
    return(d)
  }
  ggplot(data = d, aes_string(x = "xpcname", y = "ypcname", color = "group")) + 
    geom_point(size = 3) + 
    xlab(paste0(xpcname, ": ", round(percentVar[xpc] * 100), "% variance")) + 
    ylab(paste0(ypcname, ": ", round(percentVar[ypc] * 100), "% variance")) + 
    coord_fixed() + 
    ggtitle(mytitle) + 
    theme(plot.title = element_text(hjust = 0.5))
}

#### Arguments & input ####
pp<-arg_parser("Get CeNDR RNA-seq (data from Zhang et al 2022) per-strain per-gene quantification estimates - after workflow run of strain-specific salmon", 
               name = "cendrrnaseq2022paper_salm2vst.R",
               hide.opts = T)
pp<-add_argument(pp, "--sampinfo",
                  help = "Path to sample information file. Must include column SampleID and RefDescrip (isotype for each sample - this is the name each strain will carry downstream!).",
                  type = "character")
pp<-add_argument(pp, "--baseoutname",
                help = "Base name for all output files",
                type = "character",
                default = "out")
pp<-add_argument(pp, "--outdir",
                help = "Outer output directory. Sub-directories will be created internally.",
                type = "character")
pp<-add_argument(pp, "--exampquantsf",
                help = "example filepath to salmon quant.sf (or quant.sf.gz) RNA quantifiaction file for one sample. Transcripts in name-sorted order. Where each Sample ID goes, needs to have _sampid_ (e.g. path/to/file/_sampid__genecounts.txt.gz).
                For any other differences in filepath, include * for interpolation.",
                type = "character")
pp<-add_argument(pp, "--tx2genef",
                help = "Path to file mapping transcripts to genes. Two columns (transcript ID, gene ID), no header.",
                type = "character")
pp<-parse_args(pp)

#### Process data ####
# Sample info
sampinfo<-fread(pp$sampinfo, select = c("SampleID", "RefDescrip"))
setnames(sampinfo, "RefDescrip", "Isotype")
rownames(sampinfo)<-sampinfo$SampleID # For DESeq2

# Import salmon data
dat<-tximportdata(sampinfo, exampquantsf = pp$exampquantsf, tx2genef = pp$tx2genef,
                  tmpdir = paste0("tmp",sample.int(10e03, 1)))

# DESeq2 format
dds<-DESeqDataSetFromTximport(dat$txi,
                              colData = sampinfo,
                              design = ~ Isotype)
## Remove genes with too few reads; those missing in some samples (due to strain-specific transcriptome)
##    [may get added back in as 0s later]
g.few<-which(rowSums(counts(dds)) < 10)
g.miss<-which(rownames(dds) %in% dat$missingts[genemissing==T, fromgene])
dds<-dds[-(unique(c(g.few, g.miss))), ]
cat(paste(nrow(dds), "genes remain after removal of", length(g.few), "genes with < 10 reads across samples, including",
          sum(g.miss%in%g.few),"genes with < 10 reads across samples AND missing in some samples' input;",
          sum(!g.miss%in%g.few), "genes missing in some samples' input but having enough reads.", "\n"))

# ---- do variance stabilizing transform
vsd<-vst(dds, blind = T)

# ---- do PCA for my own QC

plot1<-plotPCA_givePCs(vsd, intgroup = "Isotype", ntop = 500, returnData = F, xpc = 1, ypc = 2,
                       mytitle = "All CeNDR RNA-seq samples (PCA Plot, variance stabilizing transformed, top 500 most variable genes)")
plot2<-plotPCA_givePCs(vsd, intgroup = "Isotype", ntop = 500, returnData = F, xpc = 2, ypc = 3,
                       mytitle = "All CeNDR RNA-seq samples (PCA Plot, variance stabilizing transformed, top 500 most variable genes)")
plot3<-plotPCA_givePCs(vsd, intgroup = "Isotype", ntop = 500, returnData = F, xpc = 3, ypc = 4,
                       mytitle = "All CeNDR RNA-seq samples (PCA Plot, variance stabilizing transformed, top 500 most variable genes)")
pdf(file.path(pp$outdir, paste0(pp$baseoutname, "_PCAPlots.pdf")), 16, 12)
print(plot1 + theme(legend.position = "none"))
print(plot2 + theme(legend.position = "none"))
print(plot3 + theme(legend.position = "none"))
invisible(dev.off())

# ---- look more into PC outliers
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(500, 
                                                   length(rv)))]
pca <- prcomp(t(assay(vsd)[select, ]))
pcs<-data.table(pca$x, keep.rownames = T)
# pcs[PC1>30, .(rn, PC1, PC2)]
exclsamps<-pcs[PC1>35, rn] ## going to add. Checked any of these before doing.
# DROP
vsd<-vsd[, -which(colnames(vsd)%in%exclsamps)]
sampinfodropped<-sampinfo[!SampleID%in%exclsamps, ]

# --- Get per-isotype means
datformean<-data.table(assay(vsd), keep.rownames = T)
setnames(datformean, "rn", "gene_id")

isomeans<-do.call(cbind,
                  lapply(sampinfo[,unique(Isotype)], function(iso){
                    thisdat<-datformean[, sampinfodropped[Isotype==iso, SampleID], with = F]
                    out<-data.table(tmp = rowMeans(thisdat))
                    setnames(out, "tmp", iso)
                    return(out)
                    }))
# Add in gene ID
isomeans[, gene_id:=datformean[, gene_id]]
setkey(isomeans, gene_id)
setcolorder(isomeans)

# Save
write.table(isomeans, gzfile(file.path(pp$outdir, paste0(pp$baseout, "_meanvstexpperisotype.txt.gz"))),
            sep = "\t", quote = F, row.names = F)

# ** write output with isotype not strain
# Save sample-level per gene and strain level? Not sure why I'd want sample level though, anything I'd do with that I'd probably want to be in DESeq2 framework

# *** saving missing transcripts info? (it's less that they're 0 than missing....I added in as 0...)

# If gene not here, assign it 0 [or some value from input!! but prob 0?]
#       can get full gene list from from tx2genef

#### Print out session info for reproducibility ####
cat("Processing complete. R session info:\n")
sessionInfo()