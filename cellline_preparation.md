# Environment for processing on MDA servers.

```bash
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
```

# Loading in packages and custom fuctions
(see scalemet_dcis repo for more details.)

```R
source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
workshop_outputdir="/data/rmulqueen/projects/scalebio_dcis/workshop"
system(paste("mkdir -p",workshop_outputdir))
setwd(workshop_outputdir)
dat_celllines<-readRDS(file=paste(sep="/",project_data_directory,merged_dat_folder,"01_celllines.amethyst.rds"))#previously processed from scalebio_dcis

#subset to number that can be processed reasonably in a workshop.
#using cells with shared dedup.bam

samp<-lapply(unique(dat_celllines@metadata$sample), function(x){
    sort(table(dat_celllines@metadata[dat_celllines@metadata$sample==x,]$bam_path))
    dat_celllines@metadata[sample(which(dat_celllines@metadata$sample==x),size=50),]
})

#~50 MCF10A
#/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/scalemethyl_pipeline_out/prelim1-2_scalebio_plate1-1/alignments/dedup/MCF10A.1A06/MCF10A.1A06.dedup.bam
#/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/scalemethyl_pipeline_out/prelim1-2_scalebio_plate1-1/alignments/dedup/MCF10A.1B03/MCF10A.1B03.dedup.bam

#~50 MCF7
#/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/scalemethyl_pipeline_out/prelim1-2_scalebio_plate1-1/alignments/dedup/MCF7.1D03/MCF7.1D03.dedup.bam
#/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/scalemethyl_pipeline_out/prelim1-2_scalebio_plate1-1/alignments/dedup/MCF7.1C01/MCF7.1C01.dedup.bam

#~50 MDA-MB-231
#/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/scalemethyl_pipeline_out/prelim1-2_scalebio_plate1-1/alignments/dedup/MDA-MB-231.1F03/MDA-MB-231.1F03.dedup.bam
#/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/scalemethyl_pipeline_out/prelim1-2_scalebio_plate1-1/alignments/dedup/MDA-MB-231.1H09/MDA-MB-231.1H09.dedup.bam

bam_list=c("/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/scalemethyl_pipeline_out/prelim1-2_scalebio_plate1-1/alignments/dedup/MCF10A.1A06/MCF10A.1A06.dedup.bam","/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/scalemethyl_pipeline_out/prelim1-2_scalebio_plate1-1/alignments/dedup/MCF10A.1B03/MCF10A.1B03.dedup.bam","/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/scalemethyl_pipeline_out/prelim1-2_scalebio_plate1-1/alignments/dedup/MCF7.1D03/MCF7.1D03.dedup.bam","/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/scalemethyl_pipeline_out/prelim1-2_scalebio_plate1-1/alignments/dedup/MCF7.1C01/MCF7.1C01.dedup.bam","/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/scalemethyl_pipeline_out/prelim1-2_scalebio_plate1-1/alignments/dedup/MDA-MB-231.1F03/MDA-MB-231.1F03.dedup.bam","/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1/scalemethyl_pipeline_out/prelim1-2_scalebio_plate1-1/alignments/dedup/MDA-MB-231.1H09/MDA-MB-231.1H09.dedup.bam")

filtered_cells<-row.names(dat_celllines@metadata[dat_celllines@metadata$bam_path %in% bam_list,])
dat_celllines<-subsetObject(dat_celllines,cells=filtered_cells)

#write out h5 files and bam files unique to dat_celllines after filtering. This is the raw data for the workshop.
#write out the metadata from dat_celllines as well.
system(paste0("mkdir -p ",workshop_outputdir))
system(paste0("mkdir -p ",workshop_outputdir,"/h5_files"))
#system(paste0("mkdir -p ",workshop_outputdir,"/bam_files"))

#h5 files copy
h5_paths=dat_celllines@h5paths
h5_paths$path<-paste0("/content/h5_files/",basename(h5_paths$path)) #changing h5 paths name to workshop colab notebook
h5_paths$barcode<-row.names(h5_paths) #changing h5 paths name to workshop colab notebook
h5_paths<-h5_paths[,c("barcode","path")]
write.table(h5_paths,file=paste0(workshop_outputdir,"/h5_files/h5_paths.tsv"),row.names=T,col.names=T,sep="\t")
lapply(unique(dat_celllines@h5paths$paths),function(x){system(paste("cp",x,paste0(workshop_outputdir,"/h5_files")))})

#bam files copy
#lapply(unique(dat_celllines@metadata$bam_path),function(x){system(paste("cp",x,paste0(workshop_outputdir,"/bam_files")))})
#might preprocess the bams further into just read counts 
#yep defintely do, just to make the data much smaller.

met<-dat_celllines@metadata
met<-met[,c("cell_id","tgmt_well","i5_well","i7_well","plate","batch","plate_info",
"unique_reads","total_reads","mito_reads",
"mcg_pct","cov","ch_cov","tss_enrich","cg_cov","mch_pct","sample" )]
write.table(met,file=paste0(workshop_outputdir,"/metadata.tsv"),row.names=T,col.names=T,sep="\t")
write.table(met,file="metadata.tsv",row.names=T,col.names=T,sep="\t")
```



# Prerun 5kb windows for cell lines for visualization
(running on not-downsampled by cell count)

```R
dat_celllines<-readRDS(file=paste(sep="/",project_data_directory,merged_dat_folder,"01_celllines.amethyst.rds"))#previously processed from scalebio_dcis

celltype500bpwindows <- calcSmoothedWindows(dat_celllines, 
                                        type = "CG", 
                                        threads = 10,
                                        step = 500, 
                                        smooth = 3,
                                        genome = "hg38",
                                        index = "chr_cg",
                                        groupBy = "sample",
                                        returnSumMatrix = TRUE, # save sum matrix for DMR analysis
                                        returnPctMatrix = TRUE)
saveRDS(celltype500bpwindows,file="celltype_500bpwindow_tracks.rds")
```



# Pregenerate Copykit RDS to decrease bam file transfer

```bash
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
```

# Generate CopyKit for each sample

```R
source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
library(Rsamtools)
library(copykit)
library(circlize)
detach("package:GeneNMF",unload=TRUE)
library(dendextend)
library(RColorBrewer)
library(ComplexHeatmap)

#set environment and read in data
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing
task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
obj<-readRDS(file=paste(sep="/",project_data_directory,merged_dat_folder,"01_celllines.amethyst.rds"))

#set output position
workshop_outputdir="/data/rmulqueen/projects/scalebio_dcis/workshop"
setwd(workshop_outputdir)

read_scalebio_bam<-function(obj_met,x,sample_name){
    #scalebio pipeline outputs bam files as Tn5 wells. so multiple cell IDs are in a bam. this function splits out the bam to the query cellid
    bam=obj_met[obj_met$sample %in% c(sample_name),]$bam_path[x]
    cellid=strsplit(row.names(obj_met)[x],"[+]batch|[+]prelim")[[1]][1]

    what <- c("qname","rname", "pos")
    param <- ScanBamParam(what=what,
                            flag=scanBamFlag(isPaired=TRUE,
                                            isProperPair=TRUE,
                                            isSecondaryAlignment=FALSE,
                                            isDuplicate=FALSE,
                                            isSupplementaryAlignment=FALSE))

    input_bam<-Rsamtools::scanBam(bam,param=param)
    input_bam<-do.call("DataFrame", input_bam)
    input_bam$cellid<-gsub("^.*:", "", input_bam$qname)
    input_bam<-input_bam[input_bam$cellid==cellid,]
    input_bam$end<-input_bam$pos+1
    input_bam<-makeGRangesFromDataFrame(input_bam,seqnames.field="rname",start.field="pos",end.field="end")
    return(input_bam)
}

genome = "hg38"
resolution="220kb"
remove_Y = TRUE
min_bincount = 10
cores=100
# bindings for NSE and data
Chr <- chr <- strand <- GeneID <- NULL
reads_assigned_bins <- reads_duplicates <- reads_total <- NULL

# Reading hg38 VarBin ranges
hg38_grangeslist <- hg38_grangeslist

hg38_rg <- switch(resolution,
    "55kb" = hg38_grangeslist[["hg38_50kb"]],
    "110kb" = hg38_grangeslist[["hg38_100kb"]],
    "195kb" = hg38_grangeslist[["hg38_175kb"]],
    "220kb" = hg38_grangeslist[["hg38_200kb"]],
    "280kb" = hg38_grangeslist[["hg38_250kb"]],
    "500kb" = hg38_grangeslist[["hg38_500kb"]],
    "1Mb" = hg38_grangeslist[["hg38_1Mb"]],
    "2.8Mb" = hg38_grangeslist[["hg38_2Mb"]]
)

hg38_rg <- as.data.frame(hg38_rg)

rg <- hg38_rg %>%
    dplyr::rename(chr = "seqnames") %>%
    dplyr::mutate(GeneID = 1:nrow(hg38_rg))

if (remove_Y == TRUE) {
    rg <- dplyr::filter(rg,chr != "chrY")
}

message("Counting reads for genome ",genome," and resolution: ",resolution)

#get list of bams and cellids
obj_met<-read.csv("metadata.tsv",header=T,sep="\t")
sample_name=unique(obj_met$sample)
obj_met$cellline<-obj_met$sample

#return chr start position for reads filtered in bam to cell id
varbin_counts_list_all_fields<-mclapply(
                                    1:nrow(obj_met), 
                                    function(i) 
                                    read_scalebio_bam(obj_met=obj_met,x=i,sample_name=sample_name), 
                                    mc.cores=cores)

message("Read in all bam files.")

names(varbin_counts_list_all_fields)<- row.names(obj_met)
varbin_counts_list_all_fields<-as(varbin_counts_list_all_fields, "GRangesList")
ref<-as(rg,"GRanges")

varbin_counts_list <-mclapply(varbin_counts_list_all_fields,
                                function(x) 
                                GenomicRanges::countOverlaps(
                                query=ref,
                                subject=x,
                                type="any",
                                ignore.strand=TRUE),
                                mc.cores=cores)
message("Counted reads across all bins.")

varbin_counts_list <- lapply(varbin_counts_list,as.vector)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# filtering for minimal mean bin count
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# obtaining the index of the ones that FAIL to meet the min_bincount arg
min_bc <- which(vapply(varbin_counts_list, mean, numeric(1)) < min_bincount)
# subsetting counts list and main counts list

if (length(min_bc) > 0) {
    varbin_counts_list <- varbin_counts_list[-min_bc]
    varbin_counts_list_all_fields <- varbin_counts_list_all_fields[-min_bc]
    message(
        length(min_bc), " bam files had less than ", min_bincount,
        " mean bincounts and were removed."
    )
}

# LOWESS GC normalization

message("Performing GC correction.")

varbin_counts_list_gccor <-
    mclapply(varbin_counts_list, function(x) {
        gc_cor <- lowess(rg$gc_content, log(x + 1e-3), f = 0.05)
        gc_cor_z <- approx(gc_cor$x, gc_cor$y, rg$gc_content)
        exp(log(x) - gc_cor_z$y) * median(x)
    },mc.cores=cores
    )

varbin_counts_df <- round(dplyr::bind_cols(varbin_counts_list_gccor), 2)

# filtering low read counts where the sum of bins does not reach more than 0
good_cells <- names(varbin_counts_df[which(colSums(varbin_counts_df) != 0)])

varbin_counts_df <- varbin_counts_df[good_cells]

rg <- rg %>%
    dplyr::select(-strand, -GeneID)

rg_gr <- GenomicRanges::makeGRangesFromDataFrame(rg,
    ignore.strand = TRUE,
    keep.extra.columns = TRUE)


cna_obj <- CopyKit(
    assays = list(bincounts = varbin_counts_df),
    rowRanges = rg_gr)

# Adding genome and resolution information to metadata
S4Vectors::metadata(cna_obj)$genome <- genome
S4Vectors::metadata(cna_obj)$resolution <- resolution

# saving info and removing columns from list elements
bam_metrics <- obj_met[c("unique_reads","tss_enrich","mcg_pct","cg_cov","batch","plate_info","tgmt_well","i7_well","i5_well","sample")]

# making sure metrics match varbin_counts_df
bam_metrics <- bam_metrics[good_cells,]

bam_metrics$sample <- rownames(bam_metrics)
bam_metrics$sample_name=sample_name
bam_metrics$reads_assigned_bins <- colSums(varbin_counts_df)

# adding to metadata
SummarizedExperiment::colData(cna_obj) <-
    S4Vectors::DataFrame(bam_metrics)
colnames(cna_obj) <- names(varbin_counts_df)

saveRDS(cna_obj,file="workshop.copykit.obj.rds")
```
## Now archive it all together and upload to drive

Need to update the h5 files to new format
```bash
#pip install facet
#pip install amethyst-facet
for i in ./h5_files/*h5; do /home/rmulqueen/.local/bin/facet convert -c CG ${i::-3}.workshop.h5 $i; done
tar -czf workshop_data.tar.gz *
```


