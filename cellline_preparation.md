```bash
singularity shell --bind /data/rmulqueen/projects/scalebio_dcis ~/singularity/amethyst.sif
```

```R
source("/data/rmulqueen/projects/scalebio_dcis/tools/scalemet_dcis/src/amethyst_custom_functions.R") #to load in
set.seed(111)
options(future.globals.maxSize= 80000*1024^2) #80gb limit for parallelizing

task_cpus=300
project_data_directory="/data/rmulqueen/projects/scalebio_dcis/data/250815_milestone_v1"
merged_dat_folder="merged_data"
wd=paste(sep="/",project_data_directory,merged_dat_folder)

system(paste("mkdir -p",wd))
setwd(wd)

#make merged object
plate_obj=paste(sep="/",project_data_directory,"scalemethyl_pipeline_out/amethyst_plate_obj")
amethyst_files=list.files(path=plate_obj,pattern="*.amethyst.rds",recursive=TRUE,full.names=TRUE)

########################################
## Read in all amethyst files
########################################

window_name="cg_100k_score"
dat_list<-mclapply(amethyst_files, function(x) {
    obj<-readRDS(x)
    return(obj)},mc.cores=20)

dat <- combineObject(objList = dat_list, genomeMatrices=window_name)
```

########################################
## Run nakshatri atac peaks on all samples to cluster
########################################

```R
nakshatri_peaks<-readRDS("/data/rmulqueen/projects/scalebio_dcis/ref/celltype_peaks.500bp.rds")
nakshatri_bed<-data.frame(
        seqnames=seqnames(nakshatri_peaks),
        starts=start(nakshatri_peaks),
        ends=end(nakshatri_peaks))

nakshatri_bed<-nakshatri_bed[nakshatri_bed$seqnames %in% unique(dat@ref$seqid),]
dat@genomeMatrices[["nakshatri_peaks"]] <- makeWindows(dat, 
                                                     bed = nakshatri_bed,
                                                     type = "CG", 
                                                     metric = "score", 
                                                     threads = 300, 
                                                     index = "chr_cg", 
                                                     nmin = 2) 

saveRDS(dat,file="00_merged.amethyst.rds")

dat<-readRDS(dat,file="00_merged.amethyst.rds")
#subset to cell lines
dat_celllines<-subsetObject(dat,cells=row.names(dat@metadata)[dat@metadata$sample %in% c("MCF10A","MCF7","MDA-MB-231")])
saveRDS(dat_celllines,file="01_celllines.amethyst.rds")
