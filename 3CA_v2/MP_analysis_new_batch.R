
####  Run MP analysis for new datasets not included in the Nature papaer

## (1) Prepare new samples for NMF
## (2) Run NMF - using cluster 
## (3) Generate MPs using the new NMFs together with the previous ones 


source("~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/Return_All_cells_Separated.R")



New_datasets <- c("Andrade_2019_Melanoma_scRNAseq",
                  "Azizi_2018_Breast_scRNAseq",
                  "Bassez_2021_Breast_scRNAseq",
                  "Bhan_2018_HCC_CTC_scRNAseq",
                  "Biermann_2022_Melanoma_scRNAseq",
                  "Bischoff_2021_LUAD_scRNAseq",
                  "Chan_2021_SCLC_scRNAseq",
                  "Chen_2021_Colon_dis_scRNAseq",
                  "Chen_2021_Colon_val_scRNAseq",
                  "Choudhury_2022_Meningioma_scRNAseq",
                  "Cillo_2020_HNSCC_scRNAseq",
                  "Durante_2020_Uveal_Melanoma_reduced_scRNAseq",
                  "Gaydosik_2019_CTCL_scRNAseq",
                  "Geistlinger_2020_Ovarian_scRNAseq",
                  "Gonzalez_2022_BrMs_scRNAseq",
                  "Griffiths_2021_Breast_snRNAseq",
                  "Gubin_2018_Models_Sarcoma_scRNAseq",
                  "Gulati_2020_Breast_scRNAseq",
                  "Guo_2018_Lung_scRNAseq",
                  "He_2021_Prostate_scRNAseq",
                  "Hollern_2019_TNBC_mouse_scRNAseq",
                  "Hwang_2022_PDAC_scRNAseq",
                  "Kim_2018_TNBC_snRNAseq",
                  "Krishna_2021_RCC_scRNAseq",
                  "Kuerten_2021_HNSCC_scRNAseq",
                  "Kumar_2022_Gastric_scRNAseq",
                  "Li_2019_Melanoma_scRNAseq",
                  "Mahuron_2020_Melanoma_scRNAseq",
                  "Nassiri_2021_Meningioma_snRNAseq",
                  "Neal_2018_Kidney_scRNAseq",
                  "Olbrecht_2021_Ovarian_scRNAseq",
                  "Paulson_2020_MCC_scRNAseq",
                  "Pelka_2021_CRC_scRNAseq",
                  "Raghavan_2021_PDAC_scRNAseq_Model",
                  "Raghavan_2021_PDAC_scRNAseq_Tumor",
                  "Rao_2020_NE_scRNAseq",
                  "Regner_2021_Gyne_scRNAseq",
                  "Rendeiro_2020_CLL_scRNAseq",
                  "Riether_2020_AML_scRNAseq",
                  "Roider_2020_Lymphoma_scRNAseq",
                  "Sade_Feldman_Melanoma_2018_scRNAseq",
                  "Sathe_2020_Gastric_scRNAseq",
                  "Savas_2018_TNBC_scRNAseq",
                  "Sharma_2020_HCC_scRNAseq",
                  "Shih_2018_Ovarian_scRNAseq",
                  "Song_2019_Lung_scRNAseq",
                  "Song_2022_Prostate_scRNAseq",
                  "Steen_2021_DLBCL_scRNAseq",
                  "Sun_2021_HCC_CTC_scRNAseq",
                  "Tang_Huau_2018_Ovarian_scRNAseq",
                  "Wu_2021_Breast_scRNAseq",
                  "Yost_2019_Skin_scRNAseq",
                  "Yuen_2020_Kidney_scRNAseq",
                  "Zhang_2018_CRC_scRNAseq",
                  "Zhang_2019_FL_scRNAseq",
                  "Zhang_2019_Gastric_scRNAseq",
                  "Zhang_2019_HGSOC_scRNAseq",
                  "Zhang_2019_Liver_scRNAseq",
                  "Zhang_2022_Ovarian_scRNAseq",
                  "Zheng_2017_Liver_scRNAseq",
                  "Zilionis_2019_Lung_scRNAseq"
)



Previous_studies <- list.files("~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/Output")
Previous_studies[which(Previous_studies == "Antonella_2020_lungcancer_scRNAseq")] <- "Dost_2020_lungcancer_scRNAseq"
Previous_studies <- sort(Previous_studies)


df <- data.frame(Previous_studies = rep("" , max(length(Previous_studies) , length(New_datasets) )),
                 New_studies      = rep("" , max(length(Previous_studies) , length(New_datasets) )))
df$Previous_studies[1:length(Previous_studies)] <- Previous_studies
df$New_studies[1:length(New_datasets)]          <- New_datasets

write.csv(df , file = "~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/Studies_for_new_NMF.csv")



######  Prepare data for NMF


Process_mat <- function (Mat) {   ### filter genes, log, center, assign 0 to negative
  
  gene_inds_sorted <- sort(log2(rowMeans(Mat) + 1) , decreasing = TRUE , index.return = TRUE)$ix[1:7000] 
  Mat   <- Mat[gene_inds_sorted,]
  # Mat   <- log2((Mat/10) + 1)
  Mat   <- Mat - rowMeans(Mat)
  Mat[Mat < 0]   <- 0 
  if (length(which(rowSums(Mat)==0))>0)   Mat <- Mat[-which(rowSums(Mat)==0) , ]   ## nmf does not work if all the values in one row are zero.
  if (length(which(colSums(Mat)==0))>0)   Mat <- Mat[ , -which(colSums(Mat)==0)]   ## nmf does not work if all the values in one column are zero.
  
  return(Mat)
}


for (I in New_datasets){
  
  Curr_data <-   Return_All_cells_Separated(I , "Malignant", filter_genes = FALSE , log_normalize = TRUE)  ## returns list of cancer cells per sample, after cell filtering, no gene filtering, after log-normalization, not centered

  if (!is.na(Curr_data)){
    Curr_data <- lapply(Curr_data, function(x) Process_mat(x))
    path      <- paste0("~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/Data_new_batch/" , I , "_data.RDS")
    saveRDS(Curr_data , file = path)
  }
  
  print(I)
  
}




A <- readRDS("~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/Data_new_batch/Biermann_2022_Melanoma_scRNAseq_data.RDS")

ee <- as.vector(unlist(lapply(A, function(x) dim(x)[2])))
which(ee > 4000)




####### For checking if NMFs are complete


### Check if exists 

i = "Biermann_2022_Melanoma_scRNAseq"
file1 <- paste0("~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/Data_new_batch/", i,"_data.RDS")
R1 <- readRDS(file1)
R0 <- R1
R1 <- paste0(names(R1), "_rank4_9_nrun10.RDS")

file2 <- paste0("~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/Output_new_batch/", i)
R2 <- list.files(file2)

R3 <- which(is.na(match(R1,R2)))
R3


### *** Check that all output files have 4-9 ranks 

setwd(paste0("~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/Output_new_batch/", i))
p <- 0
inds_incomplete <- c()


for (j in R1) {
  
  p <- p+1 
  
  W <- j
  
  if (file.exists(W)){
    
    curr_data <- readRDS( W )
    # if (sum(as.integer(names(curr_data$fit))) < 39)                     inds_incomplete <- c(inds_incomplete , p)
    if (length(which(is.na(curr_data$fit[as.character(4:9)]))) > 0)     inds_incomplete <- c(inds_incomplete , p)
  } 
  
}
inds_incomplete




##### Check for all datasets that NMF has completed. count number of samples and number of cells

a <- list.files("~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/Output_new_batch/")

study_not_complete <- c()
sample_num <- 0
cell_num <- 0
for (i in a){
  
  file1 <- paste0("~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/Data_new_batch/", i,"_data.RDS")
  R1 <- readRDS(file1)
  R0 <- R1
  R1 <- paste0(names(R1), "_rank4_9_nrun10.RDS")
  
  file2 <- paste0("~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/Output_new_batch/", i)
  R2 <- list.files(file2)
  
  R3 <- which(is.na(match(R1,R2)))
  
  if (length(R3)>0) { study_not_complete <- c(study_not_complete , i) }
  
  
  setwd(paste0("~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/Output_new_batch/", i))
  p <- 0
  inds_incomplete <- c()
  
  
  for (j in R2) {
    
    p <- p+1 
    
    W <- j
    
    if (file.exists(W)){
      
      curr_data <- readRDS( W )
      # if (sum(as.integer(names(curr_data$fit))) < 39)                     inds_incomplete <- c(inds_incomplete , p)
      if (length(which(is.na(curr_data$fit[as.character(4:9)]))) > 0)     inds_incomplete <- c(inds_incomplete , p)
    } 
    
  }
  
  if (length(inds_incomplete)>0) { study_not_complete <- c(study_not_complete , i) }
  
  sample_num <- sample_num + length(R0)
  cell_num   <- cell_num   + sum(as.vector(unlist(  lapply(R0, function(x) dim(x)[2] )  )))
  
}






### Generate MPs combining old and new batch 

### NMF comparison for Pan-Cancer dataset 

# ----------------------------------------------------------------------------------------------------
# Pan-Cancer :  Definying heterogeneity patterns that are shared between datasets 
# ----------------------------------------------------------------------------------------------------


# load necessary R packages and functions
setwd("~/../avishaig/R_WD/Pan_cancer/My_scripts")
library(reshape2)
library(ggplot2)
library(scales)
library(egg)
source("robust_nmf_programs.R")
source("custom_magma.R")
library(NMF)
library(RColorBrewer)

setwd("~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF") # location of function
source("sample_names_abbrev.R")




# ----------------------------------------------------------------------------------------------------

Genes_nmf_w_basis <- list() # nmf gene scores
Cells_nmf_h_coef <- list() # nmf cell scores


### Previous

Data_list <- list.files(path = "~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/Output")   ### directories with NMF data 

for (i in Data_list){
  
  Sample_list <- list.files(path = paste0("~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/Output/",i))
  q =  as.data.frame(apply(as.data.frame(Sample_list) , 1 , function(x)  strsplit(x, "[.]")[[1]][length(strsplit(x, "[.]")[[1]])])); colnames(q) <- "V1"
  Sample_list <- Sample_list[which(q$V1=="RDS")]
  
  
  for (k in Sample_list){
    
    curr_data <- readRDS( file = paste0("~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/Output/",i,"/",k) )
    
    m1 <- c() # temp genes
    m2 <- c() # temp cells 
    for (j in 4:9){  ## some ranks had NULL which give an error(only for coefficients command), so first check in not NULL 
      if (! is.null(  basis(curr_data$fit[[as.character(j)]]) )) { 
        k1 <- basis(curr_data$fit[[as.character(j)]])
        k2 <- t(coefficients(curr_data$fit[[as.character(j)]]))
        colnames(k1)  <- paste0(k, ".", j, ".", 1:j)
        colnames(k2)  <- paste0(k, ".", j, ".", 1:j)
        
        m1           <- cbind(m1 , k1 )
        m2           <- cbind(m2 , k2 )  }
      
    }
    
    Genes_nmf_w_basis[[k]] <- m1
    Cells_nmf_h_coef[[k]]  <- m2
    
  }
  
}


### New


Data_list <- list.files(path = "~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/Output_new_batch")   ### directories with NMF data 

for (i in Data_list){
  
  Sample_list <- list.files(path = paste0("~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/Output_new_batch/",i))
  q =  as.data.frame(apply(as.data.frame(Sample_list) , 1 , function(x)  strsplit(x, "[.]")[[1]][length(strsplit(x, "[.]")[[1]])])); colnames(q) <- "V1"
  Sample_list <- Sample_list[which(q$V1=="RDS")]
  
  
  for (k in Sample_list){
    
    curr_data <- readRDS( file = paste0("~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/Output_new_batch/",i,"/",k) )
    
    m1 <- c() # temp genes
    m2 <- c() # temp cells 
    for (j in 4:9){  ## some ranks had NULL which give an error(only for coefficients command), so first check in not NULL 
      if (! is.null(  basis(curr_data$fit[[as.character(j)]]) )) { 
        k1 <- basis(curr_data$fit[[as.character(j)]])
        k2 <- t(coefficients(curr_data$fit[[as.character(j)]]))
        colnames(k1)  <- paste0(k, ".", j, ".", 1:j)
        colnames(k2)  <- paste0(k, ".", j, ".", 1:j)
        
        m1           <- cbind(m1 , k1 )
        m2           <- cbind(m2 , k2 )  }
      
    }
    
    Genes_nmf_w_basis[[k]] <- m1
    Cells_nmf_h_coef[[k]]  <- m2
    
  }
  
}


saveRDS(Genes_nmf_w_basis    , file = "Genes_nmf_w_basis_New_Batch.RDS")
saveRDS(Cells_nmf_h_coef    , file = "Cells_nmf_h_coef_New_Batch.RDS")


### ************************************************************************** 


# get gene programs (top 50 genes by NMF score)
nmf_programs_sig_ccle <- lapply(Genes_nmf_w_basis, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
nmf_programs_sig_ccle <- lapply(nmf_programs_sig_ccle,toupper) ## convert all genes to uppercase 

# for each sample, select robust NMF programs (i.e. obseved using different ranks in the same sample), remove redundancy due to multiple ranks, and apply a filter based on the similarity to programs from other samples. 
nmf_filter_ccle <- robust_nmf_programs(nmf_programs_sig_ccle, intra_min = 35, intra_max = 10, inter_filter=T, inter_min = 10)
saveRDS(nmf_filter_ccle    , file = "nmf_filter_ccle_New_Batch.RDS")

nmf_programs_sig_ccle <- lapply(nmf_programs_sig_ccle, function(x) x[, is.element(colnames(x), nmf_filter_ccle),drop=F])
nmf_programs_sig_ccle <- do.call(cbind, nmf_programs_sig_ccle)

# calculate similarity between programs
nmf_intersect_ccle <- apply(nmf_programs_sig_ccle , 2, function(x) apply(nmf_programs_sig_ccle , 2, function(y) length(intersect(x,y)))) 

# hierarchical clustering of the similarity matrix 
nmf_intersect_hc_ccle <- hclust(as.dist(50-nmf_intersect_ccle), method="average") 
nmf_intersect_hc_ccle <- reorder(as.dendrogram(nmf_intersect_hc_ccle), colMeans(nmf_intersect_ccle))
nmf_intersect_ccle <- nmf_intersect_ccle[order.dendrogram(nmf_intersect_hc_ccle), order.dendrogram(nmf_intersect_hc_ccle)]

abbreviated_sample_names <- sample_names_abbrev(colnames(nmf_intersect_ccle))

saveRDS(nmf_intersect_ccle    , file = "nmf_intersect_ccle_New_Batch.RDS")
saveRDS(nmf_programs_sig_ccle , file = "nmf_programs_sig_ccle_New_Batch.RDS")




###  This function receives as input the intersection matrix nmf_intersect_ccle and nmf_programs_sig_ccle (all the NMFs)
###  Parameters are : 
# (1) Minimum overlap between NMF programs for defining initial groups. Default is 10
# (2) Minimum overlap for adding an NMF program to a cluster. Default is 10
# (3) Minimum size of initial groups. Stop the process when all initial group sizes are smaller than this parameter (or come from only one study). Default is 10. 
# (4) Minimal overlap requered with the initial 2 NMF programs to be considered for the intersection process.

### Function: 
# (1) For each NMF, calculate how many other NMFs have intersection >=10
# (2) Pick the NMF which has the largest number of intersection NMFs in 1 
# (3) Combine the genes with the highest overlapping NMF (to be explained in detail)
# (4) Create a new intersection matrix relating to the genes in 3, only with NMFs with at least 1 gene overlap to NMF in 1  
# (5) Add NMF with highest intersection (at least 10)
# (6) Create new gene list after adding new NMF
# (7) Stop process when no new NMF can be added. This defines one MP. 
# (8) Repeat 1
# (9) Stop when all groups in 1 have less than 5 NMFs (or are all coming from 1 study?) 

### Use these two from compareNMF_pan_cancer:

nmf_intersect_ccle_KEEP     <-  nmf_intersect_ccle
nmf_programs_sig_ccle_KEEP  <-  nmf_programs_sig_ccle

### Parameters (later change to function form):
Min_intersect_initial <- 10    # the minimal intersection cutoff for defining the Founder NMF program of a cluster
Min_intersect_cluster <- 10    # the minimal intersection cuttof for adding a new NMF to the forming cluster 
Min_group_size        <- 5     # the minimal group size to consider for defining the Founder_NMF of a MP 

Sorted_intersection       <-  sort(apply(nmf_intersect_ccle , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)

Cluster_list <- list()   ### Every entry contains the NMFs of a chosec cluster
k <- 1
Curr_cluster <- c()
MP_list      <- list()

while (Sorted_intersection[1]>Min_group_size) {   ### CHECK!
  
  Curr_cluster <- c(Curr_cluster , names(Sorted_intersection[1]))
  
  ### intersection between all remaining NMFs and Genes in MP 
  Genes_MP                   <- nmf_programs_sig_ccle[,names(Sorted_intersection[1])] # initial genes are those in the first NMF. Genes_MP always has only 50 genes consisting of the current MP
  nmf_programs_sig_ccle      <- nmf_programs_sig_ccle[,-match(names(Sorted_intersection[1]) , colnames(nmf_programs_sig_ccle))]  # remove selected NMF
  Intersection_with_Genes_MP <- sort(apply(nmf_programs_sig_ccle, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
  NMF_history                <- Genes_MP  # has all genes in all NMFs in the current cluster, for newly defining Genes_MP after adding a new NMF 
  
  ### Create gene list - composed of intersecting genes in descending order + genes with highest NMF scores to add up to 50 genes. Update Curr_cluster each time
  
  while ( Intersection_with_Genes_MP[1] >= Min_intersect_cluster) {   ### Define current cluster 
    
    Curr_cluster  <- c(Curr_cluster , names(Intersection_with_Genes_MP)[1])
    
    # Genes_MP_temp   <- sort(table(c(Genes_MP , nmf_programs_sig_ccle[,names(Intersection_with_Genes_MP)[1]])), decreasing = TRUE)
    Genes_MP_temp   <- sort(table(c(NMF_history , nmf_programs_sig_ccle[,names(Intersection_with_Genes_MP)[1]])), decreasing = TRUE)   ## Genes_MP is newly defined each time according to all NMFs in the current cluster 
    Genes_at_border <- Genes_MP_temp[which(Genes_MP_temp == Genes_MP_temp[50])]   ### genes with overlap equal to the 50th gene
    
    if (length(Genes_at_border)>1){
      ### Sort last genes in Genes_at_border according to maximal NMF gene scores
      ### Run over all NMF programs in Curr_cluster and extract NMF scores for each gene
      Genes_curr_NMF_score <- c()
      for (i in Curr_cluster) {
        curr_study           <- paste( strsplit(i , "[.]")[[1]][1 : which(strsplit(i , "[.]")[[1]] == "RDS")]   , collapse = "."  )
        Q                    <- Genes_nmf_w_basis[[curr_study]][ match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]])))[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))]   ,i] 
        names(Q)             <- names(Genes_at_border[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))])  ### sometimes when adding genes the names do not appear 
        Genes_curr_NMF_score <- c(Genes_curr_NMF_score,  Q )
      }
      Genes_curr_NMF_score_sort <- sort(Genes_curr_NMF_score , decreasing = TRUE)
      Genes_curr_NMF_score_sort <- Genes_curr_NMF_score_sort[unique(names(Genes_curr_NMF_score_sort))]   ## take only the maximal score of each gene - which is the first entry after sorting
      
      Genes_MP_temp <- c(names(Genes_MP_temp[which(Genes_MP_temp > Genes_MP_temp[50])]) , names(Genes_curr_NMF_score_sort))
      
    } else {
      Genes_MP_temp <- names(Genes_MP_temp)[1:50] 
    }
    
    NMF_history   <- c(NMF_history , nmf_programs_sig_ccle[,names(Intersection_with_Genes_MP)[1]]) 
    Genes_MP <- Genes_MP_temp[1:50]
    
    nmf_programs_sig_ccle      <- nmf_programs_sig_ccle[,-match(names(Intersection_with_Genes_MP)[1] , colnames(nmf_programs_sig_ccle))]  # remove selected NMF
    
    Intersection_with_Genes_MP <- sort(apply(nmf_programs_sig_ccle, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) # intersection between all other NMFs and Genes_MP  
    
  }
  
  Cluster_list[[paste0("Cluster_",k)]] <- Curr_cluster
  MP_list[[paste0("MP_",k)]]           <- Genes_MP
  k <- k+1
  
  
  nmf_intersect_ccle        <- nmf_intersect_ccle[-match(Curr_cluster,rownames(nmf_intersect_ccle) ) , -match(Curr_cluster,colnames(nmf_intersect_ccle) ) ]  # remove current chosen cluster
  
  Sorted_intersection       <-  sort(apply(nmf_intersect_ccle , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)   ### Sort intersection of remaining NMFs not included in any of the previous clusters
  
  Curr_cluster <- c()
  print(dim(nmf_intersect_ccle)[2])
}



saveRDS(Cluster_list    , file = "Cluster_list_raw_New_Batch.RDS")
saveRDS(MP_list         , file = "MP_list_raw_New_Batch.RDS")




### Remove MPs coming from a single study

num_unique_study <- function (x) {   

  N1       <- as.data.frame(apply(as.data.frame(   x   ) , 1 , function(x)  paste(strsplit(x, "[_]")[[1]][ 1 : which(strsplit(x, "[_]")[[1]] == "scRNAseq"  |   strsplit(x, "[_]")[[1]] == "snRNAseq")  ] , collapse = "_")  )) ; colnames(N1) <- "V1"; N1 <- as.vector(N1$V1)
  
  N2       <- length(unique(N1))
  
  return(N2)
}

unique_study_vec <- as.vector(unlist(   lapply(Cluster_list, function(x) num_unique_study(x))   ))

Cluster_list <- Cluster_list[which(unique_study_vec > 1)]
MP_list      <- MP_list[which(unique_study_vec > 1)]

names(Cluster_list)  <-  paste0("Cluster_" , 1:length(Cluster_list) )
names(MP_list)       <-  paste0("MP_"      , 1:length(MP_list) )



### Also remove MPs that contained NMFs from one single study, and only a single NMF from another study 

num_unique_study_2 <- function (x) {   
  
  N1       <- as.data.frame(apply(as.data.frame(   x   ) , 1 , function(x)  paste(strsplit(x, "[_]")[[1]][ 1 : which(strsplit(x, "[_]")[[1]] == "scRNAseq"  |   strsplit(x, "[_]")[[1]] == "snRNAseq")  ] , collapse = "_")  )) ; colnames(N1) <- "V1"; N1 <- as.vector(N1$V1)
  
  N2       <- length(unique(N1))
  
  N3       <- 0
  
  if (N2 == 2 & (length(which(table(N1) == 1)) > 0) ){
    N3 <- 1
  }
  
  return(N3)
}

unique_study_vec <- as.vector(unlist(   lapply(Cluster_list, function(x) num_unique_study_2(x))   ))

Cluster_list <- Cluster_list[-which(unique_study_vec > 0)]
MP_list      <- MP_list[-which(unique_study_vec > 0)]

names(Cluster_list)  <-  paste0("Cluster_" , 1:length(Cluster_list) )
names(MP_list)       <-  paste0("MP_"      , 1:length(MP_list) )


### Remove MPs we annotated as LQ or miss-annotations

inds_rm <- c(5,7,11,17,25,30,36, 44,58,68, 74) # LQ, LQ, Macro, T-cells, Endo, LQ, LQ, LQ, miss-annotation, LQ, LQ

Cluster_list <- Cluster_list[-inds_rm]
MP_list      <- MP_list[-inds_rm]

names(Cluster_list)  <-  paste0("Cluster_" , 1:length(Cluster_list) )
names(MP_list)       <-  paste0("MP_"      , 1:length(MP_list) )


### Remove small MPs with less than 5 NMFs

inds_rm <- which(as.vector(unlist(lapply(Cluster_list, function(x) length(x)))) < 5)

Cluster_list <- Cluster_list[-inds_rm]
MP_list      <- MP_list[-inds_rm]

names(Cluster_list)  <-  paste0("Cluster_" , 1:length(Cluster_list) )
names(MP_list)       <-  paste0("MP_"      , 1:length(MP_list) )



#### Sort MPs so that they fit the previous order, as in the hallmarks heatmap
#### Put the new MPs next to the MPs in the relevant hallmark

Order_MPs <- c( 1,3,20,22,16,
                2,7,14, 50,
                5,26,32,29,64,
                8,9,28,41,25,46,42,
                6,36, 60,
                4,47,
                11, 66, 10, 
                17, 62, 31,
                38, 13,
                18 ,48 ,
                15, 19, 24, 52, 53,  44,
                21,27,39 , 33, 23, 37 , 51,
                35,40, 56, 49, 59, 55, 65,
                30, 54, 12 , 43, 57  ,61, 67 , 34, 63 , 58 , 45
                
)

MP_names <- c("Cell Cycle - G2/M","Cell Cycle - G1/S","Cell Cylce HMG-rich","Chromatin","Cell cycle single-nucleus",
              "Stress 1" , "Hypoxia", "Stress (in vitro)", "Stress 2" ,
              "Proteasomal degradation", "Unfolded protein response 1", "Protein maturation", "Translation initiation", "Unfolded protein response 2",
              "EMT-I","EMT-II","EMT-III","EMT-IV", "EMT-V", "EMT-VI", "MES (glioma)",
              "Interferon/MHC-II (I)" , "Interferon/MHC-II (II)", "Interferon/MHC-II (III)",
              "Epithelial Senescence", "Epithelial Senescence (HNSCC)",
              "MYC", "NRF2 targets", "P53",
              "Respiration 1", "Respiration 2" , "Respiration (HNSCC)",
              "Secreted I" , "Secreted II",
              "Cilia 1" , "Cilia 2",
              "Astrocytes" , " NPC (glioma)" , "Oligo Progenitor" , "Oligo normal" , "NPC/OPC" ,  "Glioma single-nucleus",
              "PDAC-classical","Alveolar","Skin-pigmentation" , "CRC stemness", "Colon-related" ,"Androgen response (prostate)", "Complement and coagulation (liver)",
              "RBCs", "Platelet-activation", "Hemato-related-I", "IG", "Hemato-related-II", "Hemato-related-III", "Hemato-related-IV",
              "Glutathione", "Metal-response", "PDAC-related 1" , "PDAC-related 2" , "PDAC-related 3" , "PDAC-related 4" , "PDAC-related 5", "Adherens", "Cholesterol Homeostasis", "Unassigned 1", "Unassigned 2"
              )



Cluster_list <- Cluster_list[Order_MPs]
MP_list      <- MP_list[Order_MPs]

names(Cluster_list)  <-  MP_names
names(MP_list)       <-  MP_names


saveRDS(Cluster_list , file = "~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/Cluster_list_Final_New_Batch.RDS")
saveRDS(MP_list      , file = "~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/MP_list_Final_New_Batch.RDS")


### Now sort genes by abundance per cluster 

MP_list_sorted         <- list()
MP_list_sorted_above_1 <- list()
MP_list_sorted_above_2 <- list()


for (I in 1:length(Cluster_list)){
  
  ## unique studies 
  N1 <- as.data.frame(apply(as.data.frame(   Cluster_list[[I]]   ) , 1 , function(x)  paste(strsplit(x, "[_]")[[1]][ 1 : which(strsplit(x, "[_]")[[1]] == "scRNAseq"  |   strsplit(x, "[_]")[[1]] == "snRNAseq")  ] , collapse = "_")  )) ; colnames(N1) <- "V1"; N1 <- as.vector(N1$V1)
  N2 <- unique(N1)
  
  gene_vector <- c()
  
  for (J in 1:length(N2)) { ## run across studies, and collect unique gene from the NMF programs in each 
    
    curr_NMFs <- Cluster_list[[I]][which(N1 == N2[J])]
    R         <- match(curr_NMFs , colnames(nmf_programs_sig_ccle_KEEP))
    
    gene_vector <- c(gene_vector , unique(as.vector(nmf_programs_sig_ccle_KEEP[,R])))
    
  }
  
  sorted_gene_vector <- sort(table(gene_vector) , decreasing = T)
  R2                 <- match( names(sorted_gene_vector) , MP_list[[I]]  ) 
  
  MP_genes_sorted          <- names(sorted_gene_vector)[which(!is.na(R2))]
  MP_genes_sorted_above_1  <- names(sorted_gene_vector)[which(!is.na(R2) & sorted_gene_vector > 1)]
  MP_genes_sorted_above_2  <- names(sorted_gene_vector)[which(!is.na(R2) & sorted_gene_vector > 2)]
  
  
  MP_list_sorted[[names(MP_list)[I]]]          <- MP_genes_sorted
  MP_list_sorted_above_1[[names(MP_list)[I]]]  <- MP_genes_sorted_above_1
  MP_list_sorted_above_2[[names(MP_list)[I]]]  <- MP_genes_sorted_above_2
  
}




#### SAVE MPs : 

MP <-  do.call(cbind, MP_list_sorted)
colnames(MP) <- paste0("MP",1:dim(MP)[2], " ", colnames(MP) )

write.csv(MP, file = "Meta_Programs_generated_automatically_New_Batch.csv")



## Prepare table with NMFs per cluster

abbrev_name <- function (x) {   
  
  N1       <- as.data.frame(apply(as.data.frame(   x   ) , 1 , function(x)  paste(strsplit(x, "[_]")[[1]][ 1 : (which(strsplit(x, "[_]")[[1]] == "rank4") - 1)  ] , collapse = "_")  )) ; colnames(N1) <- "V1"; N1 <- as.vector(N1$V1)
  
  return(N1)
}

Cluster_list_abbrev <- lapply(Cluster_list, function(x) abbrev_name(x))

n_max <- max(as.vector(unlist(lapply(Cluster_list_abbrev, function(x) length(x) ))))

df <- as.data.frame(matrix(data = "" , nrow = n_max , ncol = length(Cluster_list_abbrev)))

for (i in 1:dim(df)[2]) {
  
  df[1:length(Cluster_list_abbrev[[i]]) , i] <- Cluster_list_abbrev[[i]]
  
}


write.csv(df , file = "~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/New_batch_NMFs_per_cluster.csv")








### Plot Jaccard similarity between new MPs and pan-cancer MPs of this cell type  

Cancer_MPs    <- import("~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/MP_by_cell_type/New/Cancer.xlsx") 


intersect_MPs_Jaccard  <- apply(Cancer_MPs , 2, function(x) apply(MP , 2, function(y) length(intersect(x,y)) / length(unique(c(x,y)))  )) 


intersect_MPs_Jaccard_meltI <- reshape2::melt(intersect_MPs_Jaccard) 


P1 <- ggplot(data = intersect_MPs_Jaccard_meltI, aes(x=Var1, y=Var2, fill=100*value, color=100*value)) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 10), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  # theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  # theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  # theme(axis.title.x = element_blank())+
  # theme(axis.title.y = element_blank())+
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))+
  xlab("New Meta Programs")+
  ylab("Previous Meta Programs")


path1 <- paste0("~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/")
file1 <- paste0("MP_comparison.pdf")

ggsave(  filename = file1,
         path = path1,
         width = 12, 
         height = 8,
         device='pdf',
         P1,
         dpi=700)





#### *****  Sort Jaccard similarity plot according to new clusters:

inds_sorted <- c()

for (j in 1:length(Cluster_list)){
  
  inds_sorted <- c(inds_sorted , match(Cluster_list[[j]] , colnames(nmf_intersect_ccle_KEEP)))
  
}
# inds_new <- c(inds_sorted   ,   which(is.na( match(1:dim(nmf_intersect_ccle_KEEP)[2],inds_sorted)))) ### combine inds in clusters with inds without clusters
inds_new   <- inds_sorted


# plot re-ordered similarity matrix heatmap     
nmf_intersect_meltI_ccle_NEW <- reshape2::melt(nmf_intersect_ccle_KEEP[inds_new,inds_new]) 

#### flip matrix so cell cycle is on top 
A = apply(nmf_intersect_ccle_KEEP[inds_new,inds_new], 2, rev)
A1 = reshape2::melt(t(A))


# minor bounds around each MP
limits2 <- c(1)
for (i in 1:length(Cluster_list) ){ # major groups
  
  new_limit        <- limits2[length(limits2)] + length(unlist(Cluster_list[[i]])) - 1
  
  limits2 <- c(limits2 , new_limit , new_limit+1) ## add new limit and begining of next one
  
}
limits2 <- limits2[-length(limits2)]

x1 <- limits2[seq(from=1 , to =length(limits2) , by =2)]
x2 <- limits2[seq(from=2 , to =length(limits2) , by =2)]
y1<-x1; y2<-x2
f=data.frame(x1=x1, x2=x2, y1=y1, y2=y2)



P <- ggplot(data = A1, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1)) +
  # geom_rect(data = f, aes(x = NULL,y = NULL, xmin=x1, xmax=x2, ymin=y1, ymax=y2) , color="red", linetype="dashed", fill = NA , alpha=0.5 , size = 0.3)
  geom_rect(data = f, aes(x = NULL,y = NULL, xmin=x1, xmax=x2, ymin=(length(inds_new)-y2+1), ymax=(length(inds_new)-y1+1)) , color="red", linetype="dashed", fill = NA , alpha=0.5 , size = 0.3)  #  size = 0.25
  


path1 <- paste0("~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/")
file1 <- paste0("MP_cluster_heatmap.pdf")

ggsave(  filename = file1,
         path = path1,
         width = 12, 
         height = 8,
         device='pdf',
         P,
         dpi=700)


path1 <- paste0("~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/")
file1 <- paste0("MP_cluster_heatmap.jpeg")

ggsave(  filename = file1,
         path = path1,
         width = 12, 
         height = 8,
         device='jpeg',
         P,
         dpi=300)





###### calculate proportion of snRNAseq samples per MP
# pure single nucleous studies : "Griffiths_2021_Breast_snRNAseq" , "Hwang_2022_PDAC_scRNAseq", "Jansky_2021_Neuroblastoma_snRNAseq" , "Kim_2018_TNBC_snRNAseq"  , "Unpublished_Glioma_H3G34R_snRNAseq"
# mixed studies : "Biermann_2022_Melanoma_scRNAseq" , "Nath_2021_Ovarian_scRNAseq" , "Wang_2019_Glioblastoma_scRNAseq"


nuc_seq_proportion <- c()
mixed_samples_num  <- c()
total_samples_num  <- c()
nuc_seq_num        <- c()

for (I in 1:length(Cluster_list)){
  
  x       <- Cluster_list[[I]]
  N1      <- as.data.frame(apply(as.data.frame(   x   ) , 1 , function(x)  paste(strsplit(x, "[_]")[[1]][ 1 : which(strsplit(x, "[_]")[[1]] == "scRNAseq"  |   strsplit(x, "[_]")[[1]] == "snRNAseq")  ] , collapse = "_")  )) ; colnames(N1) <- "V1"; N1 <- as.vector(N1$V1)
  inds_rm <- which(N1 == "Biermann_2022_Melanoma_scRNAseq" | N1 == "Nath_2021_Ovarian_scRNAseq" | N1 == "Wang_2019_Glioblastoma_scRNAseq") 
  if (length(inds_rm) > 0){
    N1 <- N1[-inds_rm]
  }
  
  inds_single_nuc <- which(N1 == "Griffiths_2021_Breast_snRNAseq" | N1 == "Hwang_2022_PDAC_scRNAseq" | N1 == "Jansky_2021_Neuroblastoma_snRNAseq" | N1 == "Kim_2018_TNBC_snRNAseq" | N1 == "Unpublished_Glioma_H3G34R_snRNAseq")
  
  nuc_seq_proportion <- c(nuc_seq_proportion , length(inds_single_nuc) / length(N1))
  mixed_samples_num  <- c(mixed_samples_num  , length(inds_rm))
  total_samples_num  <- c(total_samples_num  , length(N1))
  nuc_seq_num        <- c(nuc_seq_num        , length(inds_single_nuc))
  
}

df <- data.frame(MP = names(Cluster_list) , samples_removed = mixed_samples_num , total_samples_num = total_samples_num , nuc_seq_num = nuc_seq_num , NucSeq_proportion = nuc_seq_proportion)

write.csv(df , file = "~/../avishaig/R_WD/Pan_cancer/My_scripts/Pan_Cancer_NMF/nuc_seq_porportion.csv")
