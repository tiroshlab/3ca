# -------------------------------------------------------------------------------------------
# Program for depicting the MP distribution of a given cell types across and within samples 
# ------------------------------------------------------------------------------------------- 

# My_study  = a list in which each element represents a sample in the study; the elements are expression matrices of a given cell type (i.e. cell_type below) in CPM or TPM units. Row names in each matrix should be gene symbols, and column names the cell IDs. 
# cell_type = one of the following cell types : "Cancer" , "Endothelial" , "Epithelial" , "Fibroblasts" , "Macrophages" , "CD4_T_cells", "CD8_T_cells" , "B_cells". Here "Cancer" represents malignant cells
# MinGenes  = the min number of MP genes that should exist in the study (default is 25)
# MinScore  = the minimal score for assigning a cell to a MP (default is 1). A cell is assigned to the MP with the maximal score, given that it exceeded MinScore. If the maximal score is below MinScore, the cell is unassigned. 
# MinCells  = the minimal % of cells that should be assigned to a MP in order to account for the MP (default is 0.05)
# MP_list   = a list of the MPs for each cell_type

### Output: 
# (a) Pie chart of the MP distribution across the study
# (b) Bar plot per tumor with the % of cells per MP
# (c) Expression heatmap per tumor 



# ------------------------------------------------------------------------------------------- 


library(ggplot2)
library(scalop)  # see https://rdrr.io/github/jlaffy/scalop/man/sigScores.html
library(gridExtra)
library(ggpubr)



### Define the following:
# My_study  =
# cell_type = 
# MP_list   = readRDS("MP_list.RDS")   # can be found in the Github repository
# heatCols  = readRDS("heatCols.RDS")  # can be found in the Github repository


Score_cells <- function(L,
                        cell_type,
                        MinGenes = 25,
                        MinScore = 1,
                        MinCells = 0.05,
                        MP_list
                        ) 
  {
  
  MP_list <- MP_list[[cell_type]]
  
  L <- lapply(L, function(x) log2((x/10) + 1) ) # Log normalize before scoring 
  
  MP_scores_per_sample    <- lapply(L, function(x)  scalop::sigScores(x, sigs = MP_list, conserved.genes = MinGenes/50 )  )
  
  remove_cells <- function(x,MinScore){         # remove cells whose max score was below MinScore
    max_score <- apply(x, 1, function(y) max(y))
    cells_rm  <- which(max_score < MinScore)
    if (length(cells_rm) > 0 ){
      x <- x[-cells_rm , ]
    }
    return(x)
  }
  MP_scores_per_sample <- lapply(MP_scores_per_sample, function(x) remove_cells(x,MinScore))
  
  Assign_MP_per_cell   <- lapply(MP_scores_per_sample, function(x) apply(x, 1, function(y) colnames(x)[which(y==max(y))] ) ) 
  
  filter_cells <- function(x,MinCells){         # remove MP that were assassin to less than MinCells in each sample
    MP_frequency <- as.numeric(ave(x, x, FUN = length))
    MP_rm        <- which(MP_frequency/length(x) < MinCells)  # MPs to be removed
    if (length(MP_rm)>0){
      x <- x[-MP_rm]
    }
    return(x)
  }
  Assign_MP_per_cell_filtered <- lapply(Assign_MP_per_cell, function(x) filter_cells(x,MinCells)) 
  
  return(Assign_MP_per_cell_filtered)
  
}


Assign_MP_per_cell_filtered <- Score_cells(L = My_study , cell_type = "Cancer" , MP_list=MP_list)



#### Part A: pie chart of MPs across the whole study
MinCells <- 0.05 
df       <- data.frame(MPs = names(table(unlist(Assign_MP_per_cell_filtered)))  ,  frequency = as.vector(table(unlist(Assign_MP_per_cell_filtered))) )
MPs_rm   <- which(df$frequency/sum(df$frequency) < MinCells)  ## remove also MPs that were assassin to less than MinCells from the total  
if (length(MPs_rm > 0)){
  df$MPs[MPs_rm] <- paste0("Other MPs (<", MinCells*100, "%)")
}


blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )


pie <- ggplot(df, aes(x="", y=frequency, fill=MPs)) + 
       geom_bar(width = 1, stat = "identity") + 
       coord_polar("y", start=0)+
       blank_theme +
       theme(axis.text.x=element_blank())
  

pie


#### Part B: Bar plot with MP distribution per sample

sample_num <- length(Assign_MP_per_cell_filtered)
nCol       <- floor(sqrt(sample_num))
MPs        <- unique(unlist(Assign_MP_per_cell_filtered))    

### For arranging the full MP names as legend in the plot
MP_list <- MP_list[[cell_type]]
MPs     <- MPs[sort(match(MPs, names(MP_list)),index.return = T)$ix]

Assign_MP_per_cell_filtered_abbrev <- lapply(Assign_MP_per_cell_filtered, function(x) apply(as.data.frame(x) , 1 , function(x)  strsplit(x, "[ ]")[[1]][1]) )

df_list    <- lapply(Assign_MP_per_cell_filtered_abbrev, function(x) data.frame(MPs = names(table(x))  ,  frequency = as.vector(table(x))/sum(as.vector(table(x))) ))
df_list2   <- lapply(Assign_MP_per_cell_filtered, function(x) data.frame(MPs = names(table(x))  ,  frequency = as.vector(table(x))/sum(as.vector(table(x))) ))

P <- lapply(seq_along(df_list), function(I) 
     ggplot(data=df_list[[I]], aes(x=reorder(MPs,-frequency), y=frequency)) +
     geom_bar(stat="identity")+
     theme(axis.title.x=element_blank())+
     ggtitle(names(df_list)[I])
  )


P

#### Part C: an expression heatmap for each sample

df_list_ordered  <- lapply(df_list, function(x) x[sort(x$frequency,decreasing = T , index.return = T)$ix,])  ## order MPs in each sample as in the bar plots, in decreasing order
df_list_ordered2 <- lapply(df_list2, function(x) x[sort(x$frequency,decreasing = T , index.return = T)$ix,])  ## order MPs in each sample as in the bar plots, in decreasing order


L1        <- lapply(seq_along(My_study), function(I) My_study[[I]][ , names(Assign_MP_per_cell_filtered[[I]])]) ## extract relevant cells that were assigned to an MP (at least 5%)
names(L1) <- names(My_study)

L1 <- lapply(L1, function(x) log2((x/10) + 1) ) ## Log normalize 
L1 <- lapply(L1, function(x) (x - rowMeans(x)) ) ## center

L2 <- lapply(seq_along(L1), function(I) L1[[I]][ , sort(match(Assign_MP_per_cell_filtered_abbrev[[I]],df_list_ordered[[I]]$MPs), index.return = T)$ix] )  # sort cells according to MPs in df_list_ordered
names(L2) <- names(L1)

MP_genes_per_sample <- lapply(df_list_ordered2, function(x) unlist(MP_list[x$MPs]))  ### MP genes sorted by the MP order in df_list_ordered2
  
L_plot        <- lapply(seq_along(L2), function(I) L2[[I]][match(MP_genes_per_sample[[I]] , rownames(L2[[I]]))[!is.na(match(MP_genes_per_sample[[I]] , rownames(L2[[I]])))] , ])
names(L_plot) <- names(L2)

color.scheme <- colorRampPalette(c(heatCols))(n=333)

Plot_heatmap <- function(M){ 
  
  M_new2        <- M
  M_new2        <- apply(M_new2, 2, rev)
  M_meltII      <-  reshape2::melt(t(M_new2)) 
  M_meltII$Var2 <- factor(M_meltII$Var2)
  
  G <- ggplot(data = M_meltII, aes(x=Var1, y=Var2, fill=value, color=value)) + 
    geom_tile() + 
    scale_color_gradient2(limits=c(-4,4), low=color.scheme[1:111],  mid =color.scheme[112:222], high = color.scheme[223:333], midpoint = 0, oob=squish, name=NULL) +                                
    scale_fill_gradient2(limits=c(-4,4), low=color.scheme[1:111],  mid =color.scheme[112:222], high = color.scheme[223:333], midpoint = 0, oob=squish, name=NULL)  +
    theme(  panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 8), axis.title = element_text(size = 8), legend.title = element_text(size=8), legend.text = element_text(size = 8), legend.text.align = 0.5, legend.justification = "bottom" ) + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank() )
  
  return(G)
}

P1 <- lapply(L_plot, function(x) Plot_heatmap(x))

P_plot <- c()
for (i in 1:length(P)){
  P_plot <- c(P_plot , P[i] , P1[i])
}


P_plot <- do.call("ggarrange", c(P_plot, ncol=2, nrow = 3))

lapply(P_plot, function(x) 
  
  annotate_figure(x, left   = text_grob(paste(MPs , collapse = "\n"), 
                                        color = "red", face = "bold", size = 10))
  
)



