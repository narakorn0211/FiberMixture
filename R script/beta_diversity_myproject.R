
#Install the packages
install.packages("tidyverse")
install.packages("vegan")
install.packages("devtools")
library(devtools)
devtools::install_github("jbisanz/qiime2R")

#Load the packages. Everyone needs to do this.
library(tidyverse)
library(vegan)
library(qiime2R)


##############################################
#Set UP
#
#These are the things that  we need from Qiime:
#
#metadata.tsv
#core-metrics-results/bray_curtis_pcoa_results.qza
#core-metrics-results/weighted_unifrac_pcoa_results.qza
#core-metrics-results/rarefied_table.qza
#rooted-tree.qza
#taxonomy.qza
#core-metrics-results/evenness_vector.qza
#core-metrics-results/faith_pd_vector.qza
#core-metrics-results/observed_otus_vector.qza
#core-metrics-results/shannon_vector.qza
#
# These files are already in the ANSC516-repo
##############################################

getwd()
###Set your working directory path/to/ANSC516/ANSC-repo/data/moving-pictures
setwd("~/Desktop/ANSC_file/project/R")

list.files()

if(!dir.exists("output"))
  dir.create("output")

#How to load a file into R
metadata2 <- read.delim("metadata.tsv", sep = "\t", header = T, quote = "", stringsAsFactors = F)
metadata2[1,]
metadata2[,1]
# When subsetting, the first number is the row and after the comma is the column
metadata3 <- metadata2[-1,]

#Now the qiime2R method
metadata <- read.delim("metadata.tsv", sep = "\t", header = TRUE)
str(metadata)

row.names(metadata) <- metadata[,1]
row.names(metadata) <- metadata$SampleID
#metadata <- metadata[,-1]
row.names(metadata)

bc_PCoA<-read_qza("bray_curtis_pcoa_results.qza")
wUF <- read_qza("weighted_unifrac_pcoa_results.qza")

fiber <- c("Black", "Blue", "Green", "Gray","Pink", "Navy", "Yellow", "Brown", "Orange", "Purple", "Magenta", "Cyan", "Red", "Darkgreen", "Khaki")

bc_meta <- bc_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2, PC3) %>%
  inner_join(metadata, by = c("SampleID" = "sample.id"))

# Now we are going to make an ordination plot
library(ggplot2)
ggplot(bc_meta, aes(x=PC1, y=PC2, color=fiber)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  xlab("PC1 (32.27%)") +
  ylab("PC2 (22.28%)") +
  scale_color_manual(values=c("Blue", "Black", "Green", "Gray", "Pink", "Navy", "Yellow", "Brown", "Orange", "Purple", "Magenta", "Cyan", "Red", "Darkgreen", "Khaki"), name = "fiber")

# Now we are going to make our code a little more re-usable
fiber_colors <- c("Black", "Blue", "Green", "Gray","Pink", "Navy", "Yellow", "Brown", "Orange", "Purple", "Magenta", "Cyan", "Red", "Darkgreen", "Khaki")
my_column <- "fiber"
#my_column <- "FiberTreatment"

ggplot(bc_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  theme_q2r() +
  facet_grid(~donor) +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=fiber_colors, name = my_column)
ggsave(paste0("output/FB-basic_", my_column,".tiff"), height=3, width=4, device="tiff") # save a PDF 3 inches by 4 inches

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),bc_meta,mean)
colnames(centroids)[1] <- "fiber"

ggplot(bc_meta, aes(x=PC1, y=PC2, color= fiber)) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=fiber_colors, name = fiber)
ggsave(paste0("output/FB-ellipse_", fiber,".pdf"), height=7, width=7, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(bc_meta, aes(x=PC1, y=PC2, color= fiber)) +
  geom_point(aes(shape= rep)) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) 
#scale_color_manual(values=fiber_colors, name = Fiber)
ggsave(paste0("output/FB-ellipse_", fiber,"-rep.pdf"), height=7, width=7, device="pdf") # save a PDF 3 inches by 4 inches

#######NEWONE FOR UPDATE THE NAME COLUMN AS FIBER######
bc_meta$fiber <- bc_meta[[my_column]]
ggplot(bc_meta, aes(x = PC1, y = PC2, color = fiber)) +
  geom_point(aes(shape= rep)) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  #stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*bc_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*bc_PCoA$data$ProportionExplained[2], digits = 2), "%)")) 
#scale_color_manual(values=fiber_colors, name = Fiber)
ggsave(paste0("output/FB-ellipse_", my_column,"-rep.pdf"), height=3, width=4.5, device="pdf")

##################################################################################
## SAME thing but with weighted UniFrac

Wuni_PCoA<-read_qza("weighted_unifrac_pcoa_results.qza")

Wuni_meta <- Wuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "sample.id"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Wuni_meta,mean)

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 2) +
  theme_q2r() +
  theme(
    axis.title = element_text(size = 12),     # Axis titles
    axis.text = element_text(size = 12),      # Tick labels
    legend.title = element_text(size = 10),   # Legend title
    legend.text = element_text(size = 9),    # Legend labels
    strip.text = element_text(size = 14)      # For facets if any
  ) +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 1), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 1), "%)")) +
  scale_color_manual(values= fiber_colors, name = "Fiber")
ggsave(paste0("output/Wuni-ellipse_", my_column,".pdf"), height=5, width=7, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Wuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= rep), size = 1.5) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  theme(
    axis.title = element_text(size = 12),     # Axis titles
    axis.text = element_text(size = 12),      # Tick labels
    legend.title = element_text(size = 10),   # Legend title
    legend.text = element_text(size = 9),    # Legend labels
    strip.text = element_text(size = 12)      # For facets if any
  ) +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Wuni_PCoA$data$ProportionExplained[1], digits = 1.5), "%)")) +
  ylab(paste0("PC2 (", round(100*Wuni_PCoA$data$ProportionExplained[2], digits = 1.5), "%)")) +
  scale_color_manual(values=fiber_colors, name = "Fiber")
ggsave(paste0("output/Wuni-ellipse_", my_column,"-rep.pdf"), height=8, width=7, device="pdf") 

#########################Unweight Uni#####################################
UNWuni_PCoA<-read_qza("unweighted_unifrac_pcoa_results.qza")

UNWuni_meta <- UNWuni_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "sample.id"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),UNWuni_meta,mean)

ggplot(UNWuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  theme(
    axis.title = element_text(size = 12),     # Axis titles
    axis.text = element_text(size = 12),      # Tick labels
    legend.title = element_text(size = 10),   # Legend title
    legend.text = element_text(size = 9),    # Legend labels
    strip.text = element_text(size = 12)      # For facets if any
  ) +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*UNWuni_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*UNWuni_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=fiber_colors, name = "Fiber")
ggsave(paste0("output/UNWuni-ellipse_", my_column,".pdf"), height=5, width=7, device="pdf") 

ggplot(UNWuni_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= rep), size = 1.5) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  theme(
    axis.title = element_text(size = 12),     # Axis titles
    axis.text = element_text(size = 12),      # Tick labels
    legend.title = element_text(size = 10),   # Legend title
    legend.text = element_text(size = 9),    # Legend labels
    strip.text = element_text(size = 12)      # For facets if any
  ) +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*UNWuni_PCoA$data$ProportionExplained[1], digits = 1.5), "%)")) +
  ylab(paste0("PC2 (", round(100*UNWuni_PCoA$data$ProportionExplained[2], digits = 1.5), "%)")) +
  scale_color_manual(values=fiber_colors, name = "Fiber")
ggsave(paste0("output/UNWuni-ellipse_", my_column,"-rep.pdf"), height=8, width=7, device="pdf") # save a PDF 3 inches by 4 inches


##################################################################################
##########Jaccard

Jaccard_PCoA<-read_qza("jaccard_pcoa_results.qza")

Jaccard_meta <- Jaccard_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "sample.id"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Jaccard_meta,mean)

ggplot(Jaccard_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  theme(
    axis.title = element_text(size = 12),     # Axis titles
    axis.text = element_text(size = 12),      # Tick labels
    legend.title = element_text(size = 10),   # Legend title
    legend.text = element_text(size = 9),    # Legend labels
    strip.text = element_text(size = 12)      # For facets if any
  ) +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Jaccard_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Jaccard_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=fiber_colors, name = "Fiber")
ggsave(paste0("output/Jaccard-ellipse_", my_column,".pdf"), height=5, width=7, device="pdf") # save a PDF 3 inches by 4 inches

ggplot(Jaccard_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= rep), size = 1.5) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 1.5) +
  theme_q2r() +
  theme(
    axis.title = element_text(size = 12),     # Axis titles
    axis.text = element_text(size = 12),      # Tick labels
    legend.title = element_text(size = 10),   # Legend title
    legend.text = element_text(size = 9),    # Legend labels
    strip.text = element_text(size = 12)      # For facets if any
  ) +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Jaccard_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Jaccard_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=fiber_colors, name = "Fiber")
ggsave(paste0("output/Jaccard-ellipse_", my_column,"-rep.pdf"), height=8, width=7, device="pdf") # save a PDF 3 inches by 4 inches

#################################################################################
###### Bray Curtis
Bray_PCoA<-read_qza("bray_curtis_pcoa_results.qza")

Bray_meta <- Bray_PCoA$data$Vectors %>%
  select(SampleID, PC1, PC2) %>%
  inner_join(metadata, by = c("SampleID" = "sample.id"))

centroids <- aggregate(cbind(PC1,PC2)~get(my_column),Bray_meta,mean)

ggplot(Bray_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point() + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 3) +
  theme_q2r() +
  theme(
    axis.title = element_text(size = 12),     # Axis titles
    axis.text = element_text(size = 12),      # Tick labels
    legend.title = element_text(size = 10),   # Legend title
    legend.text = element_text(size = 9),    # Legend labels
    strip.text = element_text(size = 12)      # For facets if any
  ) +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Bray_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Bray_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=fiber_colors, name = "Fiber")
ggsave(paste0("output/Bray-ellipse_", my_column,".pdf"), height=5, width=7, device="pdf") 

ggplot(Bray_meta, aes(x=PC1, y=PC2, color=get(my_column))) +
  geom_point(aes(shape= rep), size = 1.5) + #alpha controls transparency and helps when points are overlapping
  #geom_point(data=centroids, size = 1.5) +
  theme_q2r() +
  theme(
    axis.title = element_text(size = 12),     # Axis titles
    axis.text = element_text(size = 12),      # Tick labels
    legend.title = element_text(size = 10),   # Legend title
    legend.text = element_text(size = 9),    # Legend labels
    strip.text = element_text(size = 12)      # For facets if any
  ) +
  stat_ellipse(level = 0.95, type = "t") +
  xlab(paste0("PC1 (", round(100*Bray_PCoA$data$ProportionExplained[1], digits = 2), "%)")) +
  ylab(paste0("PC2 (", round(100*Bray_PCoA$data$ProportionExplained[2], digits = 2), "%)")) +
  scale_color_manual(values=fiber_colors, name = "Fiber")
ggsave(paste0("output/Bray-ellipse_", my_column,"-rep.pdf"), height=8, width=7, device="pdf")
###################################
#Run some PERMANOVAs

bc_dist_mat<-read_qza("bray_curtis_distance_matrix.qza")
bc_dm <- as.matrix(bc_dist_mat$data) 
rownames(bc_dm) == metadata$SampleID ## all these values need to be "TRUE"
metadata_sub <- metadata[match(rownames(bc_dm),metadata$sample.id),]
rownames(bc_dm) == metadata_sub$sample.id ## all these values need to be "TRUE"

metadata_sub <- metadata_sub[!is.na(metadata_sub$fiber), ]
PERMANOVA_out <- adonis2(bc_dm ~ fiber, data = metadata_sub)

head(rownames(bc_dm))
head(metadata$sample.id)

write.table(PERMANOVA_out,"output/Fiber_Adonis_overall.csv",sep=",", row.names = TRUE) 

######################################################################################
##  Pairwise adonis function
##  we can also performe a pairwise comparison with the function 
##  Pairwise Adonis funtion by edro Martinez Arbizu & Sylvain Monteux
##  https://github.com/pmartinezarbizu/pairwiseAdonis/blob/master/pairwiseAdonis/R/pairwise.adonis.R
#######################################################################################

pairwise.adonis2 <- function(x, data, strata = NULL, nperm=999, ... ) {
  
  #describe parent call function 
  ststri <- ifelse(is.null(strata),'Null',strata)
  fostri <- as.character(x)
  #list to store results
  
  #copy model formula
  x1 <- x
  # extract left hand side of formula
  lhs <- x1[[2]]
  # extract factors on right hand side of formula 
  rhs <- x1[[3]]
  # create model.frame matrix  
  x1[[2]] <- NULL   
  rhs.frame <- model.frame(x1, data, drop.unused.levels = TRUE) 
  
  # create unique pairwise combination of factors 
  co <- combn(unique(as.character(rhs.frame[,1])),2)
  
  # create names vector   
  nameres <- c('parent_call')
  for (elem in 1:ncol(co)){
    nameres <- c(nameres,paste(co[1,elem],co[2,elem],sep='_vs_'))
  }
  #create results list  
  res <- vector(mode="list", length=length(nameres))
  names(res) <- nameres
  
  #add parent call to res 
  res['parent_call'] <- list(paste(fostri[2],fostri[1],fostri[3],', strata =',ststri, ', permutations',nperm ))
  
  
  #start iteration trough pairwise combination of factors  
  for(elem in 1:ncol(co)){
    
    #reduce model elements  
    if(inherits(eval(lhs),'dist')){	
      xred <- as.dist(as.matrix(eval(lhs))[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),
                                           rhs.frame[,1] %in% c(co[1,elem],co[2,elem])])
    }else{
      xred <- eval(lhs)[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),]
    }
    
    mdat1 <-  data[rhs.frame[,1] %in% c(co[1,elem],co[2,elem]),] 
    
    # redefine formula
    if(length(rhs) == 1){
      xnew <- as.formula(paste('xred',as.character(rhs),sep='~'))	
    }else{
      xnew <- as.formula(paste('xred' , 
                               paste(rhs[-1],collapse= as.character(rhs[1])),
                               sep='~'))}
    
    #pass new formula to adonis
    if(is.null(strata)){
      ad <- adonis2(xnew,data=mdat1, ... )
    }else{
      perm <- how(nperm = nperm)
      setBlocks(perm) <- with(mdat1, mdat1[,ststri])
      ad <- adonis2(xnew,data=mdat1,permutations = perm, ... )}
    
    res[nameres[elem+1]] <- list(ad[1:5])
  }
  #names(res) <- names  
  class(res) <- c("pwadstrata", "list")
  return(res)
} 

fiber_Pair <- pairwise.adonis2(bc_dm ~ fiber, data = metadata_sub)
write.table(fiber_Pair,"output/Fiber_Adonis_pairwise.csv",sep=",", row.names = TRUE) 

