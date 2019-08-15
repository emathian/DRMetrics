---
title: "Vignette_Mesomics"
author: "Emilie Mathian"
date: "6/19/2019"
output:
  html_document:
    toc: true
    toc_depth: 2
    theme: united
---

# What are `Dimensionality_reduction_comparisons` goals ?

`Dimensionality_reduction_comparisons` is a set of function designed to compared dimensionality reduction methods. In order to do this we focus on several indexes, metrics, and statistics :

* Neighborhood preservation : Because distances are no longer meaningful in high dimensional space, we focus on neighborhood preservation between the original space and the porjection spaces.This problematic gather two metrics :

  + Centrality preservation : In order to evaluate the centrality preservation ($C$) between the low dimensional space $D^l$ and the high dimensional one $D^h$ \cite{martins2015explaining}, we used the following formula: 
$$C^d_k(j) = \sum_{1\leq i \leq N} k-\rho^d_i(j)$$
For the scale $k$, and for dimension $d$, $j$'s centrality is defined as the sum of differences between $k$ and $\rho^d_i(j)$, which is the rank of $j$ in the $k-$neighborhood of $i$. Like this high $C$ values represent central points and reciprocally. This metric applied on the projection could be used as a visual tool to evaluate the consistency of the chosen $k$ level. Furthermore, to evaluated how local proprieties are affected by the projection, means of absolutes differences between all $CP^l$ and $CP^h$ for each $k$ were calculated such as $CP_k = \frac{1}{N} \sum_{i = 1}^N (|C^l_i - C^h_i|)$. Finally non parametric tests were realised to evaluated, if $CP_k$  distributions, obtained on real data, are the same than the ones expected under the random hypothesis, and to evaluate if the observed distributions are equal for the different DR methods.

  + Sequence diffence view : The preservation of the $k-$neighborhood between $D^l$ and $D^h$ was calculated according the Sequence Difference view metric ($SD$) \cite{martins2015explaining}. For each point $i$ at the level $k$, $SD$'s formula is: 
$$SD_k(i) = \frac{1}{2} \sum_{j \in V^l_k(i)}(k-\rho^l_i(j)).|\rho^l_i(j)-\rho^h_i(j)|+  \frac{1}{2} \sum_{j \in V^h_k(i)}(k-\rho^h_i(j)).|\rho^l_i(j)-\rho^h_i(j)| $$
where $\rho^d_i(j)$ was previously defined and $V^d_k(i)$ is the $k-$neighborhood of $i$ in the $d$ dimension. This metric penalizes the non preservation of $j$'s rank between $D^l$ and $D^h$, and weights this penalization according $ij$'s  closeness. As previously, non parametric tests were done on the distributions means $SD$ values by level $k$, to tests if observed distributions are the same than the ones expected under the random hypothesis, and to evaluated if the distributions resulting of different DR methods are equal.  

* The consevation variables' spatial autocorreltation :
  + Moran Index : Under the hypothesis that proximal points on the projection share a similar molecular profile, spatial auto-correlations were measured according Moran's Index metric such as:
$$I = \frac{n \sum_{i=1}^n \sum_{j=1}^n W_{ij}(x_i - \bar{x})(x_j - \bar{x})}{\sum_{i=1}^n \sum_{j=1}^n W_{ij} \sum_{i=1}^n (x_i - \bar{x})^2}$$
where $W_{ij}$ is the spatial weight between the the points $i$ and $j$, which belong to the set $n$, $x_i$ is the value of the variable in the unit $i$, and $\bar{x}$ the mean of $x$. The binary matrix $W$ was built according the $k-$nearest neighbor method (KNN), like this the metric is still meaningful in the high dimensional space, since it relies on a relative definition of distances. This observation implies that global spatial auto-correlations are measured for high $k$ values and that $k$ must be must smaller than $n$. Significance of features' Moran's Index  were tested according the Monte Carlo procedure to avoid the Gaussian hypothesis. Finally to compare DR methods according their performances to preserve  features' the spatial distribution, Moran's Index were calculated for several variables and multiple comparison tests were performed.\\

# Installation 

```{r install1 , message=FALSE}
source('../Dimensionality_reduction_comparison_metrics.R')
```
```{r install2, message=FALSE}
library(ade4)
library(umap)
library(dplyr)
library(data.table)
library(rspatial)
library(gridExtra)
library(ggpubr)
library(reshape2)
```

```{r install3, eval=TRUE, echo=TRUE}
sessionInfo()
```
This scpit includes the required library and `CP.R`, `SEQ_DIFF.R` and `MORAN_I.R`.

# Example : Mesomics 

## Importation of data 

```{r import1, eval=TRUE, echo=TRUE}
gene_expr_df <- read.table("gene_expr_mesomics_sample284_genes_7146.txt", header = T)
samples_att <- read.table("Attributes_PCA_f1.csv", header = T, sep=";")
gene_expr_df[1:5,1:4]
head(samples_att[1:5, 1:4])
```
The `gene_expr_df` data frame gathers the expression profile of 284 malignant pleural mesothelioma, with the normalized expression of 7145 genes. The `samples_att` contains the main clinical and physiologic features of these 284 samples.


## Dimension reductions :  methods and projections 

Projection of samples on to the first two principal components.
```{r RD1, eval=TRUE, echo=FALSE}
theme_set(theme_bw())

Histo_type_df <- data.frame("Sample_ID" = as.character(samples_att$sample), "Type" = samples_att$Type )
gene_expr_df_for_PCA <- gene_expr_df[ ,2:dim(gene_expr_df)[2]]
rownames(gene_expr_df_for_PCA) <- gene_expr_df[ ,1]
acp_fig1 <- dudi.pca( gene_expr_df_for_PCA,center = T , scale = F , scannf = F , nf=2) #
acp_fig1_li_df =  as.data.frame(acp_fig1$li)
acp_fig1_li_df  = setDT(acp_fig1_li_df , keep.rownames = TRUE)[]
colnames(acp_fig1_li_df )[1]<-'Sample_ID'
acp_fig1_li_type_df = merge(acp_fig1_li_df , Histo_type_df , by="Sample_ID")

p1 <- ggplot(acp_fig1_li_type_df, aes(x=Axis1 *-1, y=Axis2,  color=Type)) +  geom_point()+
        scale_colour_manual(values=c('#FF9300','#000000','#418F8E', '#FC2B2D'  ))
p1 <- p1 +  labs(title="PCA projection",  y="PC2", x="PC1") +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                                             plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),
                                                             plot.caption =element_text(size=10, color="#17202A", hjust=0.5), 
                                                             axis.title.x=element_text(size=12, face="bold"),  # X axis title
                                                             axis.title.y=element_text(size=12, face="bold"),  # Y axis title
                                                             axis.text.x=element_text(size=12),  # X axis text
                                                             axis.text.y=element_text(size=12))  # Y axis text

```


Projection of samples according UMAP method, with the default parameter (`min_dist` = 0.1, `n_neighbors` = 15).

```{r RD1.2, eval=TRUE, echo=FALSE}
gene_expr_umap_df <- gene_expr_df[order(Histo_type_df$Sample_ID),]
gene_expr_umap_df <- gene_expr_umap_df[,-1] 
gene_expr_umap_df <- apply(gene_expr_umap_df, 2, as.numeric)
umap1 = umap(gene_expr_umap_df)
umap1_res_df <-  data.frame("Sample_ID" = Histo_type_df$Sample_ID, "Axis1" = umap1$layout[, 1], "Axis2" = umap1$layout[, 2],   "Type" = Histo_type_df$Type)
  
p2 <- ggplot(umap1_res_df, aes(x=Axis1, y=Axis2,  color=Type)) +  geom_point()+
        scale_colour_manual(values=c('#FF9300','#000000','#418F8E', '#FC2B2D'  ))
p2 <- p2 +  labs(title="UMAP projection (default parameters)",  y="x", x="y") +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                                             plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),
                                                             plot.caption =element_text(size=10, color="#17202A", hjust=0.5), 
                                                             axis.title.x=element_text(size=12, face="bold"),  # X axis title
                                                             axis.title.y=element_text(size=12, face="bold"),  # Y axis title
                                                             axis.text.x=element_text(size=12),  # X axis text
                                                             axis.text.y=element_text(size=12))  # Y axis text
p1 
p2
```

We follow the same process to generate UMAP projections modifying the `min_dist`parameter such as $min\_dist = [0.1, 0.3, 0.5, 0.7, 0.9]$.

```{r RD2_LNEN, eval=TRUE, echo=FALSE}

umap_md01 = umap(gene_expr_umap_df, min_dist = 0.1)
umap_md01_res_df <-  data.frame("Sample_ID" = Histo_type_df$Sample_ID, "Axis1" = umap_md01$layout[, 1], "Axis2" = umap_md01$layout[, 2])

umap_md03 = umap(gene_expr_umap_df, min_dist = 0.3)
umap_md03_res_df <-  data.frame("Sample_ID" = Histo_type_df$Sample_ID, "Axis1" = umap_md03$layout[, 1], "Axis2" = umap_md03$layout[, 2])

umap_md05 = umap(gene_expr_umap_df, min_dist = 0.5)
umap_md05_res_df <-  data.frame("Sample_ID" = Histo_type_df$Sample_ID, "Axis1" = umap_md05$layout[, 1], "Axis2" = umap_md05$layout[, 2])

umap_md07 = umap(gene_expr_umap_df, min_dist = 0.7)
umap_md07_res_df <-  data.frame("Sample_ID" = Histo_type_df$Sample_ID, "Axis1" = umap_md07$layout[, 1], "Axis2" = umap_md07$layout[, 2])

umap_md09 = umap(gene_expr_umap_df, min_dist = 0.9)
umap_md09_res_df <-  data.frame("Sample_ID" = Histo_type_df$Sample_ID, "Axis1" = umap_md09$layout[, 1], "Axis2" = umap_md09$layout[, 2])

```


## Centrality preservation 

`CP_main` function could be use to generate main analysis on centrality. In this example centrality preservation is computed to compare the PCA projection and UMAP projections.

```{r CP1, eval=FALSE, echo=TRUE}

List_projection <- list(data.frame(acp_fig1_li_df), data.frame(umap_md01_res_df), data.frame(umap_md03_res_df), data.frame(umap_md05_res_df), data.frame(umap_md07_res_df), data.frame(umap_md09_res_df))

gene_expr_df_filter <- merge(gene_expr_df, acp_fig1_li_df[,1], by = "Sample_ID")

Main_CP_res <- CP_main(l_data = List_projection , list_K = seq(from= 1, to = 284, by = 5) , dataRef = gene_expr_df_filter , colnames_res_df = c("pca", "umap_md=0.1", "umap_md=0.3", "umap_md=0.5","umap_md=0.7", "umap_md=0.9") , filename = NULL , graphics = TRUE, stats = TRUE)
saveRDS(Main_CP_res, "Main_CP_res.rds" )

```
```{r CP_res, eval = TRUE, echo=FALSE}
Main_CP_res <- readRDS("Main_CP_res.rds")
```

```{r CP2, eval = TRUE, echo=TRUE}
theme_set(theme_bw())
Main_CP_res$Paired_wilocoxon_test
Main_CP_res$Graphic
```

The mean ranks of the means of the absolutes differences between CP values by k level ($DCP_k$), differs significantly between PCA and UMAP. Furthermore $DCP_k$ values differs for all different values of `min-dist`, this statement was statistically confirmed, even if it is graphically not visible. Indeed `min_dist` produces a constant effect on CP values which implies that even small differences allow to rejet the null hypothesis of wilcoxon tests. According to this index the dimensionality reduction method that respect the most the centrality is the PCA, since $DCP_k$ values are lower. 

**Function : ** `CP_graph_by_k`

The previously graphic could be obtained using the following function; and a logarithmic scale could be choose.

```{r CP4, eval=TRUE, echo=TRUE}
CP_df_res = Main_CP_res[[1]][, 1:(dim(Main_CP_res[[1]])[2]-1 )]
CP_ref_df_res  = Main_CP_res[[1]][1:2]
CP_ref_df_res = cbind(CP_ref_df_res, Main_CP_res[[1]]$REF)
CP_graph_by_k(CP_df_res,  CP_ref_df_res, Names=NULL, list_col=NULL, log=TRUE)
```


**Function : ** `CP_calcul`
If users just want to compute the data frame of CP values for different k levels the function `CP_calcul` should be used.

```{r Calcul, eval=TRUE, echo=TRUE}
CP_calcul_test <- CP_calcul(acp_fig1_li_df , list_K = c(100,150) )
head(CP_calcul_test)
```

### Significance test

If the previous statistic tests allow to know if the mean ranks of $DCP_k$ values are equal between the different methods. We have to prove if $DC_Pk$ values calculated using dimensionality reduction methods are better than those expected if projection methods put points randomly. In order to test this hypothesis, for $n$ repetitions we permute the 2D coordinates of samples, then $CP$ values where calculated. The means distribution of $CP$ values obtained on random data was calculated. Finally a wilcoxon test was effectued to compare the mean ranks of $DCP_k$ values obtained on the real projection coordinates and the random projection coordinates.

```{r Permut1, eval=FALSE, echo=TRUE}
permut_test = CP_permutation_test(data = data.frame(umap_md01_res_df), data_ref = gene_expr_df_filter, seq(1,284,10),  n=40, graph = TRUE)
saveRDS(permut_test, "permut_test_meso_0708.rds")
```



```{r Permut2, eval=TRUE, echo=FALSE}
permut_test2 = readRDS("permut_test_meso_0708.rds")
permut_test2[[1]]
permut_test2[[2]]
```

The means rank of $DCP_k$ values resulting from the PCA projection are statically different from the one obtained on random data.

The significance test could be realized for a level $k$ using a monte carlo procedure using the `CP_monte_carlo` function.
```{r Permut3, eval=FALSE, echo=TRUE}
CP_UMAP_MC_stat <- CP_monte_carlo(data = data.frame(umap_md01_res_df), data_ref = gene_expr_df_filter, k_val =30, n=100)
saveRDS(CP_UMAP_MC_stat,"CP_UMAP_MC_stat.rds" )
```

```{r Permut4, eval=TRUE, echo=FALSE}
CP_UMAP_MC_stat <- readRDS("CP_UMAP_MC_stat.rds")
```

```{r Permut5, eval=TRUE, echo=TRUE}
CP_UMAP_MC_stat
```

### Centrality map 

An useful representation associated to the centrality preservation, could be done with two plot. The first one is a plot of centrality values projected on the 2 dimensional space, and the second one is the projection of centrality values calculated on the high dimensional space and projected on the 2 dimensional space. The first plot just show if the chosen $k$ level is consistent and the second one show the preservation. These representations could be done using 
`CP_map` function.

```{r CP_MAP_PCA, eval=TRUE, echo=FALSE}
CP_Map_calcul_pca <- CP_calcul(acp_fig1_li_df , list_K = c(100) )
CP_K_PCA = CP_Map_calcul_pca[which(  CP_Map_calcul_pca$K == 100),]
CP_K_PCA = data_frame("Sample_ID" = CP_K_PCA$Sample_ID, "K" = CP_K_PCA$K, "CP"= scale(CP_K_PCA$CP))
CP_map(CP_K_PCA , acp_fig1_li_df, list(100), Title = 'CP values projected on the PCA layout')

```

```{r CPMAP2, eval=TRUE, echo=TRUE}
CP_Map_calcul_R <- CP_calcul(as.data.table(gene_expr_df) , list_K = c(100) )
CP_K_R = CP_Map_calcul_R[which( CP_Map_calcul_R$K == 100),]
CP_K_R  = data_frame("Sample_ID" = CP_K_R$Sample_ID, "K" = CP_K_R$K, "CP"= scale(CP_K_R$CP))
CP_map(CP_K_R , acp_fig1_li_df, list(100), Title = 'CP values calculated on genomics data and projected on the PCA layout')

```

We could repeat this process according the UMAP projection.


```{r CP_MAP_UMAP, eval=TRUE, echo=FALSE}
CP_Map_calcul_umap<- CP_calcul(as.data.table(umap_md07_res_df) , list_K = c(100) )
CP_K_Umap = CP_Map_calcul_umap[which(  CP_Map_calcul_umap$K == 100),]
CP_K_Umap = data_frame("Sample_ID" = CP_K_Umap$Sample_ID, "K" = CP_K_Umap$K, "CP"= scale(CP_K_Umap$CP))
CP_map(CP_K_Umap  , umap_md07_res_df, list(100), Title = 'CP values projected on the UMAP layout')
```

```{r CPMAP_UMAP2, eval=TRUE, echo=TRUE}
CP_Map_calcul_R <- CP_calcul(as.data.table(gene_expr_df) , list_K = c(100) )
CP_K_R = CP_Map_calcul_R[which( CP_Map_calcul_R$K == 100),]
CP_K_R  = data_frame("Sample_ID" = CP_K_R$Sample_ID, "K" = CP_K_R$K, "CP"= scale(CP_K_R$CP))
CP_map(CP_K_R , umap_md07_res_df, list(100), Title = 'CP values calculated on genomics data and projected on the UMAP layout')

```


## Sequence difference view :

The sequence difference (SD) view metric is a metric that takes into account :

* The number of true neighbors that a projected point has,

* and the potentially reordering of its neighborhood.

### Main function :

The main function could be use to calculate SD values, to get graphical representation and to obtain statistics.

```{r SQ1, eval=FALSE, echo=TRUE}

List_projection <- list(data.frame(acp_fig1_li_df),data.frame(umap_md01_res_df), data.frame(umap_md03_res_df), data.frame(umap_md05_res_df), data.frame(umap_md07_res_df), data.frame(umap_md09_res_df)) # , 

gene_expr_df_filter <- merge(gene_expr_df, acp_fig1_li_df[,1], by = "Sample_ID")


Main_SQ_res <- Seq_main(l_data = List_projection , dataRef = gene_expr_df_filter, listK = seq(from= 1, to = 280, by = 25), colnames_res_df = c("pca",  "umap_md=0.1", "umap_md=0.3", "umap_md=0.5","umap_md=0.7","umap_md=0.9")  , filename = NULL , graphics = TRUE, stats = TRUE) #``
saveRDS(Main_SQ_res, "Main_SQ_res.rds")
```

For small k values, we cannot observe differences, a logarithmic scale should be more appropriated.
For large k values the PCA seems to be more conservative, since SD values are lower.

```{r sq_res, eval = TRUE, echo=FALSE}
Main_SQ_res <- readRDS("Main_SQ_res.rds")
```

```{r sq_res2, eval = TRUE, echo=FALSE}
Main_SQ_res$pWT
Main_SQ_res$graphics
```


**Zoom in on small k values : **

```{r SQ_smallk, eval=FALSE, echo=TRUE}

List_projection <- list(data.frame(acp_fig1_li_df),data.frame(umap_md01_res_df), data.frame(umap_md03_res_df), data.frame(umap_md05_res_df), data.frame(umap_md07_res_df), data.frame(umap_md09_res_df)) # , 

gene_expr_df_filter <- merge(gene_expr_df, acp_fig1_li_df[,1], by = "Sample_ID")


Main_SQ_res3 <- Seq_main(l_data = List_projection , dataRef = gene_expr_df_filter, listK = seq(from= 1, to = 40, by = 5), colnames_res_df = c("pca",  "umap_md=0.1", "umap_md=0.3", "umap_md=0.5","umap_md=0.7","umap_md=0.9")  , filename = NULL , graphics = FALSE, stats = TRUE) #

saveRDS(Main_SQ_res3, "Main_sq_res_smallK_meso.rds")
```


**Statistic tests : **

```{r SQ_smallk3, eval = TRUE, echo=FALSE}
Main_SQ_res3 <- readRDS("Main_sq_res_smallK_meso.rds")
Main_SQ_res3$pWT
```

We could reprensent the data on a logarithmic scale using the function `Seq_graph_by_k`.
```{r sqcalcul, eval = TRUE, echo=FALSE}
Seq_graph_by_k_MESO <- Seq_graph_by_k(data_Seq = Main_SQ_res3[[1]], Names=NULL, list_col=NULL, data_diff_mean_K = NULL, log = TRUE)
``` 

Using a logarithmic scale we observed than for small k values, UMAP is more conservative, although it is not statiscally significatif.

### Significance test

Finally for each DR technic we have to check if sequence difference values are lower than those expected if the low dimensional representation were random. In order to do this as previously we permute samples' coordinates and we calculate the sequence difference values. After n simulation we calculated the mean distribution and we realize a wilcoxon test with an unilateral hypothesis.

```{r sq_significance_LNEN, eval = FALSE, echo=TRUE}
seq_permut <- seq_permutation_test(data = data.frame(umap_md07_res_df), data_ref = gene_expr_df_filter, list_K = seq(1,284,10), n=4, graph = TRUE)
saveRDS(seq_permut,"seq_permut_meso_2506.rds" )
```

```{r sq_significance2, eval = TRUE, echo=FALSE}
seq_permut <- readRDS("seq_permut_meso_2506.rds")
```

```{r sq_significance3, eval = TRUE, echo=TRUE}
seq_permut
```

On this graph the red curve represents the sequence difference values  resulting from the UMAP projection, and the gray ones result from sequence difference values calculated on randomized low dimensional coordinates. Because UMAP sequence difference values are lower than the ones that are expected on random data, we could confirm that UMAP preserve point's neighborhood. This observation is statiscally confirm by the Wilcoxon test. 

## Spatial Autocorrelation 
Dimensionnality reduction techniques allow to extract the key feataures. Like this it is inetresting to check if correlations are preserved. The Moran Index allow to measure spatial autocorrelation, this statistic is defined according the notion of neighborhood and so for different k level.

### Calculation of Moran Indexes

For different variables, we could use the main function associated to Moran index, to calculated this statistic, generated statistics and graphics.
```{r Moran_I1, eval = TRUE, echo=TRUE, warning=FALSE}
spatial_att_2 <- data.frame("Sample_ID" = as.character(samples_att$sample),"VISTA" = samples_att$VISTA , "VEGFR1" =  samples_att$VEGFR1,  "PDL1" =  samples_att$PDL1, "B_cells" = samples_att$B.Cells)

List_coords_2 <- list('pca' = data.frame(acp_fig1_li_df),'umap_md=01'= data.frame(umap_md01_res_df), 'umap_md=03'=data.frame(umap_md03_res_df), 'umap_md=05' =data.frame(umap_md05_res_df), 'umap_md=07'=data.frame(umap_md07_res_df), 'umap_md=09'=data.frame(umap_md09_res_df), 'gene_expr_data'=data.frame(gene_expr_df))

MI_meso <- moran_I_main(l_coords_data = List_coords_2  , spatial_data =spatial_att_2, listK= seq(50), nsim = 10, Stat=FALSE, Graph = TRUE, methods_name = c('pca','umap_md=01','umap_md=03','umap_md=05', 'umap_md=07','umap_md=09','gene_expr_data'))


```

### Representation of Moran Indexes by level k
Then we could calculate Moran indexes for all k values and for different variables, like for the differents DR technics we get the following graphics : 

```{r Miran_I_plotK, eval = TRUE, echo=FALSE}
List_coords_2 <- list('pca' = data.frame(acp_fig1_li_df),'umap_md=01'= data.frame(umap_md01_res_df),'umap_md=09'=data.frame(umap_md09_res_df), 'gene_expr_data'=data.frame(gene_expr_df))

MI_meso <- moran_I_main(l_coords_data = List_coords_2  , spatial_data =spatial_att_2, listK= seq(50), nsim = 10, Stat=FALSE, Graph = TRUE, methods_name = c('pca','umap_md=01','umap_md=09','gene_expr_data'))

```


```{r Miran_I, eval = TRUE, echo=FALSE}
MI_meso_by_k = moran_I_scatter_plot_by_k(data = as.array(MI_meso)[[1]], Xlab = NULL, Ylab=NULL, Title= NULL)

```
Globally the different technics tend to maintain the same spatial auto-correlations than those observe in the high dimensional space. We also observe that the chosen `min_dist` parameter could influenced the  spatial auto-correlations.

# Example : Lung NeuroEndocrin Neoplasm (LNEN)

**For the MESOMICS dataset a continuum of molocular profiles appeared, like this neither the PCA nor UMAP create clusters. At the opposite for the LNEN dataset revelant molecular group have been identified. This second example allow to show one of the PCA limitation. The PCA is based on a linear transformation of data to find a basis which preserves the maximal amount of variability. On a two dimensional projection no more than three equidistant group could be display, that evidence justify the UMAP use. Nevertheless it is important to keep in mind that if we could represent more than three groups in a two dimensional projection large distances are no longer meaningful.**

## Importation of data 

* Importation of the genes expression data frame, that includes the genes which explained more than 50% of the variation.
* Importation of the features data frame, which gathers the clinical, the histopathological features ... etc

```{r import1_lnen, eval=TRUE, echo=TRUE}
gene_expr_type_LNEN_df <- read.table("gene_expr_LNEN_SCLC_LCNEC.tsv", header = F)
samples_att_LNEN <- read.table("Attributes_LNEN_SCLC.tsv", header = T, sep = '\t')
Histo_6A <- read.table("Histo6A.tsv",  header = T)
```

## Dimensionality reduction

As previously we are going to compare the PCA and the UMAP dimensionality reduction technics : 

```{r RD1_LNEN, eval=TRUE, echo=FALSE}

Clusters_LNEN_df <- data.frame("Sample_ID" = as.character(samples_att_LNEN$Sample_ID), "Clusters" = samples_att_LNEN$Cluster_LNET_LCNEC_SCLC )
gene_expr_LNEN_df_for_PCA <- gene_expr_type_LNEN_df[ ,2:dim(gene_expr_type_LNEN_df)[2]]
rownames(gene_expr_LNEN_df_for_PCA) <- gene_expr_type_LNEN_df[ ,1]
acp_fig1 <- dudi.pca( gene_expr_LNEN_df_for_PCA,center = T , scale = F , scannf = F , nf=2) #
acp_fig1_li_df =  as.data.frame(acp_fig1$li)
acp_fig1_li_df  = setDT(acp_fig1_li_df , keep.rownames = TRUE)[]
colnames(acp_fig1_li_df )[1]<-'Sample_ID'
acp_fig1_li_type_df = merge(acp_fig1_li_df , Clusters_LNEN_df , by="Sample_ID")
acp_fig1_li_type_df = merge(acp_fig1_li_type_df, Histo_6A, by="Sample_ID")
p1 <- ggplot(acp_fig1_li_type_df, aes(x=Axis1 *-1, y=Axis2,  color=Clusters, shape = Histopathology)) +  geom_point()+
      scale_color_viridis(discrete=TRUE) 
p1 <- p1 +  labs(title="PCA projection",  y="PC2 (9.7%)", x="PC1 (27%)") +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                                             plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),
                                                             plot.caption =element_text(size=10, color="#17202A", hjust=0.5), 
                                                             axis.title.x=element_text(size=12),  # X axis title
                                                             axis.title.y=element_text(size=12),  # Y axis title
                                                             axis.text.x=element_text(size=12),  # X axis text
                                                             axis.text.y=element_text(size=12))  # Y axis text


gene_expr_umap_df <- gene_expr_type_LNEN_df[ ,2:dim(gene_expr_type_LNEN_df)[2]]    
gene_expr_umap_df <-   gene_expr_type_LNEN_df[order(Clusters_LNEN_df$Sample_ID),]
gene_expr_umap_df <- gene_expr_umap_df[,-1] 
gene_expr_umap_df <- apply(gene_expr_umap_df, 2, as.numeric)
umap1 = umap(gene_expr_umap_df)
umap1_res_df <-  data.frame("Sample_ID" = Clusters_LNEN_df$Sample_ID, "Axis1" = umap1$layout[, 1], "Axis2" = umap1$layout[, 2],   "Type" = Clusters_LNEN_df$Clusters)
  
p2 <- ggplot(umap1_res_df, aes(x=Axis1, y=Axis2,  color=Type)) +  geom_point()+scale_color_viridis(discrete=TRUE) 
p2 <- p2 +  labs(title="UMAP projection (default parameters)",  y="x", x="y") +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                                             plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),
                                                             plot.caption =element_text(size=10, color="#17202A", hjust=0.5), 
                                                             axis.title.x=element_text(size=12),  # X axis title
                                                             axis.title.y=element_text(size=12),  # Y axis title
                                                             axis.text.x=element_text(size=12),  # X axis text
                                                             axis.text.y=element_text(size=12))  # Y axis text

umap2 = umap(gene_expr_umap_df, min_dist = 0.9,  n_neighbors =100)
umap2_res_df <-  data.frame("Sample_ID" = Clusters_LNEN_df$Sample_ID, "Axis1" = umap2$layout[, 1], "Axis2" = umap2$layout[, 2],   "Type" = Clusters_LNEN_df$Clusters)
  
p3 <- ggplot(umap2_res_df, aes(x=Axis1, y=Axis2,  color=Type)) +  geom_point()+ scale_color_viridis(discrete=TRUE) 
p3 <- p3 +  labs(title="UMAP projection (default parameters)",  y="x", x="y") +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                                             plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),
                                                             plot.caption =element_text(size=10, color="#17202A", hjust=0.5), 
                                                             axis.title.x=element_text(size=12, face="bold"),  # X axis title
                                                             axis.title.y=element_text(size=12, face="bold"),  # Y axis title
                                                             axis.text.x=element_text(size=12),  # X axis text
                                                             axis.text.y=element_text(size=12))  # Y axis text


ggarrange(p1,p2,p3, nrow = 3, ncol = 1,widths = c(2,2,2))
```


We follow the same process to generate other UMAP projections modifying the `min_dist`parameter such as $min\_dist = [0.1, 0.3, 0.5, 0.7, 0.9]$.

```{r RD2, eval=TRUE, echo=FALSE}

umap_md01 = umap(gene_expr_umap_df, min_dist = 0.1)
umap_md01_res_df <-  data.frame("Sample_ID" = Clusters_LNEN_df$Sample_ID, "Axis1" = umap_md01$layout[, 1], "Axis2" = umap_md01$layout[, 2])


umap_md03 = umap(gene_expr_umap_df, min_dist = 0.3)
umap_md03_res_df <-  data.frame("Sample_ID" = Clusters_LNEN_df$Sample_ID, "Axis1" = umap_md03$layout[, 1], "Axis2" = umap_md03$layout[, 2])

umap_md05 = umap(gene_expr_umap_df, min_dist = 0.5)
umap_md05_res_df <-  data.frame("Sample_ID" = Clusters_LNEN_df$Sample_ID, "Axis1" = umap_md05$layout[, 1], "Axis2" = umap_md05$layout[, 2])

umap_md07 = umap(gene_expr_umap_df, min_dist = 0.7)
umap_md07_res_df <-  data.frame("Sample_ID" = Clusters_LNEN_df$Sample_ID, "Axis1" = umap_md07$layout[, 1], "Axis2" = umap_md07$layout[, 2])

umap_md09 = umap(gene_expr_umap_df, min_dist = 0.9)
umap_md09_res_df <-  data.frame("Sample_ID" = Clusters_LNEN_df$Sample_ID, "Axis1" = umap_md09$layout[, 1], "Axis2" = umap_md09$layout[, 2])
```

## Centrality preservation 

`CP_main` function allows to get the main informations about centrality preservation, according the different DR methods.

```{r CP1_LNEN, eval=FALSE, echo=TRUE}

List_projection <- list(data.frame(acp_fig1_li_df), data.frame(umap_md01_res_df), data.frame(umap_md03_res_df), data.frame(umap_md05_res_df), data.frame(umap_md07_res_df), data.frame(umap_md09_res_df))

colnames(gene_expr_type_LNEN_df)[1] <- "Sample_ID"
gene_expr_df_filter <- merge(gene_expr_type_LNEN_df, acp_fig1_li_df[,1], by = "Sample_ID")

Main_CP_res_LNEN <- CP_main(l_data = List_projection , list_K = seq(from= 1, to = 210, by = 5) , dataRef = gene_expr_df_filter , colnames_res_df = c("pca", "umap_md=0.1", "umap_md=0.3", "umap_md=0.5","umap_md=0.7", "umap_md=0.9") , filename = NULL , graphics = TRUE, stats = TRUE)
saveRDS(Main_CP_res_LNEN, 'Main_CP_res_LNEN.rds')
```


```{r CP1_LNEN2, eval=TRUE, echo=FALSE}
Main_CP_res_LNEN <- readRDS("Main_CP_res_LNEN.rds")
```

```{r CP1_LNEN2_4, eval=TRUE, echo=TRUE}
Main_CP_res_LNEN$Paired_wilocoxon_test
```

```{r CP1_LNEN3, eval=TRUE, echo=TRUE}
Main_CP_res_LNEN$Graphic
```

The means ranks of $DCP_k$ values differ between the PCA and UMAP, this was statically confirmed and graphically evident. As previously the `min_dist` parameter affect the centrality preservation, however for this data set the local structure is better preserve with UMAP setting `min_dist` to
$0.9$ or $0.1$. For $k$ level greater than $100$ the PCA better preserves the structure. Finally on LNEN data set high values of `min_dist` seem to improve the preservation of the global structure. 

**Centrality preservation using the first four eigenvalues of the PCA.**
Using the first two axis of the PCA we cannot display on a two dimensional projection more than three equidistant cluster. Because the LNEN dataset presents 5 clusters it could be interesting to evaluate the this metric according the first four axis of the PCA ($n-1$).

```{r ACP_4D, eval=TRUE, echo=TRUE}

acp_4D <- dudi.pca( gene_expr_LNEN_df_for_PCA,center = T , scale = F , scannf = F , nf=4) #
acp_4D_li_df =  as.data.frame(acp_4D$li)
acp_4D_li_df  = setDT(acp_4D_li_df , keep.rownames = TRUE)[]
colnames(acp_4D_li_df )[1]<-'Sample_ID'
acp_4D_li_type_df = merge(acp_4D_li_df , Clusters_LNEN_df , by="Sample_ID")
acp_4D_li_type_df = merge(acp_4D_li_type_df, Histo_6A, by="Sample_ID")

```

```{r ACP_4D_2, eval=TRUE, echo=FALSE}
p1 <- ggplot(acp_4D_li_type_df, aes(x=Axis1 , y=Axis2,  color=Clusters, shape = Histopathology)) +  geom_point()+scale_color_viridis(discrete=TRUE) 
p1 <- p1 +  labs(title="PCA projection",  y="PC2 (9.7%)", x="PC1 (27%)") +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                                             plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),
                                                             plot.caption =element_text(size=10, color="#17202A", hjust=0.5), 
                                                             axis.title.x=element_text(size=12),  # X axis title
                                                             axis.title.y=element_text(size=12),  # Y axis title
                                                             axis.text.x=element_text(size=12),  # X axis text
                                                             axis.text.y=element_text(size=12))  # Y axis text





p2 <- ggplot(acp_4D_li_type_df, aes(x=Axis1 , y=Axis3,  color=Clusters, shape = Histopathology)) +  geom_point()+scale_color_viridis(discrete=TRUE) 
p2 <- p2 +  labs(title="PCA projection",  y="PC3", x="PC1") +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                                             plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),
                                                             plot.caption =element_text(size=10, color="#17202A", hjust=0.5), 
                                                             axis.title.x=element_text(size=12),  # X axis title
                                                             axis.title.y=element_text(size=12),  # Y axis title
                                                             axis.text.x=element_text(size=12),  # X axis text
                                                             axis.text.y=element_text(size=12))  # Y axis text



p3 <- ggplot(acp_4D_li_type_df, aes(x=Axis2 , y=Axis3,  color=Clusters, shape = Histopathology)) +  geom_point()+scale_color_viridis(discrete=TRUE)
p3 <- p3 +  labs(title="PCA projection",  y="PC3", x="PC2") +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                                             plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),
                                                             plot.caption =element_text(size=10, color="#17202A", hjust=0.5), 
                                                             axis.title.x=element_text(size=12),  # X axis title
                                                             axis.title.y=element_text(size=12),  # Y axis title
                                                             axis.text.x=element_text(size=12),  # X axis text
                                                             axis.text.y=element_text(size=12))  # Y axis text




p4 <- ggplot(acp_4D_li_type_df, aes(x=Axis1 , y=Axis4,  color=Clusters, shape = Histopathology)) +  geom_point()+scale_color_viridis(discrete=TRUE)
p4 <- p4 +  labs(title="PCA projection",  y="PC4", x="PC1") +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                                             plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),
                                                             plot.caption =element_text(size=10, color="#17202A", hjust=0.5), 
                                                             axis.title.x=element_text(size=12),  # X axis title
                                                             axis.title.y=element_text(size=12),  # Y axis title
                                                             axis.text.x=element_text(size=12),  # X axis text
                                                             axis.text.y=element_text(size=12))  # Y axis text



ggarrange(p1,p2,p3,p4, nrow = 2, ncol = 2,widths = c(1,1))
```


```{r CP_LNEN_ACP4D0, eval=FALSE, echo=FALSE}
umap_md09 = umap(gene_expr_umap_df, min_dist = 0.9, n_neighbors = 100)
umap_md09_res_df <-  data.frame("Sample_ID" = Clusters_LNEN_df$Sample_ID, "Axis1" = umap_md09$layout[, 1], "Axis2" = umap_md09$layout[, 2])

umap_md01 = umap(gene_expr_umap_df)
umap_md01_res_df <-  data.frame("Sample_ID" = Clusters_LNEN_df$Sample_ID, "Axis1" = umap_md01$layout[, 1], "Axis2" = umap_md01$layout[, 2])


umap_md01_res_df$Sample_ID <- as.character(umap_md01_res_df$Sample_ID)
acp_fig1_li_df$Sample_ID <- as.character(acp_fig1_li_df$Sample_ID)

colnames(gene_expr_type_LNEN_df)[1] <- "Sample_ID"
gene_expr_type_LNEN_df$Sample_ID <- as.character(gene_expr_type_LNEN_df$Sample_ID)
gene_expr_df_filter <- merge(gene_expr_type_LNEN_df, acp_fig1_li_df[,1], by = "Sample_ID")
```

```{r CP_LNEN_ACP4D, eval=FALSE, echo=TRUE}
List_projection <- list( data.frame(acp_4D_li_df),data.frame(acp_fig1_li_df), data.frame(umap_md01_res_df), data.frame(umap_md09_res_df))

Main_CP_res_LNEN_PCA_4D_vfig <- CP_main(l_data = List_projection , list_K = seq(from= 1, to = 208, by = 5) , dataRef = gene_expr_df_filter , colnames_res_df = c("PCA_4D", "PCA", "UMAP_md01_nn15", "UMAP_md09_nn100") , filename = NULL , graphics = TRUE, stats = FALSE)
saveRDS(Main_CP_res_LNEN_PCA_4D_vfig, "Main_CP_res_LNEN_PCA_4D.rds")
```


```{r CP_LNEN_ACP4D2, eval=TRUE, echo=FALSE}
Main_CP_res_LNEN_PCA_4D <- readRDS("Main_CP_res_LNEN_PCA_4D.rds")
```


```{r CP_LNEN_ACP4D3, eval=TRUE, echo=FALSE}
Main_CP_res_LNEN_PCA_4D$Graphic
```


**Logarithmic representation**


```{r CP4D2, eval=TRUE, echo=FALSE}
CP_df_res = Main_CP_res_LNEN_PCA_4D[[1]][, 1:(dim(Main_CP_res_LNEN_PCA_4D[[1]])[2]-1 )]
CP_ref_df_res  = Main_CP_res_LNEN_PCA_4D[[1]][1:2]
CP_ref_df_res = cbind(CP_ref_df_res, Main_CP_res_LNEN_PCA_4D[[1]]$REF)
CP_graph_by_k(CP_df_res,  CP_ref_df_res, Names=NULL, list_col=NULL, log=TRUE)
```



### Significance test according umap projection 

As previously we check if $CP$ values obtained using UMAP stastically different from those calculated on random data, using the `CP_permutation_test`.

```{r Permut1_LNEN, eval=FALSE, echo=TRUE}
permut_test_LNEN_UMAP <- CP_permutation_test(data = data.frame(umap_md01_res_df), data_ref = data.frame(gene_expr_df_filter), seq(1,210,12),  n=10, graph = TRUE)
saveRDS(permut_test_LNEN_UMAP, "CP_ermut_test_LNEN_UMAP.rds")
```


```{r Permut1_LNEN31, eval=TRUE, echo=FALSE}
permut_test_LNEN_UMAP <- readRDS("CP_ermut_test_LNEN_UMAP.rds")
```

```{r Permut1_LNEN4, eval=TRUE, echo=TRUE}
permut_test_LNEN_UMAP
```

Graphically and statiscally UMAP allows a better presrvation of data structure than ones 
expected on randomised projections.


### Centrality preservation map 

As previously we could represent $CP$ values calculated on the 2 dimensional projection and on the high dimensional space on a 2D projection, using `CP_map` function.

```{r CP_MAP_PCA_LNEN, eval=TRUE, echo=FALSE}
CP_Map_calcul_pca <- CP_calcul(acp_fig1_li_df , list_K = c(140) )
CP_K_PCA = CP_Map_calcul_pca[which(  CP_Map_calcul_pca$K == 140),]
CP_K_PCA = data_frame("Sample_ID" = as.character(CP_K_PCA$Sample_ID), "K" = CP_K_PCA$K, "CP"= scale(CP_K_PCA$CP))
CP_map(CP_K_PCA , acp_fig1_li_df, list(140), Title = 'CP values projected on the PCA layout')

```


```{r CPMAP2_LNEN, eval=TRUE, echo=TRUE}
colnames(gene_expr_type_LNEN_df)[1] <- "Sample_ID"
gene_expr_df_filter <- merge(gene_expr_type_LNEN_df, acp_fig1_li_df[,1], by = "Sample_ID")

CP_Map_calcul_R <- CP_calcul(as.data.table(gene_expr_df_filter) , list_K = c(15) )
CP_K_R = CP_Map_calcul_R[which( CP_Map_calcul_R$K == 15),]
CP_K_R  = data_frame("Sample_ID" = as.character(CP_K_R$Sample_ID), "K" = CP_K_R$K, "CP"= scale(CP_K_R$CP))
CP_map(CP_K_R , acp_fig1_li_df, list(15), Title = 'CP values calculated on genomics data and projected on the PCA layout')

```




**According UMAP layout**

We could repeat this process according the UMAP projection.


```{r CP_MAP_UMAP_LNEN, eval=TRUE, echo=FALSE}
CP_Map_calcul_umap<- CP_calcul(as.data.table(umap_md07_res_df) , list_K = c(140) )
CP_K_Umap = CP_Map_calcul_umap[which(  CP_Map_calcul_umap$K == 140),]
CP_K_Umap = data_frame("Sample_ID" = CP_K_Umap$Sample_ID, "K" = CP_K_Umap$K, "CP"= scale(CP_K_Umap$CP))
CP_map(CP_K_Umap  , umap_md07_res_df, list(140), Title = 'CP values projected on the UMAP layout')
```

```{r CPMAP_UMAP_LNEN2, eval=TRUE, echo=TRUE}
CP_Map_calcul_R <- CP_calcul(as.data.table(gene_expr_df_filter) , list_K = c(140) )
CP_K_R = CP_Map_calcul_R[which( CP_Map_calcul_R$K == 140),]
CP_K_R  = data_frame("Sample_ID" = CP_K_R$Sample_ID, "K" = CP_K_R$K, "CP"= scale(CP_K_R$CP))
CP_map(CP_K_R , umap_md07_res_df, list(140), Title = 'CP values calculated on genomics data and projected on the UMAP layout')

```


## Sequence difference view :

As previously explained once the consistancy of points position is checked we have to measure the preservation of points' neighborhood.

### Main function :

As previously we first could apply the main function to calculate the sequence difference values to have the main graphic, and to obtain statistics that compare the different technics. 

```{r SQ1_LNEN, eval=FALSE, echo=TRUE}

List_projection <- list(data.frame(acp_fig1_li_df),data.frame(umap_md01_res_df), data.frame(umap_md03_res_df), data.frame(umap_md05_res_df), data.frame(umap_md07_res_df), data.frame(umap_md09_res_df)) # , 

colnames(gene_expr_type_LNEN_df)[1] <- "Sample_ID"
gene_expr_df_filter <- merge(gene_expr_type_LNEN_df, acp_fig1_li_df[,1], by = "Sample_ID")

Main_SQ_res_LNEN <- Seq_main(l_data = List_projection , dataRef = gene_expr_df_filter, listK = seq(from= 1, to = 208, by = 5), colnames_res_df = c("pca",  "umap_md=0.1", "umap_md=0.3", "umap_md=0.5","umap_md=0.7","umap_md=0.9")  , filename = NULL , graphics = TRUE, stats = TRUE) #
saveRDS( Main_SQ_res_LNEN,  "SQ_MAIN_LNEN.rds")

```


```{r SQ3_LNEN, eval=TRUE, echo=FALSE}
Main_SQ_res_LNEN <- readRDS("SQ_MAIN_LNEN.rds")
```

```{r SQ3.2_LNEN, eval=TRUE, echo=FALSE}
Main_SQ_res_LNEN$pWT
```
The values that result from the PCA are statiscally different from those obtained using UMAP and this whatever the `min_dist` parameter value.  Furthemore, the `min_dist` parameter significally impacts the sequence difference values since most of the time the null hypothesis is rejected.

```{r SQ4_LNEN3, eval=TRUE, echo=TRUE}
Main_SQ_res_LNEN$graphics
```



**Zoom in on  small k values.**

```{r SQ_smallk_LNEN, eval=FALSE, echo=TRUE}
Main_SQ_res_smallK_LNEN <- Seq_main(l_data = List_projection , dataRef = gene_expr_df_filter, listK = seq(from= 1, to = 40, by = 5), colnames_res_df = c("pca",  "umap_md=0.1", "umap_md=0.3", "umap_md=0.5","umap_md=0.7","umap_md=0.9")  , filename = NULL , graphics = FALSE, stats = TRUE) #
saveRDS(Main_SQ_res_smallK_LNEN, "Main_SQ_res_smallK_LNEN.rds")
```
```{r SQ2_LNEN, eval=TRUE, echo=FALSE}
Main_SQ_res_smallK_LNEN <- readRDS("Main_SQ_res_smallK_LNEN.rds")
```

```{r SQ3_LNEN3_SK, eval=TRUE, echo=TRUE}
Main_SQ_res_smallK_LNEN$pWT
```


According this representation we observe that UMAP is localy more conservative.

**As could previously it could be interesting to evaluate the Sequence difference metric for the first four eigenvectors of the PCA.**

```{r sq_main_analysis4D, eval = FALSE, echo=FALSE, warning=FALSE}
List_projection <- list( data.frame(acp_4D_li_df),data.frame(acp_fig1_li_df),data.frame(umap_md01_res_df), data.frame(umap_md09_res_df))
umap_md01_res_df$Sample_ID <- as.character(umap_md01_res_df$Sample_ID)
acp_fig1_li_df$Sample_ID <- as.character(acp_fig1_li_df$Sample_ID)
colnames(gene_expr_type_LNEN_df)[1] <- "Sample_ID"
gene_expr_type_LNEN_df$Sample_ID <- as.character(gene_expr_type_LNEN_df$Sample_ID)
gene_expr_df_filter <- merge(gene_expr_type_LNEN_df, acp_fig1_li_df[,1], by = "Sample_ID")

Main_SQ_res_4D_VFIG <- Seq_main(l_data = List_projection , dataRef = gene_expr_df_filter, listK = seq(from= 1, to = 200, by = 15), colnames_res_df = c( "PCA_4D", "PCA_2D","umap_md=0.1","umap_md=0.9")  , filename = NULL , graphics = TRUE, stats = TRUE) #
saveRDS(Main_SQ_res_4D_VFIG, "Main_SQ_res_4D_VFIG.rds")
```


```{r sq_main_analysis4D2, eval = TRUE, echo=FALSE, warning=FALSE}
Main_SQ_res_4D_VFIG <- readRDS("Main_SQ_res_4D_VFIG.rds")
Main_SQ_res_4D_VFIG$graphics
Main_SQ_res_4D_VFIG$pWT
```


### Significance test

Finally we have to check if sequence difference values obtained are statistically different from those expected under the random hypothesis.

```{r sq_significance_LNEN1, eval = FALSE, echo=TRUE}
seq_permut_Umap_LNEN <- seq_permutation_test(data = data.frame(umap_md07_res_df), data_ref = gene_expr_df_filter, list_K = seq(1,210,10), n=4, graph = TRUE)
saveRDS(seq_permut_Umap_LNEN ,"seq_permut_Umap_LNEN .rds" )
```


```{r sq_significance_LNEN2, eval = TRUE, echo=FALSE}
seq_permut_Umap_LNEN <- readRDS("seq_permut_Umap_LNEN .rds")

```


```{r sq_significance_LNEN3, eval = TRUE, echo=TRUE}
seq_permut_Umap_LNEN
```

According the previous statistic and the sttistic test UMAP preserve the points' neighborhood. 

## Spatial Autocorrelation 

As previously we could use the Moran's statistic to measure spatial auto-correlations, in order to know if the different dimensionality reduction methods allow to extract the same features.

### Calculation of Moran Indexes

Firstly we could the Moran indexes values for different features and for the different DR methods using the `moran_I_main` function.

```{r Moran_I_LNEN1, eval = TRUE, echo=TRUE, warning= FALSE}
spatial_att <- data.frame("Sample_ID" = as.character(samples_att_LNEN$Sample_ID),"CYP2C8" = samples_att_LNEN$CYP2C8 , "OTP" =  samples_att_LNEN$OTP,  "Neutrophils" =  samples_att_LNEN$Neutrophils)

spatial_att <- spatial_att[complete.cases(spatial_att),] 
spatial_att_kepp_ids <- data_frame("Sample_ID" = spatial_att$Sample_ID)

acp_fig1_li_df <- merge(acp_fig1_li_df, spatial_att_kepp_ids, by = 'Sample_ID' )
umap_md01_res_df <-  merge(umap_md01_res_df , spatial_att_kepp_ids, by = 'Sample_ID' )
umap_md03_res_df <-  merge(umap_md03_res_df , spatial_att_kepp_ids, by = 'Sample_ID' )
umap_md05_res_df <-  merge(umap_md05_res_df , spatial_att_kepp_ids, by = 'Sample_ID' )
umap_md07_res_df <-  merge(umap_md07_res_df , spatial_att_kepp_ids, by = 'Sample_ID' )
umap_md09_res_df <-  merge(umap_md09_res_df , spatial_att_kepp_ids, by = 'Sample_ID' )
gene_expr_df_filter  <- merge(gene_expr_df_filter, spatial_att_kepp_ids, by = 'Sample_ID' )

List_coords <- list('pca' = data.frame(acp_fig1_li_df),'umap_md=01'= data.frame(umap_md01_res_df), 'umap_md=03'=data.frame(umap_md03_res_df), 'umap_md=05' =data.frame(umap_md05_res_df), 'umap_md=07'=data.frame(umap_md07_res_df), 'umap_md=09'=data.frame(umap_md09_res_df), 'gene_expr_data'=data.frame(gene_expr_df_filter))

MI_LNEN <- moran_I_main(l_coords_data = List_coords  , spatial_data =spatial_att, listK= seq(5,200,10), nsim = 10, Stat=FALSE, Graph = TRUE, methods_name = c('pca','umap_md=01','umap_md=03','umap_md=05', 'umap_md=07','umap_md=09','gene_expr_data'))

```

According to these results we observe that the spatial auto-correlations calculated in the low dimensional space are consistant with those obtain in the high dimensional space. For the chosen features it is difficult to affirm if the DR methods tend to increase or decrease the spatial autocorrelations observed in the high dimensional space. Finally, we observe that the choice of `min_dist` parameter highly interfer with the spatial auto-correlations values. 

```{r Moran_I_LNEN2, eval = TRUE, echo= TRUE, warning=FALSE}

umap_md09 = umap(gene_expr_umap_df, min_dist = 0.9, n_neighbors = 100)
umap_md09_res_df <-  data.frame("Sample_ID" = Clusters_LNEN_df$Sample_ID, "Axis1" = umap_md09$layout[, 1], "Axis2" = umap_md09$layout[, 2])
umap_md09_res_df$Sample_ID <- as.character(umap_md09_res_df$Sample_ID)
umap_md09_res_df <- umap_md09_res_df[-which(umap_md09_res_df$Sample_ID == "S02322.R2"| umap_md09_res_df$Sample_ID == "S00716_B"),]

umap_md01 = umap(gene_expr_umap_df)
umap_md01_res_df <-  data.frame("Sample_ID" = Clusters_LNEN_df$Sample_ID, "Axis1" = umap_md01$layout[, 1], "Axis2" = umap_md01$layout[, 2])
umap_md01_res_df$Sample_ID <- as.character(umap_md01_res_df$Sample_ID)
umap_md01_res_df <- umap_md01_res_df[-which(umap_md01_res_df$Sample_ID == "S02322.R2"| umap_md01_res_df$Sample_ID == "S00716_B"),]

acp_4D_li_df <- acp_4D_li_df[-which(acp_4D_li_df$Sample_ID == "S02322.R2" | acp_4D_li_df$Sample_ID == "S00716_B" ),]

acp_fig1_li_df <- acp_fig1_li_df[-which(acp_fig1_li_df$Sample_ID == "S02322.R2" | acp_fig1_li_df$Sample_ID == "S00716_B" ),]


MI_LNEN <- moran_I_main(l_coords_data = List_coords  , spatial_data =spatial_att, listK= seq(5,200,10), nsim = 10, Stat=FALSE, Graph = TRUE, methods_name = c('pca','umap_md=01','umap_md=03','umap_md=05', 'umap_md=07','umap_md=09','gene_expr_data'))

MI_meso_by_k = moran_I_scatter_plot_by_k(data = as.array(MI_LNEN)[[1]], Xlab = NULL, Ylab=NULL, Title= NULL)
```


### Overlap the features that have the highest Moran index according the different methods 

For a given k level we could determined the Moran Indexes values for 
all genes and for different methods (PCA, UMAP, high dimensional data), using the funcion `moran_ranking`. This function allows to calculated Moran indexes of a list of features according *a* level k and It allows to return a data frame containing the results. The function also allows to order the lists of genes according their Moran Index  and to determine the overlapping between the different methods taking into account the first $n$ genes.

```{r MoranRanks1, eval = FALSE, echo= TRUE}
all_genes <- read.table("VST_nosex_TCACLCNECSCLC.txt", header = T)
#colnames(all_genes) 
all_genes <- all_genes[,- which(colnames(all_genes) == "S02322_B")]
colnames(all_genes)[which(colnames(all_genes) == "S02322_A")] = "S02322.R1"

vst_mat <- all_genes[,c(which(colnames(all_genes) %in% as.character(samples_att_LNEN$Sample_ID)))]
rvdm <- apply(vst_mat,1,stats::var)
ntopdm = rev(which( cumsum( sort(rvdm,decreasing = T) )/sum(rvdm) <= 0.5 ))[1] + 1
selectdm <- order(rvdm, decreasing = TRUE)[seq_len(min(ntopdm, length(rvdm)))]
dim(vst_mat)
vst_mat <- vst_mat[selectdm,]
dim(vst_mat)
vst_mat  <- t(vst_mat)

acp_fig1_li_df <- acp_fig1_li_df[order(acp_fig1_li_df$Sample_ID),]
umap_md09_res_df <- umap_md09_res_df[order(umap_md09_res_df$Sample_ID),]
gene_expr_df_filter$Sample_ID <- as.character(gene_expr_df_filter$Sample_ID)
gene_expr_df_filter <- gene_expr_df_filter[order(gene_expr_df_filter$Sample_ID),]

List_coords <- list('pca' = data.frame(acp_fig1_li_df), 'umap_nn=09'=data.frame(umap_md09_res_df), 'gene_expr_data'=data.frame(gene_expr_df_filter))
vst_mat_sample <- cbind("Sample_ID" = acp_fig1_li_df$Sample_ID , vst_mat)
vst_mat_sample <- apply(vst_mat_sample[,2:dim(vst_mat_sample )[2]] , 2, as.numeric)
vst_mat_sample <- as.data.frame(vst_mat_sample)
vst_mat_sample <- cbind('Sample_ID' = as.character(spatial_att$Sample_ID), vst_mat_sample )
vst_mat_sample$Sample_ID <- as.character(vst_mat_sample$Sample_ID)
vst_mat_sample <- vst_mat_sample[order(vst_mat_sample$Sample_ID),]

MI_LNEN_rank <- moran_ranking(l_coords_data = List_coords, spatial_data = vst_mat_sample, K_value =30, methods_name = NULL, N = 10)
```
```{r MoranRanks2, eval = TRUE, echo= FALSE }
MI_LNEN_rank <- readRDS("MI_LNEN_rank.rds")
```


According to the venn table a venn diagram could be built. For $n=1000$ we obtained the following diagram :

[Venn_diagram](moran_I_venn.pdf) 



### Moran Index overlapping between PCA 4D and UMAP nn =121
```{r MoranRanks_PCA4D_UMAP, eval = FALSE, echo= TRUE}
vst50_TCACLCNECSCLC<- read.table("VST_nosex_50pc_TCACLCNECSCLC.txt", header = T)
dim(vst50_TCACLCNECSCLC)

# SELECTION
vst50_TCACLCNECSCLC <- data.frame(t(vst50_TCACLCNECSCLC ))
vst50_TCACLCNECSCLC_designRD <- vst50_TCACLCNECSCLC 
vst50_TCACLCNECSCLC_designRD <- vst50_TCACLCNECSCLC_designRD[- which(rownames(vst50_TCACLCNECSCLC_designRD) == "S00716_A"|rownames(vst50_TCACLCNECSCLC_designRD) == "S02322_B"),]
rownames(vst50_TCACLCNECSCLC_designRD)[which(rownames(vst50_TCACLCNECSCLC_designRD) == "S02322_A")] <- "S02322.R1"

vst50_TCACLCNECSCLC <- setDT(vst50_TCACLCNECSCLC , keep.rownames = TRUE)[]
colnames(vst50_TCACLCNECSCLC)[1] <- "Sample_ID"
vst50_TCACLCNECSCLC <- vst50_TCACLCNECSCLC[-which(vst50_TCACLCNECSCLC$Sample_ID == "S02322_B"| vst50_TCACLCNECSCLC$Sample_ID == "S00716_A"), ]
vst50_TCACLCNECSCLC$Sample_ID[which(vst50_TCACLCNECSCLC$Sample_ID == "S02322_A")] <- "S02322.R1"
dim(vst50_TCACLCNECSCLC)

vst50_TCACLCNECSCLC$Sample_ID <- as.character(vst50_TCACLCNECSCLC$Sample_ID)
vst50_TCACLCNECSCLC <- vst50_TCACLCNECSCLC[order(vst50_TCACLCNECSCLC$Sample_ID),]
acp_4D_li_df <- acp_4D_li_df[order(acp_4D_li_df$Sample_ID ),]

List_coords <- list('PCA 4D' = data.frame(acp_4D_li_df), 'umap_nn=121'=data.frame(umap121_res_df[,1:3]))

MI_LNEN_rank_pca4D_umap <- moran_ranking(l_coords_data = List_coords, spatial_data = data.frame(vst50_TCACLCNECSCLC), K_value =121, methods_name = NULL, N = 1000)

saveRDS(MI_LNEN_rank_pca4D_umap, "MI_LNEN_rank_pca4D_umap.rds")
```



```{r MoranRanks2_PCA4D_UMAP, eval = TRUE, echo= FALSE }
MI_LNEN_rank <- readRDS("MI_LNEN_rank.rds")
```
 
**Finaly we could determine the overlapping for all the different n values** (where $n$ are the number of genes considered). In order to show if these overlapping values are different from those that could be obtained by hazard, we repeat the procedure permuting the lists. Finally we plot the overlapping values in function of $n$ for the real observation for on the simulated data.

```{r MoranRanks4, eval = FALSE, echo= TRUE }
overlap <- c()
for (N in 5:dim(MI_LNEN_rank$MoranIndex )[1]){#
  first_cor_list <- c()
  for (i in 1:dim(MI_LNEN_rank$MoranIndex)[2]){
    sort_df <- MI_LNEN_rank$MoranIndex[order(MI_LNEN_rank$MoranIndex[,i]),]
    sortN_l <- rownames(sort_df)[1:N]
    first_cor_list[[i]] <- sortN_l
  }
  names(first_cor_list) <- c("PCA", "UMAP", "Gene_expr")
  v.table<- length(intersect(  intersect(first_cor_list[[1]],   first_cor_list[[2]]),   first_cor_list[[3]]))
  overlap <- c(overlap, v.table/N )
}

overlap_simul <- unlist(lapply(5:6576, function(i){
  (i *(i/6576)^2) / i
}))
overlap_df <- data.frame("N" = seq(5,6576,1),"Real_Data" = overlap , "Simulation" = overlap_simul )
saveRDS(overlap_df, "overlap_df.rds")
```

```{r MoranRanks4.1, eval = TRUE, echo= FALSE }
overlap_df <- readRDS("overlap_df.rds")
```


```{r MoranRanks6, eval = TRUE, echo= FALSE }

overlap_df_2 <- melt(overlap_df , id=c("N"))

theoric <- unlist(lapply(5:dim(MI_LNEN_rank$MoranIndex )[1], function(i){
  (i *(i/dim(MI_LNEN_rank$MoranIndex )[1])^2) / i
}))

theme_set(theme_bw())
p <- ggplot()
p <- p + geom_point(data = overlap_df_2 , aes(x=N, y=value,  colour=variable))+ 
  scale_colour_manual(values=c("#1f1f99", '#848484'))
p <- p + geom_line(aes(x=seq(5:dim(MI_LNEN_rank$MoranIndex )[1]),  y=theoric))
p <- p +  labs(title="Percentage of overlapping between lists of genes ranked by their Moran Index", 
               y="Overlapping %", x="number of genes", caption = "Lists compared : Genes expression data frame, PCA, UMAP") +theme(plot.title=element_text(size=18, face="bold", color="#17202A", hjust=0.5,lineheight=1.2),  # title
                                                                                                                                        plot.subtitle =element_text(size=13, color="#17202A", hjust=0.5),  # caption
                                                                                                                                        plot.caption =element_text(size=10, color="#17202A", hjust=0.5),  # caption
                                                                                                                                        axis.title.x=element_text(size=12, face="bold"),  # X axis title
                                                                                                                                        axis.title.y=element_text(size=12, face="bold"),  # Y axis title
                                                                                                                                        axis.text.x=element_text(size=12),  # X axis text
                                                                                                                                        axis.text.y=element_text(size=12) ,
                                                                                                                                        legend.text = element_text( size = 16),
                                                                                                                                        legend.title = element_text( face = 'bold',size=18)  )  # Y axis text
print(p) 
```


## Tuning `n_neighbors` 

```{r RD1_NN_0, eval=TRUE, echo=FALSE}

vst50_TCACLCNECSCLC<- read.table("VST_nosex_50pc_TCACLCNECSCLC.txt", header = T)
samples_att_LNEN <- read.table("Attributes_LNEN_SCLC.tsv", header = T, sep = '\t')
Histo_6A <- read.table("Histo6A.tsv",  header = T)
attributes_TCACLCNECSCLC <- read.table("Attributes_fig6A.tsv", header = T, sep = '\t', fill =TRUE)

# SELECTION
vst50_TCACLCNECSCLC <- data.frame(t(vst50_TCACLCNECSCLC ))
vst50_TCACLCNECSCLC_designRD <- vst50_TCACLCNECSCLC # Structure match with the one require for PCA and UMAP
vst50_TCACLCNECSCLC_designRD <- vst50_TCACLCNECSCLC_designRD[- which(rownames(vst50_TCACLCNECSCLC_designRD) == "S00716_A"|rownames(vst50_TCACLCNECSCLC_designRD) == "S02322_B"),]
rownames(vst50_TCACLCNECSCLC_designRD)[which(rownames(vst50_TCACLCNECSCLC_designRD) == "S02322_A")] <- "S02322.R1"

vst50_TCACLCNECSCLC <- setDT(vst50_TCACLCNECSCLC , keep.rownames = TRUE)[]
colnames(vst50_TCACLCNECSCLC)[1] <- "Sample_ID"
vst50_TCACLCNECSCLC <- vst50_TCACLCNECSCLC[-which(vst50_TCACLCNECSCLC$Sample_ID == "S02322_B"| vst50_TCACLCNECSCLC$Sample_ID == "S00716_A"), ]
vst50_TCACLCNECSCLC$Sample_ID[which(vst50_TCACLCNECSCLC$Sample_ID == "S02322_A")] <- "S02322.R1"


# Cleaning Sample ID
samples_att_LNEN$Sample_ID <- as.character(samples_att_LNEN$Sample_ID)
samples_att_LNEN <- samples_att_LNEN[-which(samples_att_LNEN$Sample_ID == "S02322.R2"| samples_att_LNEN$Sample_ID == "S00716_A"),]

# Cluster Data Frame
Clusters_df <- data.frame("Sample_ID" = samples_att_LNEN$Sample_ID, "Clusters" = samples_att_LNEN$Cluster_LNET_LCNEC_SCLC )
Clusters_df$Sample_ID <- as.character(Clusters_df$Sample_ID)

# UMAP md = 0.1 NN= 15
umap1 = umap(vst50_TCACLCNECSCLC_designRD)
umap1_res_df <- as.data.frame(umap1$layout)
umap1_res_df  = setDT(umap1_res_df , keep.rownames = TRUE)[]
colnames(umap1_res_df)[1] <- "Sample_ID"
umap1_res_df <- merge(umap1_res_df, Clusters_df, by="Sample_ID")

umap121 = umap(vst50_TCACLCNECSCLC_designRD, n_neighbors = 121)
umap121_res_df <- as.data.frame(umap121$layout)
umap121_res_df  = setDT(umap121_res_df , keep.rownames = TRUE)[]
colnames(umap121_res_df)[1] <- "Sample_ID"
umap121_res_df <- merge(umap121_res_df, Clusters_df, by="Sample_ID")

umap208 = umap(vst50_TCACLCNECSCLC_designRD, n_neighbors = 208)
umap208_res_df <- as.data.frame(umap208$layout)
umap208_res_df  = setDT(umap208_res_df , keep.rownames = TRUE)[]
colnames(umap208_res_df)[1] <- "Sample_ID"
umap208_res_df <- merge(umap208_res_df, Clusters_df, by="Sample_ID")

############################# SUPPLEMENTARY TABLE LYNNETTE  ########################

attributes_TCACLCNECSCLC$Sample_ID = as.character(attributes_TCACLCNECSCLC$Sample_ID)
attributes_TCACLCNECSCLC <- attributes_TCACLCNECSCLC[-which(attributes_TCACLCNECSCLC$Sample_ID == "S02322.R2"|attributes_TCACLCNECSCLC$Sample_ID == "S00716_A" ),]
Molecular_clusters_df <- read.table("histo_details_lynnette.csv", header = T, sep = ';')


Molecular_clusters <- c()
for (i in 1:dim(attributes_TCACLCNECSCLC)[1]){
  if (attributes_TCACLCNECSCLC$Sample_ID[i] %in% Molecular_clusters_df$SCLC.LCNEC.like){
    Molecular_clusters[i] <- "SCLC/LCNEC-like"
  }
  else if(attributes_TCACLCNECSCLC$Sample_ID[i] %in% Molecular_clusters_df$LCNEC.SCLC.like ){
    Molecular_clusters[i] <- "LCNEC/SCLC-like"
  }
  else if(attributes_TCACLCNECSCLC$Sample_ID[i] %in% Molecular_clusters_df$LCNEC.type.1){
    Molecular_clusters[i] <- "LCNEC/TypeI"
  }
  else if(attributes_TCACLCNECSCLC$Sample_ID[i] %in% Molecular_clusters_df$LCNEC.type.2){
    Molecular_clusters[i] <- "LCNEC/TypeII"
  }
  else{
    Molecular_clusters[i] <- as.character(Clusters_df$Clusters[i])
  }
}

Molecular_clusters[Molecular_clusters=="SCLC"] = "SCLC/SCLC-like"
Molecular_clusters[Molecular_clusters=="CaA1"] = "Carcinoid-A1"
Molecular_clusters[Molecular_clusters=="CaA2"] = "Carcinoid-A2"
Molecular_clusters[Molecular_clusters=="CaB"] = "Carcinoid-B"
Molecular_clusters[Molecular_clusters=="LCNEC"] = "LCNEC/NA"
Molecular_clusters[which(attributes_TCACLCNECSCLC$Histpopathology_simplified == "Supra_carcinoid")] = "Supra_carcinoid"

attributes_TCACLCNECSCLC <- cbind(attributes_TCACLCNECSCLC, "Molecular_clusters" = Molecular_clusters) 
umap1_res_df <-cbind(umap1_res_df, "Molecular_clusters" = Molecular_clusters )
umap208_res_df <-cbind(umap208_res_df, "Molecular_clusters" = Molecular_clusters )
umap121_res_df <-cbind(umap121_res_df, "Molecular_clusters" = Molecular_clusters )



library(RColorBrewer)

p2_standard <- ggplot(umap1_res_df, aes(x=V1, y=-V2,  color=Molecular_clusters )) +  geom_point(size=4, alpha =.8) + scale_color_brewer(palette="Spectral")
p2_standard <- p2_standard +  labs(title="UMAP n == 15", y="dim2", x="dim1")  +
  theme( plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=16),  # X axis title
         axis.title.y=element_text(size=16),  # Y axis title
         axis.text.x=element_text(size=14),  # X axis text
         axis.text.y=element_text(size=14),
         legend.text = element_text(size = 14) ,
         legend.title = element_blank())  # Y axis text

p2_standard





p2_standard_nn200 <- ggplot(umap208_res_df, aes(x=V1, y=-V2,  color=Molecular_clusters )) +  geom_point(size=4, alpha =.8) + scale_color_brewer(palette="Spectral")
p2_standard_nn200 <- p2_standard_nn200 +  labs(title="UMAP nn = 200", y="dim2", x="dim1")  +
  theme( plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=16),  # X axis title
         axis.title.y=element_text(size=16),  # Y axis title
         axis.text.x=element_text(size=14),  # X axis text
         axis.text.y=element_text(size=14),
         legend.text = element_text(size = 14) ,
         legend.title = element_blank())  # Y axis text

p2_standard_nn200


p2_standard_nn121 <- ggplot(umap121_res_df, aes(x=V1, y=-V2,  color=Molecular_clusters )) +  geom_point(size=4, alpha =.8) + scale_color_brewer(palette="Spectral")
p2_standard_nn121 <- p2_standard_nn121 +  labs(title="UMAP nn = 121", y="dim2", x="dim1")  +
  theme( plot.title=element_text(size=16, face="bold", hjust=0.5,lineheight=1.2),  # title
         plot.subtitle =element_text(size=14, hjust=0.5),
         plot.caption =element_text(size=12,  hjust=0.5),
         axis.title.x=element_text(size=16),  # X axis title
         axis.title.y=element_text(size=16),  # Y axis title
         axis.text.x=element_text(size=14),  # X axis text
         axis.text.y=element_text(size=14),
         legend.text = element_text(size = 14) ,
         legend.title = element_blank())  # Y axis text

p2_standard_nn121


```




```{r RD1_NN_1, eval=TRUE, echo=FALSE}

gene_expr_type_LNEN_df <- read.table("gene_expr_LNEN_SCLC_LCNEC.tsv", header = F)
colnames(gene_expr_type_LNEN_df)[1] <- "Sample_ID"
gene_expr_df_filter <- gene_expr_type_LNEN_df


gene_expr_LNEN_df_PCA <- gene_expr_df_filter[,2:dim(gene_expr_df_filter)[2]]
Clusters_LNEN_df_no_rep <- merge(Clusters_LNEN_df, gene_expr_type_LNEN_df, by = "Sample_ID")
rownames(gene_expr_LNEN_df_PCA) <- Clusters_LNEN_df_no_rep$Sample_ID
acp_fig1 <- dudi.pca( gene_expr_LNEN_df_PCA,center = T , scale = F , scannf = F , nf=2) #
acp_fig1_li_df =  as.data.frame(acp_fig1$li)
acp_fig1_li_df  = setDT(acp_fig1_li_df , keep.rownames = TRUE)[]
colnames(acp_fig1_li_df )[1]<-'Sample_ID'
acp_fig1_li_type_df = merge(acp_fig1_li_df ,  Clusters_LNEN_df_no_rep, by="Sample_ID")




acp_4D <- dudi.pca( gene_expr_LNEN_df_PCA,center = T , scale = F , scannf = F , nf=5)
acp_4D_li_df =  as.data.frame(acp_4D$li)
acp_4D_li_df  = setDT(acp_4D_li_df , keep.rownames = TRUE)[]
colnames(acp_4D_li_df )[1]<-'Sample_ID'
acp_4D_li_type_df = merge(acp_4D_li_df ,  Clusters_LNEN_df_no_rep, by="Sample_ID")


```



### Centrality preservation 


```{r CP1_NN_1, eval=FALSE, echo=TRUE}
acp_fig1_li_df$Sample_ID <- as.character(acp_fig1_li_df$Sample_ID)
acp_4D_li_df$Sample_ID <- as.character(acp_4D_li_df$Sample_ID)
umap1_res_df$Sample_ID <- as.character(umap1_res_df$Sample_ID)
umap121_res_df$Sample_ID <- as.character(umap121_res_df$Sample_ID)
umap208_res_df$Sample_ID <- as.character(umap208_res_df$Sample_ID)
List_projection <- list(data.frame(acp_fig1_li_df), data.frame(acp_4D_li_df), data.frame(umap1_res_df[,2:3]), data.frame(umap121_res_df[,2:3]), data.frame(umap208_res_df[,2:3]))

colnames(gene_expr_type_LNEN_df)[1] <- "Sample_ID"
gene_expr_df_filter <- merge(gene_expr_type_LNEN_df, acp_fig1_li_df[,1], by = "Sample_ID")

Main_CP_res_NN <- CP_main(l_data = List_projection , list_K = seq(from= 1, to = 210, by = 5) , dataRef = gene_expr_df_filter , colnames_res_df = c("pca_2D", "pca_4D","umap_md=0.1_nn=15", "umap_md=0.1_nn=122", "umap_md=0.1_nn=208") , filename = NULL , graphics = TRUE, stats = TRUE)

saveRDS(Main_CP_res_NN, "Main_CP_res_NN.rds")
```


```{r CP1_NNR, eval=TRUE, echo=FALSE}
Main_CP_res_NN<- readRDS("Main_CP_res_NN.rds")
```

```{r CP1_LNEN2_2, eval=TRUE, echo=TRUE}
Main_CP_res_NN$Paired_wilocoxon_test
```


```{r CP_cNNGRAPH, eval=TRUE, echo=FALSE}
CP_df_res = Main_CP_res_NN[[1]][, 1:(dim(Main_CP_res_NN[[1]])[2]-1 )]
CP_ref_df_res  =Main_CP_res_NN[[1]][1:2]
CP_ref_df_res = cbind(CP_ref_df_res, Main_CP_res_NN[[1]]$REF)
CP_graph_by_k(CP_df_res,  CP_ref_df_res, Names=NULL, list_col=NULL, log=TRUE)
```


### Sequence difference view
```{r SQ1_consistancy3, eval=FALSE, echo=TRUE}
acp_fig1_li_df <- acp_fig1_li_df[-which(acp_fig1_li_df$Sample_ID == "S00716_A"| acp_fig1_li_df$Sample_ID == "S02322.R2"),]

acp_4D_li_df <- acp_4D_li_df[-which(acp_4D_li_df$Sample_ID == "S00716_A" | acp_4D_li_df$Sample_ID == "S02322.R2"),]
                 
gene_expr_df_filter <- gene_expr_df_filter[-which(gene_expr_df_filter$Sample_ID == "S00716_A"| gene_expr_df_filter$Sample_ID == "S02322.R2"),]
            
List_projection <- list(data.frame(acp_fig1_li_df), data.frame(acp_4D_li_df),data.frame(umap1_res_df[,1:3]), data.frame(umap121_res_df[,1:3]), data.frame(umap208_res_df[,1:3]))
str(List_projection )
str(gene_expr_df_filter)
Main_SQ_res_NN <- Seq_main(l_data = List_projection , dataRef = gene_expr_df_filter, listK = seq(from= 1, to = 208, by = 5), colnames_res_df = c("pca_2D", "pca_4D","umap_md=0.1_nn=15", "umap_md=0.1_nn=121", "umap_md=0.1_nn=208")  , filename = NULL , graphics = TRUE, stats = TRUE) #
saveRDS(Main_SQ_res_NN, "Main_SQ_res_NN.rds")
```





```{r SQ1_consistancy2_C6g, eval=TRUE, echo=FALSE}
Main_SQ_res_NN <- readRDS("Main_SQ_res_NN.rds")

Main_SQ_res_NN$pWT
p1 <- Seq_graph_by_k(data_Seq =Main_SQ_res_NN[[1]], Names=NULL, list_col=NULL, data_diff_mean_K = NULL, log = FALSE)

p2 <- Seq_graph_by_k(data_Seq =Main_SQ_res_NN[[1]], Names=NULL, list_col=NULL, data_diff_mean_K = NULL, log = TRUE)
```

