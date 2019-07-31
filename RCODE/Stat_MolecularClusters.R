##########################################################################################################
##################### STAT Lynnette question distance ####################################################
##########################################################################################################

############################# LIBRAIRIES #####################################
library(ggplot2)
library(umap)
library(RColorBrewer)
library(viridis)
###############################################################################


############################# IMPORTATIONS & PREPROCESSING #####################################
set.seed(1564404882)
umap_coords <- read.table("Coords_umap_nn208.tsv")
colnames(umap_coords) <- c("Sample_ID", "V1", "V2")
umap_coords <- umap_coords[order(umap_coords$Sample_ID),]
AttributesUmap <- read.table("Attributes_UMAP_TCACLCNECSCLC.tsv", header = T, sep= '\t')
AttributesUmap <- AttributesUmap[order(AttributesUmap$Sample_ID),]
###############################################################################

p2 <- ggplot(umap_coords, aes(x=V1, y=V2,  color=AttributesUmap$Molecular_clusters )) +  geom_point(size=4, alpha =.5)+ scale_color_brewer(palette="Spectral")
p2

Molecular_clusters <- AttributesUmap$Molecular_clusters
summary(as.factor(Molecular_clusters))
Mol_clusters_df <- data.frame("Sample_ID" = attributes_TCACLCNECSCLC$Sample_ID, 'Molecular_clusters'=  AttributesUmap$Molecular_clusters)
##################### Remove outliers

umap_coords_type <- cbind(umap_coords , Molecular_clusters) 
umap_coords_type <- umap_coords_type[-which(umap_coords_type$Sample_ID == "S00602"| umap_coords_type$Sample_ID == "S02297")]

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Test 1 : Est ce que les SCLC/LCNEC-like cluster plus avec les LCNEC qu'avec les small cells ?

# Etape 1 : calcul du centroid LCNECT1 LCNECT2 LCNEC_SCLC-like et LCNEC/NA


T1x_umap1_df <- umap_coords_type[which( umap_coords_type$Molecular_clusters== "LCNEC/TypeI"|umap_coords_type$Molecular_clusters== "LCNEC/TypeII"  ),]
T1x_centroid <- apply(T1x_umap1_df[,2:3], 2,mean)

T1_x_ref <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "SCLC/LCNEC-like"),]

dist_x <- unlist(lapply( 1:dim(T1_x_ref)[1], function(i){
  sqrt( (T1_x_ref[i,2] - T1x_centroid[1])^2 + (T1_x_ref[i,3] - T1x_centroid[2])^2 )
}))



T1y_umap1_df <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "SCLC/SCLC-like"),]
T1y_centroid <- apply(T1y_umap1_df[,2:3], 2,mean)

T1_x_ref <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "SCLC/LCNEC-like"),]

dist_y <- unlist(lapply( 1:dim(T1_x_ref)[1], function(i){
  sqrt( (T1_x_ref[i,2] - T1y_centroid[1])^2 + (T1_x_ref[i,3] - T1y_centroid[2])^2 )
}))


wilcox.test(dist_x , dist_y, 'less')

p2 <- ggplot(umap_coords, aes(x=V1, y=V2,  color=Molecular_clusters )) +  geom_point(size=1, alpha =.5)+scale_color_viridis(discrete=TRUE) + geom_point(aes(x =T1y_centroid[1], y = T1y_centroid[2] )) +
  geom_jitter(size=1.3, alpha=0.6, position = position_jitter(width = 1, height = 1))+   geom_point(size=1, alpha =.5)+geom_point(aes(x =T1x_centroid[1], y = T1x_centroid[2]))

p2

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Test 2 : Est ce que les LCNEC/SCLC-like cluster plus avec les LCNEC qu'avec les small cells ?
summary(as.factor(Molecular_clusters))

T2x_umap1_df <- umap_coords_type[which( umap_coords_type$Molecular_clusters== "SCLC/SCLC-like" ),]
T2x_centroid <- apply(T2x_umap1_df[,2:3], 2,mean)

T2_x_ref <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "LCNEC/SCLC-like"),]

T2_dist_x <- unlist(lapply( 1:dim(T2_x_ref)[1], function(i){
  sqrt( (T2_x_ref[i,2] - T2x_centroid[1])^2 + (T2_x_ref[i,3] - T2x_centroid[2])^2 )
}))



T2y_umap1_df <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "LCNEC/TypeI"|umap_coords_type$Molecular_clusters== "LCNEC/TypeII"),]
T2y_centroid <- apply(T2y_umap1_df[,2:3], 2,mean)

T2_x_ref <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "LCNEC/SCLC-like"),]

T2_dist_y <- unlist(lapply( 1:dim(T2_x_ref)[1], function(i){
  sqrt( (T2_x_ref[i,2] - T2y_centroid[1])^2 + (T2_x_ref[i,3] - T2y_centroid[2])^2 )
}))


wilcox.test(T2_dist_x , T2_dist_y, alternative = "less")

p2 <- ggplot(umap_coords, aes(x=V1, y=V2,  color=Molecular_clusters )) +  geom_point(size=1, alpha =.5)+scale_color_viridis(discrete=TRUE) + geom_point(aes(x =T2y_centroid[1], y = T2y_centroid[2] )) +
  geom_jitter(size=1.3, alpha=0.6, position = position_jitter(width = 1, height = 1))+   geom_point(size=1, alpha =.5) + geom_point(aes(x =T2x_centroid[1], y = T2x_centroid[2]))

p2
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Test 3
# Est ce que les SCLC/LCNEC-like sont plus avec les LCNEC-TypeI que les LCNEC-TypeII
summary(as.factor(Molecular_clusters))
T3x_umap1_df <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "LCNEC/TypeII"),]
T3x_centroid <- apply(T3x_umap1_df[,2:3], 2,mean)

T3_x_ref <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "SCLC/LCNEC-like"),]

T3_dist_x <- unlist(lapply( 1:dim(T3_x_ref)[1], function(i){
  sqrt( (T3_x_ref[i,2] - T3x_centroid[1])^2 + (T3_x_ref[i,3] - T3x_centroid[2])^2 )
}))



T3y_umap1_df <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "LCNEC/TypeI"),]
T3y_centroid <- apply(T3y_umap1_df[,2:3], 2,mean)

T3_x_ref <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "SCLC/LCNEC-like"),]

T3_dist_y <- unlist(lapply( 1:dim(T3_x_ref)[1], function(i){
  sqrt( (T3_x_ref[i,2] - T3y_centroid[1])^2 + (T3_x_ref[i,3] - T3y_centroid[2])^2 )
}))


wilcox.test(T3_dist_x , T3_dist_y, alternative = "less")
p2 <- ggplot(umap_coords, aes(x=V1, y=V2,  color=Molecular_clusters )) +  geom_point(size=1, alpha =.5)+scale_color_viridis(discrete=TRUE) + geom_point(aes(x =T3y_centroid[1], y = T3y_centroid[2] )) +
  geom_jitter(size=1.3, alpha=0.6, position = position_jitter(width = 1, height = 1))+   geom_point(size=1, alpha =.5)+ geom_point(aes(x =T3x_centroid[1], y = T3x_centroid[2]))

p2
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Test 4
# Est ce que les LCNEC TypeI sont plus proche des LCNEC/TypeII ou des SCLC/SCLC?
summary(as.factor(Molecular_clusters))
T4x_umap1_df <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "LCNEC/TypeII"),]
T4x_centroid <- apply(T4x_umap1_df[,2:3], 2,mean)

T4_x_ref <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "LCNEC/TypeI"),]

T4_dist_x <- unlist(lapply( 1:dim(T4_x_ref)[1], function(i){
  sqrt( (T4_x_ref[i,2] - T4x_centroid[1])^2 + (T4_x_ref[i,3] - T4x_centroid[2])^2 )
}))



T4y_umap1_df <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "SCLC/SCLC-like"),]
T4y_centroid <- apply(T4y_umap1_df[,2:3], 2,mean)

T4_x_ref <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "LCNEC/TypeI"),]

T4_dist_y <- unlist(lapply( 1:dim(T4_x_ref)[1], function(i){
  sqrt( (T4_x_ref[i,2] - T4y_centroid[1])^2 + (T4_x_ref[i,3] - T4y_centroid[2])^2 )
}))


wilcox.test(T4_dist_x , T4_dist_y, 'less')

p2 <- ggplot(umap_coords, aes(x=V1, y=V2,  color=Molecular_clusters )) +  geom_point(size=1, alpha =.5)+scale_color_viridis(discrete=TRUE) + geom_point(aes(x =T4y_centroid[1], y = T4y_centroid[2] )) +
  geom_jitter(size=1.3, alpha=0.6, position = position_jitter(width = 1, height = 1))+   geom_point(size=1, alpha =.5) + geom_point(aes(x =T4x_centroid[1], y = T4x_centroid[2]))

p2


# %%%%%%%%%%%%%%%%%%%%%%%%%%

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%n%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Test 5
# Est ce que les supra carcinoids clusters + avec les G1 ou avec les G2?
summary(as.factor(Molecular_clusters))
T6x_umap1_df <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "LCNEC/TypeII"| umap_coords_type$Molecular_clusters== "LCNEC/TypeI"|umap_coords_type$Molecular_clusters== 'SCLC/LCNEC-like' ),]
T6x_centroid <- apply(T6x_umap1_df[,2:3], 2,mean)

T6_x_ref <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "Supra_carcinoid"),]

T6_dist_x <- unlist(lapply( 1:dim(T6_x_ref)[1], function(i){
  sqrt( (T6_x_ref[i,2] - T6x_centroid[1])^2 + (T6_x_ref[i,3] - T6x_centroid[2])^2 )
}))



T6y_umap1_df <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "SCLC/SCLC-like"| umap_coords_type$Molecular_clusters== "LCNEC/SCLC-like"),]
T6y_centroid <- apply(T6y_umap1_df[,2:3], 2,mean)

T6_x_ref <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "Supra_carcinoid"),]

T6_dist_y <- unlist(lapply( 1:dim(T6_x_ref)[1], function(i){
  sqrt( (T6_x_ref[i,2] - T6y_centroid[1])^2 + (T6_x_ref[i,3] - T6y_centroid[2])^2 )
}))


wilcox.test(T6_dist_x , T6_dist_y, 'less')

p2 <- ggplot(umap_coords, aes(x=V1, y=V2,  color=Molecular_clusters )) +  geom_point(size=1, alpha =.5)+scale_color_viridis(discrete=TRUE) + geom_point(aes(x =T6y_centroid[1], y = T6y_centroid[2] )) +
  geom_jitter(size=1.3, alpha=0.6, position = position_jitter(width = 1, height = 1))+   geom_point(size=1, alpha =.5)+scale_color_viridis(discrete=TRUE) + geom_point(aes(x =T6x_centroid[1], y = T6x_centroid[2]))

p2


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Test 5
# Est ce que les Supra-Carcinoid sont plus proche des LCNEC/TypeII ou des SCLC/SCLC?
summary(as.factor(Molecular_clusters))
T6x_umap1_df <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "LCNEC/TypeII"),]
T6x_centroid <- apply(T6x_umap1_df[,2:3], 2,mean)

T6_x_ref <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "Supra_carcinoid"),]

T6_dist_x <- unlist(lapply( 1:dim(T6_x_ref)[1], function(i){
  sqrt( (T6_x_ref[i,2] - T6x_centroid[1])^2 + (T6_x_ref[i,3] - T6x_centroid[2])^2 )
}))



T6y_umap1_df <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "LCNEC/TypeI"),]
T6y_centroid <- apply(T6y_umap1_df[,2:3], 2,mean)

T6_x_ref <- umap_coords_type[which(umap_coords_type$Molecular_clusters== "Supra_carcinoid"),]

T6_dist_y <- unlist(lapply( 1:dim(T6_x_ref)[1], function(i){
  sqrt( (T6_x_ref[i,2] - T6y_centroid[1])^2 + (T6_x_ref[i,3] - T6y_centroid[2])^2 )
}))


wilcox.test(T6_dist_x , T6_dist_y, 'less')

p2 <- ggplot(umap_coords, aes(x=V1, y=V2,  color=Molecular_clusters )) +  geom_point(size=1, alpha =.5)+scale_color_viridis(discrete=TRUE) + geom_point(aes(x =T6y_centroid[1], y = T6y_centroid[2] )) +
  geom_jitter(size=1.3, alpha=0.6, position = position_jitter(width = 1, height = 1))+   geom_point(size=1, alpha =.5)+scale_color_viridis(discrete=TRUE) + geom_point(aes(x =T6x_centroid[1], y = T6x_centroid[2]))

p2



