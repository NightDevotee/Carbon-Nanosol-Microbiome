alpha_diversity <- function(x, tree = NULL) {
  observed_species <- estimateR(x)[1, ]
  Chao1 <- estimateR(x)[2, ]
  ACE <- estimateR(x)[4, ]
  Shannon <- diversity(x, index = 'shannon',base = 2)
  Simpson <- diversity(x, index = 'simpson')    #Gini-Simpson
  goods_Coverage <- 1 - rowSums(x == 1) / rowSums(x)
  
  Shannon <- sprintf("%0.4f", Shannon)
  Simpson <- sprintf("%0.4f", Simpson)
  goods_Coverage <- sprintf("%0.4f", goods_Coverage)
  result <- data.frame(observed_species, ACE,Chao1, Shannon, Simpson, goods_Coverage)
  
  if (!is.null(tree)) {
    PD_whole_tree <- pd(x, tree, include.root = FALSE)[1]
    names(PD_whole_tree) <- 'PD_whole_tree'
    result <- cbind(result, PD_whole_tree)
    result <- data.frame(observed_species, ACE,Chao1, Shannon, Simpson,
                         PD_whole_tree ,goods_Coverage)
  }
  result
}
forceMatrixToInteger <- function(m){
  apply (m, c (1, 2), function (x) {
    (as.integer(x))
  })
}
data <- t(KO)
data<- forceMatrixToInteger(data)

#alpha diversity
alpha <- alpha_diversity (data)
alpha <-merge(alpha,metadata,by="row.names")
alpha$Site <- factor(alpha$Site,levels=c("BS","RS","RP","R"))
title="KO diversity"

p1<-ggplot(alpha, aes(x=Site, y=Chao1,fill=Treatment,color=Treatment)) +
  geom_boxplot(width = 0.6,position=position_dodge(0.8), alpha=0.7, size=1)+
  scale_color_manual(values =color)+
  scale_fill_manual(values =color)+
  geom_jitter(shape=16,size=2,position=position_jitter(0.2))+
  labs(x="",y= "Chao1 index",sep="",title=title)+
  theme_light()+ 
  theme(legend.position = c('none'))


metadata = read.table("metadata.txt", header=T, sep="\t",row.names=1)
KO<-read.table("KO.txt",header=T, sep="\t",row.names=1)
data<-t(KO)
#Bray-curtis
bray_curtis <- as.matrix(vegan::vegdist(data, method = 'bray'))

#NMDS
df_nmds <- metaMDS(bray_curtis, k = 2)
df_nmds$stress
#KO:0.06960607
stressplot(df_nmds)
points <- as.data.frame(df_nmds$points)
points$samples <- row.names(points)
names(points)[1:2] <- c('x', 'y')
points = cbind(points, metadata[match(rownames(points), rownames(metadata)), ])
points$Site <-factor(points$Site,level=c('R', 'RP', 'RS', 'BS'))

color <-c("#E2C054","#CE96E2","#842EBC")
title="NMDS of KEGG Orthology (KO)"
p<-ggplot(points, aes(x=x, y=y, fill=Treatment,shape=Site)) +
  geom_point(alpha=.9, size=4) + 
  labs(x="NMDS1",y="NMDS2", sep="",title=title)+
  scale_shape_manual(values = c(21, 22, 24, 25))+
  scale_fill_manual(values =color)+ 
  theme_bw()+
  theme(panel.grid.minor = element_blank())

