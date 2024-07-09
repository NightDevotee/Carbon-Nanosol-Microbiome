library(ggplot2)
library(reshape2)
library(vegan)
library(car)
library(ggalluvial)

rm(list=ls())

#Alpha diversity
metadata = read.table("metadata.txt", header=T, sep="\t")
rownames(design)<-design$Sample_Name

otu <- otu[,which(colnames(otu)%in%rownames(design))]
summary(colSums(otu))
otu<-otu[,colSums(otu) > 8500]

sums<-colSums(otu)
otu<-t(rrarefy(t(otu),min(sums)))
otu<-as.data.frame(otu)

data <-t(otu)
alpha <- alpha_diversity(data)
alpha <-apply(alpha,2,as.numeric)
rownames(alpha)<- rownames(otu)
alpha <-merge(alpha,design,by="row.names")
rownames(alpha)<- alpha[,1]
alpha<- alpha[,-1]
alpha<-tidyr::unite(alpha,"Group",Site,Treatment,sep="",remove=FALSE)

#画图
title="Bacteria"
#生态位
variance<-aov(Chao1 ~ Site, data=alpha) 
LSD <- LSD.test(variance,"Site", p.adj="none")
LSD<-merge(as.data.frame(LSD$means) ,as.data.frame(LSD$groups),by="row.names")
LSD$Group <- LSD$Row.names
p1<- ggplot(alpha, aes(x=Site, y=Chao1,fill= NULL,color= Site)) +
  geom_boxplot(width = 0.6,position=position_dodge(0.8), alpha=0.7, linewidth=1,outlier.shape = NA)+
  geom_jitter(shape=16,size=1,position=position_jitter(0.2))+
  stat_boxplot(geom="errorbar",width=0.16,size=1,position=position_dodge(0.7))+
  geom_text(data=LSD,aes(x=Group,y=Max*1.05,label=groups),color="black",fontface="bold")+
  scale_color_manual(values = c('#91C953', '#94B2D7', '#2FACC9', '#356195'))
  labs(x="",y= "Chao1 index",sep="",title=title)+
  theme_bw()+ 
  theme(legend.position = 'none') 


bray_curtis = read.table("distance-matrix.tsv", sep="\t", header=T, row.names= 1, check.names=F)
pcoa = cmdscale(bray_curtis, k=3, eig=T) # k is dimension, 3 is recommended; eig is eigenvalues
points = as.data.frame(pcoa$points) # get coordinate string, format to dataframme
colnames(points) = c("x", "y", "z") 
eig = pcoa$eig
points = cbind(points, design[match(rownames(points), rownames(design)), ])
points <- points[points$Date!= "T3",]
points <-na.omit(points)
points$Site <-factor(points$Site,level=c('R', 'RP', 'RS', 'BS'))

#PCoA
ggplot(points, aes(x=x, y=y, fill=Site,shape=Site)) +
  geom_point(alpha=.9, size=5) + 
  labs(x= "",
       y=paste("PCoA 2 (", format(100 * eig[2] / sum(eig), digits=4), "%)", sep=""),
       title="Bacteria")+
  scale_shape_manual(values = c(21, 22, 24, 23))+ 
  scale_fill_manual(values =c('#91C953', '#94B2D7', '#2FACC9', '#356195'))+
  theme(
    panel.grid = element_blank(),
    legend.position = c(0.95,0.85),
    legend.title = element_blank(),
    legend.text = element_text(color = 'black',size = 12,face = 'plain')
  )+
  theme(panel.grid.minor = element_blank())+
  theme_classic()+
  theme(axis.text.x = element_blank())+
  theme(axis.text.y=element_text(vjust=1,size=13))+
  theme(axis.title.x=element_text(vjust=2, size=14))+
  theme(axis.title.y=element_text(vjust=2, size=14))




#Composition at phylum level
sample = "asv_table.p.absolute.txt"
data=read.table(paste("./Absolute/",sample,sep=""), sep="\t", header=T, row.names= 1, check.names=F)

otu<-as.data.frame(lapply(data,function(x) as.numeric(as.character(x))))
rownames(otu)<-rownames(data)

sums<-colSums(otu)
otu<-t(rrarefy(t(otu),min(sums)))
otu<-as.data.frame(otu)

n<-length(colnames(otu))
n<-as.numeric(n)

tax <- apply(otu[,1:n], 2, function(x){(x/sum(x))*100})
a<-length(rownames(tax[,-1]))
a<-as.numeric(a)

tax<-tax[order(rowSums(tax),decreasing=T),]
se<-rbind(colSums(tax[c(1,17:a),]),tax[16:2,])
rownames(se)[1]<-"Others"

b<-cbind(t(se), sub_design)
tax<-aggregate(b[,1:16],by=list(sub_design$Group),FUN=mean)
rownames(tax)<-tax[,1]
tax<-t(tax[,2:17])

tax_top15=melt(tax ,value.name ="Relative abundance",variable.name = "Group")
colnames(tax_top15)<-c("Genus","Group","Relative abundance")

tax_top15$`Relative abundance`<-as.numeric(as.character(tax_top15$`Relative abundance`))


ggplot(data=tax_top30,aes(x=Group,y=`Relative abundance`,fill=Genus))+
  geom_bar(position = "stack",stat = 'identity')+
  labs(x=NULL, y="Relative abundance (%)",title=level)+
  scale_fill_manual(values= Color3)+
  theme(legend.position = "right",
        panel.grid =element_blank())+
  guides(fill=guide_legend(title="",color="black",reverse=TRUE))+
  theme_classic()+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(axis.text=element_text(colour='black',size=9))+
  theme(axis.text.x=element_text(size=10,angle = 45, hjust = 1, vjust = 1))
#Save
ggsave(file=paste(level,"pdf", sep="."), p, width = 7, height = 6)



