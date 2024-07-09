tree = data.frame(culture], stringsAsFactors = F)
# clarify taxonomy
tree[,1] = paste("p__",tree[,1],sep = "")
tree[,2] = paste("c__",tree[,2],sep = "")
tree[,3] = paste("o__",tree[,3],sep = "")
#tree[,4] = paste("f__",tree[,4],sep = "")
tree[,5] = paste("g__",tree[,5],sep = "")
tree[idx,4] = paste(tree[idx,4], 1:length(tree[idx,4]))
#idx = tree[,4] %in% "Unassigned"
#tree = tree[!idx,]
tree[,4] = gsub('_\\w*',"",tree[,4])
# save tree backbone
write.table (tree, file="tree1_backbone.txt", sep=".", col.names=F, row.names=F, quote=F)

Phylum = unique(tree[,1]) 
Class = unique(tree[,2])
Order = unique(tree[,3])
Family = unique(tree[,4])
Genus = unique(tree[,5])

pro = tree[tree[,1]=="p__Proteobacteria",4]
act = tree[tree[,1]=="p__Actinobacteriota",4] 
bac = tree[tree[,1]=="p__Bacteroidota",4]
fir = tree[tree[,1]=="p__Firmicutes",4]


label_color = data.frame(stringsAsFactors = F)
for (element in Family)
{
  # element
  anno = data.frame(stringsAsFactors = F)
  anno[1,1] = element
  anno[1,2] = "annotation"
  anno[1,3] = "*"
  anno[2,1] = element
  anno[2,2] = "annotation_rotation"
  anno[2,3] = "90"
  anno[3,1] = element
  anno[3,2] = "annotation_background_color" 
  
  if (element %in% pro)
  {
    anno[3,3] = "#85F29B"
  } else if (element %in% act)
  {
    anno[3,3] = "#F58D8D"   
  } else if (element %in% fir)
  {
    anno[3,3] = "#F7C875"  
  } else if (element %in% bac)
  {
    anno[3,3] = "#91DBF6"   
  } else {
    anno[3,3] = "grey"   
  }
  label_color = rbind(label_color,anno)
}
write.table(label_color, "tree2_label_color.txt", sep = "\t", quote = F,col.names = F,row.names = F, na="")

conda activate culturome
# mamba install graphlan
rm -rf track*
  cat cfg/global.cfg tree2_label_color.txt > track0
cat track* > graphlan_annotate.txt
graphlan_annotate.py --annot graphlan_annotate.txt tree1_backbone.txt graphlan.xml
graphlan.py graphlan.xml graphlan0_tree.pdf --size 4.5 --dpi 300 --external_legends
