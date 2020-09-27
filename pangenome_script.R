library(micropan)
library(ggplot2)
library(ggdendro)
library(UpSetR)
library(reshape2)
library(tidyverse)
# dependencias 
#archivo micropan-source.R
#debe de correrse antes de correr este script 
source('~/Escritorio/tutorial_pangenome/micropan-source.R')
x<-read.delim("~/Escritorio/tutorial_pangenome/Spneumonie.proteinortho.tsv", header = T, sep = '\t') # aqui va el archivo de salida que da proteinortho5

##### make pangenome analysis ####
panmat <- as.panmat(x) #hace la matrix de pangenoma
summary(panmat) # hace un resumen de cuantos clusters encontro por cada genoma (basic stats)
rarefy <- rarefaction(panmat, n.perm=100) #hacer la rarefraction curve 
summary(rarefy) #summary of the rarefaction results
jpeg("rarefaction.jpg") #refraction plot
plot(rarefy, ylim=c(0, max(rarefy)))
dev.off()
#Once we have the pan-matrix, we can make a bar-chart over how many clusters are found in 1,2,...,all genomes:
tibble(Clusters = as.integer(table(factor(colSums(panmat > 0),
                                          levels = 1:nrow(panmat)))), 
       Genomes = 1:nrow(panmat)) %>% 
  ggplot(aes(x = Genomes, y = Clusters)) +
  geom_col() 

heap <- heaps(panmat, 1000) #pan-genome openness  more than 1.0 is close
chao(panmat) #how bigger is the pan-genome 
fitted <- binomixEstimate(panmat)
print(fitted$BIC.table)
#compute distances between genomes and dendrogram tree to visualize their relations
d.man <- dist(panmat,method ="manhattan")
ggdendrogram(dendro_data(hclust(d.man,method = "average")),
             rotate = TRUE, theme_dendro = FALSE) +
labs(x = "Genomes", y = "Euclidean distance", title = "Pan-genome dendrogram")

#### make binary matrix and upset ####
y<-as.matrix(t(panmat))
y[y>1] = 1
m<-as.data.frame(y)
set= list(colnames(m))

upset(m,sets = colnames(m),
      nintersects =30,
      nsets = 6,
      keep.order = T,
      set_size.show = F,
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.4, 0.6),
      number.angles = 0, 
      text.scale = 1.2, 
      point.size = 2.8, 
      line.size = 1,
      mainbar.y.label = "Intersections", 
      sets.x.label = "Gene Clusters",
      queries = list(list(query = intersects, params = set, active=T, color="#b2182b")
      )
)

##### make tables and save ####
m$count<-rowSums(m)
core<-as.numeric(length(rownames(subset(m,count==10))))
acc<-as.numeric(length(rownames(subset(m,count >1 &count <10 ))))
uniq<-subset(m,count ==1)
mu<-melt(uniq)
totals<-mu%>% 
  group_by(variable) %>% 
  mutate(id = 1:n()) %>% 
  summarise(value = sum(value))
totals <- totals[-c(11), ]
totals<-totals%>% add_row(variable="Accessory", 
                          value= acc[1])
totals<-totals%>% add_row(variable="Core", 
                          value= core[1])

colnames(totals)<- c("strain", "count")
write.table(totals, 
            file = "~/Escritorio/tutorial_pangenome/counts.list",
            sep = "\t",
            row.names = F)

core_subset<-subset(x,X..Species==10)
core_subset<-as.data.frame(core_subset$GCF_000007045.1)
write.table(core_subset, 
            file = "~/Escritorio/tutorial_pangenome/core.list",
            sep = "\t",
            row.names = F,
            col.names = F,
            quote = F)

write.table(y,"~/Escritorio/tutorial_pangenome/binary_matrix.tab",
            sep = '\t',
            row.names = T, 
            quote = F)

write.table(y,"~/Escritorio/tutorial_pangenome/binary_panGP.tab",
            sep = '',
            row.names = F,
            col.names = F,
            quote = F)
