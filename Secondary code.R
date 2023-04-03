#Data in file"Secondary data.Rdata" is used to run the code in file"Secondary code.R" 
#Loading data
load("./Secondary data.Rdata")


#Clustering of time series protein expression (Mfuzz) for HCC-related differentially abundant proteins
#The quantitative profiles in file "Mfuzz_data" was used in the part of the code and it obtained by log2 conversion of the average quantification of 34 HCC-related differentially abundant proteins.
#install.packages("Mfuzz")
library(Mfuzz)
rownames(Mfuzz_data)<-Mfuzz_data[,1]
Mfuzz_data<-Mfuzz_data[,-1]
Mfuzz_data<-as.matrix(Mfuzz_data)
#Building objects
Mfuzz_data <- new('ExpressionSet',exprs = Mfuzz_data)
#Processing NA value
Mfuzz_data <- filter.NA(Mfuzz_data, thres = 0.25)
Mfuzz_data <- fill.NA(Mfuzz_data, mode = 'mean')
#Remove proteins with small differences between samples according to SD
Mfuzz_data <- filter.std(Mfuzz_data, min.std = 0)
#Standardization
Mfuzz_data <- standardise(Mfuzz_data)
#Evaluating the best m value to prevent random data clustering
m <- mestimate(Mfuzz_data)
m
#Manually defining the number of clusters
n <- 4
#Mfuzz clustering and drawing
cl <- mfuzz(Mfuzz_data, c = n, m = m)
p<-mfuzz.plot(Mfuzz_data, cl = cl, mfrow = c(2, 2), time.labels = c("AsC", "BLD", "LC","HCC"))
#Integrating and outputing the list of protein clusters
cl$size
head(cl$cluster)
cl$membership
protein_cluster <- cbind(cl$cluster, cl$membership)
colnames(protein_cluster)[1] <- 'cluster'
write.table(protein_cluster, "Differential abundance protein_cluster4.txt", sep = '\t', col.names = NA, quote = FALSE)


#The pearson’s correlation coefficient distribution of MS data of Hela standards
#install.packages("reshape2")
#install.packages("corrplot")
library(reshape2)
library(corrplot)
rownames(Hela_DDA)<-Hela_DDA[,1]
Hela_DDA<-Hela_DDA[,-1]
a1<-cor(Hela_DDA,method = "pearson")
corrplot.mixed(a1,upper = "ellipse",number.cex=0.4,tl.cex=0.2)
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
col1 <- colorRampPalette(c("#7F0000", "red", "#FF7F00", "yellow", "white",
                           "cyan", "#007FFF", "blue","#00007F"))
color_1<-colorRampPalette(c("cyan","magenta"))
color_2<-colorRampPalette(c("magenta","cyan"))
palette_1<-RColorBrewer::brewer.pal(n=11,name="RdYlGn")
palette_2<-rev(palette_1)
tiff(height=1800, width=2100,pointsize = 2,res=500, file="Hela_DDA_corrplot.tiff")
corrplot(a1, method = "ellipse", col = palette_2,
         addCoefasPercent = TRUE,
         diag = TRUE, mar = c(1,1,1,1),type = "lower",tl.pos="d")
corrplot(a1, method = "color",type = "upper",add = T ,diag = FALSE,tl.pos="n",
         addCoef.col = "black", 
         addCoefasPercent = F)
dev.off()
#The image of the correlation coefficient matrix is saved in the current working directory.

#Calculation of grouping pearson’s correlation
#install.packages("reshape2")
#install.packages("tidyverse")
library(reshape2)
library(tidyverse)
dataTN<-read.csv("DIA定量原始分组数据.csv")
rownames(dataTN)<-dataTN[,1]
dataTN<-dataTN[,-1]
Cor_AsC <- cor(dataTN[1:40]) %>%
  melt()%>%
  (function(x){
    x <- x[which(x[,1]!=x[,2]),]
  })
Cor_BLD <- cor(dataTN[41:104]) %>%
  melt()%>%
  (function(x){
    x <- x[which(x[,1]!=x[,2]),]
  })
Cor_LC <- cor(dataTN[105:157]) %>%
  melt()%>%
  (function(x){
    x <- x[which(x[,1]!=x[,2]),]
  })
Cor_HCC <- cor(dataTN[158:320]) %>%
  melt()%>%
  (function(x){
    x <- x[which(x[,1]!=x[,2]),]
  })
write.csv(Cor_AsC,"Cor_AsC.csv")
write.csv(Cor_BLD,"Cor_BLD.csv")
write.csv(Cor_LC,"Cor_LC.csv")
write.csv(Cor_HCC,"Cor_HCC.csv")
#The above saved files can get the correlation coefficient file"cor_data" by deleting duplicates and integrating them.

#Ridge diagram of correlation coefficient
#install.packages("ggridges")
#install.packages("viridis")
library(ggridges)
library(viridis)
cor_data<-read.csv("cor.data.4.csv")
ggplot(cor_data, aes(x=`exp`, y=`group`, fill=..x..))+
  geom_density_ridges_gradient(scale=2, rel_min_height=0.01, gradient_lwd = 2.)+
  scale_x_continuous(expand = c(0.01, 0),limits = c(0.85,1))+ # 扩展下横轴和纵轴
  scale_y_discrete(expand = c(0.01,0))+
  scale_fill_viridis(name="R^2", option = "C")+
  theme_ridges(font_size = 20, grid = FALSE)+
  theme(axis.title.y = element_blank())
# Names of group a-j correspond to the row named "group.name" in cor_data. 