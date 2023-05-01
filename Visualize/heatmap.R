library(pheatmap)
library(grid)

df<-read.table("celltype_communication.txt",header=T,row.names = 1)
df1<-data.matrix(df)
max1<-round(max(df1),2)
min1<-round(min(df1),2)
mean1<-round(mean(df1),2)
step1 = 0.001
bk <- c(seq(min1,mean1,by=step1),seq(mean1+step1,max1,by=step1))
bk1 = seq(min1,mean1,by=step1)
bk2 = seq(mean1+step1,max1,by=step1)
pheatmap(df,
         scale = "none",
         color = c(colorRampPalette(colors = c("Blue","white"))(length(bk1)),colorRampPalette(colors = c("white","red"))(length(bk2))),
         legend = TRUE,
         legend_breaks=seq(min1,max1,mean1-min1),   
         legend_labels=c(min1,mean1,max1),
         breaks=bk,
         cluster_rows=FALSE,
         cluster_cols=FALSE,
         width = 7,height = 6.5,    # users can choose another size of the image to be saved
         filename = 'plot.pdf')   

