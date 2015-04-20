data<-read.table("labels.txt", header=FALSE, sep="\t")
write.table(t(data), file="labels2.txt", row.names=FALSE, col.names=FALSE, quote=FALSE)

