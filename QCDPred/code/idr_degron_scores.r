
# Load IDR's
idr = read.csv("idr.csv", sep=";")

# write.table(idr[,c("name","sequence")], quote=F, row.names=F, col.names=F, file="idr.seq")
# cmd: ./QCDpred.py idr.seq > idr_qcdpred.txt

qcdpred = read.table("idr_qcdpred.txt")
colnames(qcdpred) = c("name","seq","score","aa","resi")

# assume they are ordered!
agg = aggregate(qcdpred$score, by=list(qcdpred$name), paste0, collapse=":")
idr$QCDpred = agg[match(idr$name,agg[,1]),2]

idr$QCDavg = sapply(idr$QCDpred, function(sl){ val=lapply(strsplit(sl,":"), as.numeric); sapply(val,mean) })
idr$QCDmed = sapply(idr$QCDpred, function(sl){ val=lapply(strsplit(sl,":"), as.numeric); sapply(val,median) })
idr$QCDmax = sapply(idr$QCDpred, function(sl){ val=lapply(strsplit(sl,":"), as.numeric); sapply(val,max) })

write.csv(idr[,c("uniprot","nres_unip","nres_seg","first","last","QCDavg","QCDmax","QCDmed","sequence","QCDpred")], row.names=F, file="idr_qcdpred.csv")

