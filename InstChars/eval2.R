library(plyr)
library(xtable)

d=read.table("all.txt",h=T)

d$Fullnam=ifelse(is.na(d$no),
          ifelse(is.na(d$group),paste(d$Name,d$base,sep=""),paste(d$Name,d$base,"_",d$group,sep="")),
          ifelse(is.na(d$group),paste(d$Name,d$base,"-",d$no,sep=""),paste(d$Name,d$base,"_",d$no,"_",d$group,sep="")))

## (1) compile table, w/ MSS/Liu instances separately
p1=ddply(subset(d,Source!="MSSliu"),.(Source,N,M),summarize,n=length(M),OS=mean(OS))
## (1.1) for MSS/Liu: show base number of jobs and number of open shop jobs
p2=ddply(subset(d,Source=="MSSliu"),.(Source,N-no,no,M),summarize,n=length(M),OS=mean(OS))
p2[,2]=paste(p2[,2],p2[,3],sep="+")
p2$no=NULL

## (2) join the tables
colnames(p1)=c("Ref.","n","m","#","OS")
colnames(p2)=c("Ref.","n","m","#","OS")
a=rbind(p1,p2)

## (2.1) normalize the names
## TBD
normname = function(name) {
    switch(as.numeric(name),
        "\\textcite{NasiriGA}", # GSSNasAdmu
        "\\textcite{NasiriGA}", # GSSNasTai
        "\\textcite{blum2004ant}", # GSSblum
        "\\textcite{Liu2004metaMixed}", # MSSliu
        "\\textcite{NasiriHyb}", # NasiriPSS
        "\\textcite{Taillard1993bench}", # OSSTai
        "\\textcite{brucker1997branch}", # OSSBrucker
        "\\textcite{gueret1999new}" # OSSgueret
##        "", # GSSNasAdmu
##        "\\textcite{NasiriGA}", # GSSNasTai
##        "\\textcite{blum2004ant}", # GSSblum
##        "\\textcite{Liu2004metaMixed}", # MSSliu
##        "PSSata", # NasiriPSS
##        "\\textcite{Taillard1993bench}", # OSSTai
##        "\\textcite{brucker1997branch}", # OSSBrucker
##        "\\textcite{gueret1999new}", # OSSgueret
        )
}
a$Ref.=sapply(a$Ref.,FUN=normname)

## (3) print LaTeX table in two-column format
a2=cbind(a[1:31,],a[32:62,])
#print(xtable(a2,digits=c(0,0,rep(c(1,1,0),4)),caption=table1.caption,label="tab:exp1"),file="./table1.tex",hline.after=c(-1,0,nrow(table1)-1,nrow(table1)),size="\\small\n\\medskip\\setlength{\\tabcolsep}{1ex}",floating.environment="table*",add.to.row=table1.addtorow,format.args=list(big.mark=","))
## TBD
