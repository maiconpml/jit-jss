rm(list = ls()) # clean workspace


library(plyr)
library(ggplot2)
library(stringr)
library(xtable)
library(reshape2)

rd = function(v,o) { 100*(v-o)/o }

# basic settings for producing LaTeX tables
options(width=240)
options(xtable.caption.placement="top")
options(xtable.booktabs=T)
options(xtable.sanitize.text.function=function(x){x})
options(xtable.table.placement="t")
options(xtable.include.rownames=F)

##################################################################
# READ DATA
##################################################################
psp = read.table("Data/pspChars.dat",h=T)
msp = read.table("Data/mspChars.dat",h=T)
gspBlum = read.table("Data/gspBlumChars.dat",h=T)
gspNasAdmu = read.table("Data/gspNasAdmuChars.dat",h=T)
gspNasTai = read.table("Data/gspNasTaiChars.dat",h=T)
ospBrucker = read.table("Data/ospBruckerChars.dat",h=T)
ospGueret = read.table("Data/ospGueretChars.dat",h=T)
ospTai = read.table("Data/ospTaiChars.dat",h=T)

##################################################################
# TABLES
##################################################################
#PSP
psp$num = 1
pspT = ddply(psp, .(M, N), summarize, OS=mean(OS), num=sum(num))
#pspT$O = pspT$M * pspT$N

#MSP
msp$num=1
mspT = ddply(msp, .(M, N), summarize, OS=mean(OS), num=sum(num))
#mspT$O = mspT$M * mspT$N

#GSP
gspBlum$num = 1
gspBlumT = ddply(gspBlum, .(M, N), summarize, OS=mean(OS), num=sum(num))

gspNasTai$num = 1 
gspNasTaiT = ddply(gspNasTai, .(M, N), summarize, OS=mean(OS), num=sum(num))

gspNasAdmu$num = 1 
gspNasAdmuT = ddply(gspNasAdmu, .(M, N), summarize, OS=mean(OS), num=sum(num))

#OSP
ospBrucker$num = 1
ospBruckerT = ddply(ospBrucker, .(M, N), summarize, OS=mean(OS), num=sum(num))

ospGueret$num = 1
ospGueretT = ddply(ospGueret, .(M, N), summarize, OS=mean(OS), num=sum(num))

ospTai$num = 1
ospTaiT = ddply(ospTai, .(M, N), summarize, OS=mean(OS), num=sum(num))

