# analyse avec Raphael 22/8/2014

library(rms)
library(gdata)
library(gmodels)
library(chron)
library(psy)
library(poLCA)
#
#-------------------------------------------------------------------------------
#-- Import données
#-------------------------------------------------------------------------------
#dat <- read.csv("dpool3.csv", sep=";", dec=".", na.strings=c(""," ","NA","X"))

# RASS manquant des 72 premiers
#adat <- read.csv("MissedDATA-Pogne1.csv", sep=";", dec=".", na.strings=c(""," ","NA","X"))
#dat$rass[dat$obs %in% adat$obs] <-  adat$rass[adat$obs %in% dat$obs]

getwd()
setwd("/Users/rohaut/Dropbox/POSGNES_4")


#
dat <- read.csv("dpool-def.csv", sep=";", dec=".", na.strings=c(""," ","NA","X"))
dat$delirium <- as.numeric(as.character(dat$Delirium))

#-------------------------------------------------------------------------------
#-- Classes latentes
#-------------------------------------------------------------------------------

dlat <- dat[!is.na(dat$toux)&!is.na(dat$Rphot)&!is.na(dat$Rcorn)&!is.na(dat$Rocg)&
             !is.na(dat$myosis)&!is.na(dat$PMF),]
dlat$Rphot <- as.numeric(dlat$Rphot)+1
dlat$Rcorn <- as.numeric(dlat$Rcorn)+1
dlat$Rocg <- as.numeric(dlat$Rocg)+1
dlat$toux <- as.numeric(dlat$toux)+1
dlat$myosis <- as.numeric(dlat$myosis)+1
dlat$PMF <- as.numeric(dlat$PMF)+1
dlat$Goc <- dlat$Goc+1      # Codé au départ: 0=1, 1=2+
dlat$Gmot <- dlat$Gmot+1    # Codé au départ: 0=1, 1=2+
dlat$Gm2 <- as.numeric(cut(dlat$repmot,c(0,1,3,6)))
dlat$rassc <- as.numeric(dlat$rass > -5)+1

# sans variables avec elts du glasgow binarisés 1/2+
f <- cbind(Rphot,Rcorn,Rocg,toux,myosis,PMF,Goc,Gmot)~1
lcg1 <- poLCA(f,dlat,nclass=1)
lcg2 <- poLCA(f,dlat,nclass=2,nrep=10)
lcg3 <- poLCA(f,dlat,nclass=3,nrep=10)
lcg4 <- poLCA(f,dlat,nclass=4,nrep=10)
lcg5 <- poLCA(f,dlat,nclass=5,maxiter=3000,nrep=10)

crits <- matrix(nrow=5,ncol=2)
crits[1,] <- c(lcg1$aic,lcg1$bic)
crits[2,] <- c(lcg2$aic,lcg2$bic)
crits[3,] <- c(lcg3$aic,lcg3$bic)
crits[4,] <- c(lcg4$aic,lcg4$bic)
crits[5,] <- c(lcg5$aic,lcg5$bic)
plot(crits[,1],type="o",ylim=c(1500,2000))
lines(crits[,2],type="o", col="red")

poLCA(f,dlat,nclass=2, graphs=T)
poLCA(f,dlat,nclass=3,nrep=10, graphs=T)

dlat$predclass3 <- factor(lcg3$predclass) # Retenu sur la base d'AIC (et de Tarek)
dlat$predclass2 <- lcg2$predclass

# elts du glasgow mieux détaillés
f <- cbind(Rphot,Rcorn,Rocg,toux,myosis,PMF,Goc,Gm2)~1
lcgc1 <- poLCA(f,dlat,nclass=1)
lcgc2 <- poLCA(f,dlat,nclass=2,nrep=10)
lcgc3 <- poLCA(f,dlat,nclass=3,nrep=10)
lcgc4 <- poLCA(f,dlat,nclass=4,nrep=10)
lcgc5 <- poLCA(f,dlat,nclass=5,maxiter=3000,nrep=10)

crits <- matrix(nrow=5,ncol=2)
crits[1,] <- c(lcgc1$aic,lcgc1$bic)
crits[2,] <- c(lcgc2$aic,lcgc2$bic)
crits[3,] <- c(lcgc3$aic,lcgc3$bic)
crits[4,] <- c(lcgc4$aic,lcgc4$bic)
crits[5,] <- c(lcgc5$aic,lcgc5$bic)
plot(crits[,1],type="o",ylim=c(1500,2000))
lines(crits[,2],type="o", col="red")

dlat$predgc4 <- lcgc4$predclass
#

# elts du glasgow : que oculaire
f <- cbind(Rphot,Rcorn,Rocg,toux,myosis,PMF,Goc)~1
lcgo1 <- poLCA(f,dlat,nclass=1)
lcgo2 <- poLCA(f,dlat,nclass=2,nrep=10)
lcgo3 <- poLCA(f,dlat,nclass=3,nrep=10)
lcgo4 <- poLCA(f,dlat,nclass=4,nrep=10)
lcgo5 <- poLCA(f,dlat,nclass=5,maxiter=3000,nrep=10)

crits <- matrix(nrow=5,ncol=2)
crits[1,] <- c(lcgo1$aic,lcgo1$bic)
crits[2,] <- c(lcgo2$aic,lcgo2$bic)
crits[3,] <- c(lcgo3$aic,lcgo3$bic)
crits[4,] <- c(lcgo4$aic,lcgo4$bic)
crits[5,] <- c(lcgo5$aic,lcgo5$bic)
plot(crits[,1],type="o",ylim=c(1200,2000))
lines(crits[,2],type="o", col="red")

dlat$predgo4 <- lcgo4$predclass

# sans variables avec rass
f <- cbind(Rphot,Rcorn,Rocg,toux,myosis,PMF,rassc)~1
lcr1 <- poLCA(f,dlat,nclass=1)
lcr2 <- poLCA(f,dlat,nclass=2,nrep=10)
lcr3 <- poLCA(f,dlat,nclass=3,maxiter=3000,nrep=10)
lcr4 <- poLCA(f,dlat,nclass=4,maxiter=3000,nrep=10)
lcr5 <- poLCA(f,dlat,nclass=5,maxiter=3000,nrep=10)

crits <- matrix(nrow=5,ncol=2)
crits[1,] <- c(lcg1$aic,lcg1$bic)
crits[2,] <- c(lcg2$aic,lcg2$bic)
crits[3,] <- c(lcg3$aic,lcg3$bic)
crits[4,] <- c(lcg4$aic,lcg4$bic)
crits[5,] <- c(lcg5$aic,lcg5$bic)
plot(crits[,1],type="o",ylim=c(1500,2000))
lines(crits[,2],type="o", col="red")

poLCA(f,dlat,nclass=3,maxiter=3000,nrep=10, graphs=T)

# sans variables ni elts du glasgow
f <- cbind(Rphot,Rcorn,Rocg,toux,myosis,PMF)~1
lc1 <- poLCA(f,dlat,nclass=1)
lc2 <- poLCA(f,dlat,nclass=2,nrep=10)
lc3 <- poLCA(f,dlat,nclass=3,maxiter=3000,nrep=10)
lc4 <- poLCA(f,dlat,nclass=4,maxiter=10000,nrep=10)
lc5 <- poLCA(f,dlat,nclass=5,maxiter=3000,nrep=10)

crits <- matrix(nrow=5,ncol=2)
crits[1,] <- c(lc1$aic,lc1$bic)
crits[2,] <- c(lc2$aic,lc2$bic)
crits[3,] <- c(lc3$aic,lc3$bic)
crits[4,] <- c(lc4$aic,lc4$bic)
crits[5,] <- c(lc5$aic,lc5$bic)
plot(crits[,1],type="o",ylim=c(1200,1600))
lines(crits[,2],type="o", col="red")

poLCA(f,dlat,nclass=3,maxiter=3000,nrep=10, graphs=T)


dlat$predsg3 <- factor(lc3$predclass)

dlat$predr3 <- factor(lcr3$predclass)


#dlat$predclass3 <- reorder(dlat$predclass3, new.order=c(2,3,1))

lrm(dcJ28 ~ predclass3,data=dlat)
lrm(dcJ28 ~ predsg3,data=dlat)
lrm(dcJ28 ~ predr3,data=dlat)

lrm(dcJ28 ~ factor(toux),data=dlat)

lrm(Coma_.Delirium ~ predclass3,data=dlat)
lrm(Coma_.Delirium ~ predsg3,data=dlat)
lrm(Coma_.Delirium ~ predr3,data=dlat)

lrm(delirium ~ predclass3,data=dlat)
lrm(delirium ~ predsg3,data=dlat)
lrm(delirium ~ predr3,data=dlat)

#source("C:/Users/Raphael Porcher/Documents/Raphael/Scripts/utils.r")
source("/Users/rohaut/Dropbox/POSGNES_3/analyses Aout 2013/utils.r")

temp <- NULL
temp <- rbind(c("No. patients",table(dlat$predclass3),""),temp)
temp <- descr("repoc", y="predclass3", dat=dlat, res=temp, type="c", desc="med", prec=0, test="kruskal")
temp <- descr("repmot", y="predclass3", dat=dlat, res=temp, type="c", desc="med", prec=0, test="kruskal")
temp <- descr("Rphot", y="predclass3", dat=dlat, res=temp, type="b",categ="2", prec=0, test="fisher")   # 2 = présent
temp <- descr("Rcorn", y="predclass3", dat=dlat, res=temp, type="b",categ="2", prec=0, test="fisher")
temp <- descr("Rocg", y="predclass3", dat=dlat, res=temp, type="b",categ="2", prec=0, test="fisher")
temp <- descr("toux", y="predclass3", dat=dlat, res=temp, type="b",categ="2", prec=0, test="fisher")
temp <- descr("myosis", y="predclass3", dat=dlat, res=temp, type="b",categ="2", prec=0, test="fisher")
temp <- descr("PMF", y="predclass3", dat=dlat, res=temp, type="b",categ="2", prec=0, test="fisher")
temp <- descr("hypno", y="predclass3", dat=dlat, res=temp, type="c", desc="msd", prec=1, test="kruskal")
temp <- descr("hypnok", y="predclass3", dat=dlat, res=temp, type="c", desc="msd", prec=1, test="kruskal")
temp <- descr("SUFenta", y="predclass3", dat=dlat, res=temp, type="c", desc="msd", prec=1, test="kruskal")
temp <- descr("SUFentak", y="predclass3", dat=dlat, res=temp, type="c", desc="msd", prec=1, test="kruskal")
temp <- descr("age", y="predclass3", dat=dlat, res=temp, type="c", desc="med", prec=0, test="kruskal")
temp <- descr("sexe..1.f.", y="predclass3", dat=dlat, res=temp, type="b",categ="1", prec=0, test="fisher")     # 1=F, 2=M
temp <- descr("sapsII", y="predclass3", dat=dlat, res=temp, type="c", desc="med", prec=0, test="kruskal")
temp <- descr("rass", y="predclass3", dat=dlat, res=temp, type="c", desc="med", prec=0, test="kruskal") 
temp <- descr("sepsis", y="predclass3", dat=dlat, res=temp, type="b",categ="1", prec=0, test="fisher")
temp <- descr("Coma", y="predclass3", dat=dlat, res=temp, type="b",categ="1", prec=0, test="fisher")
temp <- descr("Delirium", y="predclass3", dat=dlat, res=temp, type="b",categ="1", prec=0, test="fisher")
temp <- descr("Coma_.Delirium", y="predclass3", dat=dlat, res=temp, type="b",categ="1", prec=0, test="fisher")
temp <- descr("dcICU", y="predclass3", dat=dlat, res=temp, type="b",categ="1", prec=0, test="fisher")
temp <- descr("dcJ28", y="predclass3", dat=dlat, res=temp, type="b",categ="1", prec=0, test="fisher") 

disp.epi("predclass3","dcJ28",dlat,NULL)

lrm(dcJ28 ~ predclass3,data=dlat)
lrm(dcJ28 ~ sapsII + predclass3,data=dlat)
lrm(dcJ28 ~ sapsII,data=dlat)

dnat <- dlat[!is.na(dlat$dcJ28)&!is.na(dlat$age)&!is.na(dlat$sapsII),]
lrm(dcJ28 ~ predclass3,data=dnat)
lrm(dcJ28 ~ sapsII + predclass3,data=dnat)
lrm(dcJ28 ~ sapsII,data=dnat)

# Attention, le RASS est très différent entre les classes latentes!!!!
lrm(dcJ28 ~ sapsII + rass +  predclass3,data=dnat)
lrm(dcJ28 ~ sapsII + rass,data=dnat)

# poLCA.posterior pour les probas postérieures d'appartenance aux classes
odat <- read.csv("dpool3.csv", sep=";", dec=".", na.strings=c(""," ","NA","X"))[1:72,]
odat <- odat[!is.na(odat$toux)&!is.na(odat$Rphot)&!is.na(odat$Rcorn)&!is.na(odat$Rocg)&
             !is.na(odat$myosis)&!is.na(odat$PMF),]
odat$Rphot <- as.numeric(odat$Rphot)+1
odat$Rcorn <- as.numeric(odat$Rcorn)+1
odat$Rocg <- as.numeric(odat$Rocg)+1
odat$toux <- as.numeric(odat$toux)+1
odat$myosis <- as.numeric(odat$myosis)+1
odat$PMF <- as.numeric(odat$PMF)+1
odat$Goc <- odat$Goc+1      # Codé au départ: 0=1, 1=2+
odat$Gmot <- odat$Gmot+1    # Codé au départ: 0=1, 1=2+
odat$Gm2 <- as.numeric(cut(odat$repmot,c(0,1,3,6)))

preds <- poLCA.posterior(lcg3,y=odat[,c("Rphot","Rcorn","Rocg","toux","myosis","PMF","Goc","Gmot")])
odat$predclass3 <- factor(apply(preds,1,function(x){which(x==max(x))}))


lrm(dcJ28 ~ predclass3,data=odat)
lrm(dcJ28 ~ sapsII + predclass3,data=odat)

disp.epi("predclass3","dcJ28",dlat,NULL)
disp.epi("predclass3","dcJ28",odat,NULL)

tempo <- c("No. patients",table(odat$predclass3),"")
tempo <- descr("dcJ28", y="predclass3", dat=odat, res=tempo, type="b",categ="1", prec=0, test="fisher")     # 1=F, 2=M
 binom.test(1,19,0.13)
 binom.test(10,19,0.56)
 binom.test(9,32,0.27)