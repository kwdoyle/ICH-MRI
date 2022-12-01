library(openxlsx)
library(reshape2)
library(ggplot2)
library(dplyr)
library(abind)
library(corrplot)

### Analysis 2


setwd("E:/Experiments/BrainFlowMap Stuff/Results_JanClaassen_ICH_DTI")



#### Parameters to set ####
## set data types to look at
datatypes <- c("FA", "MD", "AD", "RD", "fibercount")

## uncomment whichever time of interest to analyze:
#conscious_time <- "unconscious.on.discharge"
conscious_time <- "unconsious.at.time.of.MRI"

## Set to 1 if want to check this for only the patients unconscious at time of MRI
checkUnconscMRI <- 1

## Set to 1 if want to actually save the plots
savePlots <- 0






#### Preprocessing steps ####
txtfiles <- list.files(pattern=".txt")
analysis2 <- grep("Analysis2", txtfiles, value=T)
conscdat <- read.xlsx("../MRI DTI set.xlsx", 1)


# extract "columns" from the txtfiles list; get unique patient IDs from files
split_txt <- strsplit(analysis2, "_")
ID <- sapply(split_txt, "[", 1)
ID.u <- unique(ID)

### behold; the messiest way to replace all instances of "DTI" with "ICH":
# in order to match the names in the text files.
# this does three replacements and uses the output of the prior replacement
# as the input for the next one. it first removes everything after any instance of an "A"
# (because some of the names are like DTI 8A/8B, etc.)
# and then it replaces the DTIs with a space between the number,
# and then replaces the DTIs without a space.
# THEN it removes any blank spaces, which would cause an identical match from the mainTable IDs to appear different
conscdat$DTI.Number <- gsub(" ", "", gsub("DTI", "ICH", gsub("DTI ", "ICH", gsub("A.*", "", conscdat$DTI.Number))))

# convert conscious columns to factor
conscdat[,c("unconsious.on.admission", "unconsious.at.time.of.MRI", "unconscious.on.discharge", "Recovery.of.consciousness")] <- lapply(
  conscdat[,c("unconsious.on.admission", "unconsious.at.time.of.MRI", "unconscious.on.discharge", "Recovery.of.consciousness")],
  as.factor
)






##### make separate arrays for unconscious at dch or not
## dimensions are now 28x28 because the L&R Midbrain locations will be removed
fa_unconsc <- array(data=NA, dim=c(28,28,0))
fa_consc <- array(data=NA, dim=c(28,28,0))
md_unconsc <- array(data=NA, dim=c(28,28,0))
md_consc <- array(data=NA, dim=c(28,28,0))
ad_unconsc <- array(data=NA, dim=c(28,28,0))
ad_consc <- array(data=NA, dim=c(28,28,0))
rd_unconsc <- array(data=NA, dim=c(28,28,0))
rd_consc <- array(data=NA, dim=c(28,28,0))
fiber_unconsc <- array(data=NA, dim=c(28,28,0))
fiber_consc <- array(data=NA, dim=c(28,28,0))





# number of patients
n_unconsc_MRI <- nrow(filter(conscdat, DTI.Number %in% ID.u & unconsious.at.time.of.MRI==1))
n_unconsc_Dch <- nrow(filter(conscdat, DTI.Number %in% ID.u & unconscious.on.discharge==1))

n_consc_MRI <- nrow(filter(conscdat, DTI.Number %in% ID.u & unconsious.at.time.of.MRI==0))
n_consc_Dch <- nrow(filter(conscdat, DTI.Number %in% ID.u & unconscious.on.discharge==0))


# look at these patients who are unconscious at MRI and see what they look like at discharge
n_unconsc_Dch_MRI <- nrow(filter(conscdat, DTI.Number %in% ID.u & unconsious.at.time.of.MRI==1 & unconscious.on.discharge==1))
n_consc_Dch_MRI <- nrow(filter(conscdat, DTI.Number %in% ID.u & unconsious.at.time.of.MRI==1 & unconscious.on.discharge==0))




## Now go through each ID,
## then go through each datatype,
## then calculate the matrix,
## then abind it to the corresponding array
## then take average on 3rd dimension of these arrays
## with apply: apply(fa_array, c(1,2), mean, na.rm=T)
## apparently need to specify the 1st and 2nd to average over the 3rd 
# (I guess it means 'for each element in 1st and 2nd dimensions, average over the 3rd)
#### Go through each patient, find their files, then create matrices for each datatype ####
for (i in 1:length(ID.u)) {
  print(paste("Processing ID: ", ID.u[i], sep=""))
  
  patfiles <- grep(ID.u[i], analysis2, value=T)
  
  
  for (j in 1:length(datatypes)) {
    
    typefile <- grep(datatypes[j], patfiles, value=T)
    
    # check if patient had file for this datatype
    if (length(typefile) == 0) {
      print(paste(ID.u[i], "doesn't have file for type:", datatypes[j]))
      next
    }
    
    dat <- read.table(typefile)[-1,-1]
    dat2 <- dat[-1,-1]
    # very complicated way to get the names first row, turn from data.frame to a vector, then turn only these names to factor levels (because this keeps the levels for all the different data types in each column)
    colnames(dat2) <- as.factor(as.character(unname(unlist(dat[1,][-1]))))
    rownames(dat2) <- dat[,1][-1]
    
    datmat <- as.matrix(dat2)
    # turn to numeric
    datmat.n <- apply(datmat, c(1,2), as.numeric)
    # remove lower half
    #datmat.n[lower.tri(datmat.n)] <- NA
    
    # change names with series of chained gsubs using regular expressions
    # to move the "_L/R" to the front of the strings
    colnames(datmat.n) <- gsub("Right", "mR",  
                                  gsub("Left", "mL", 
                                        gsub("^(.*)_R$", "R_\\1", 
                                            gsub("^(.*)_L$", "L_\\1", colnames(datmat.n)))))
    
    # change rownames the same way
    rownames(datmat.n) <- gsub("Right", "mR",  
                               gsub("Left", "mL", 
                                    gsub("^(.*)_R$", "R_\\1", 
                                         gsub("^(.*)_L$", "L_\\1", rownames(datmat.n)))))
    
    
    ## Remove the mR_Brainstem and mL_Brainstem columns
    datmat.n <- datmat.n[!rownames(datmat.n) %in% c("mL_Brainstem", "mR_Brainstem"), 
             !colnames(datmat.n) %in% c("mL_Brainstem", "mR_Brainstem")]
    
    
    ## reorder the matrix alphabetically
    ord <- corrMatOrder(datmat.n, order="alphabet")
    datmat.n <- datmat.n[ord, ord]
    
    ## or, instead, order by a set list of variable names
    #datmat.n[order(rownames(datmat.n))]
    var.order <- colnames(datmat.n) #VECTOR OF NAMES SHOULD GO HERE
    #View(datmat.n[var.order, var.order])
    
    # remove upper half (for heatmap plot. if plot in ggplot, then this removes the lower half, as expected.)
    #datmat.n[lower.tri(datmat.n)] <- NA
    # for ggplot plotting
    datmat.n[upper.tri(datmat.n)] <- NA
    
    # place into array depending on datatype and conscious status
    if (checkUnconscMRI == 1) {
      conscstatMRI <- filter(conscdat, DTI.Number == ID.u[i])[,"unconsious.at.time.of.MRI"] # Want this to equal 1
      conscstatDch <- filter(conscdat, DTI.Number == ID.u[i])[,"unconscious.on.discharge"] # Then want to look at patients with 0 and 1 for this
    } else {
      conscstat <- filter(conscdat, DTI.Number == ID.u[i])[,conscious_time]#$unconscious.on.discharge
    }
    
    
    
    # Do the below, but only for patients who were unconsciuous at time of MRI
    if (checkUnconscMRI == 1) {
      
      if (conscstatMRI == 1 & conscstatDch == 0) {  # Unconscious who recovered
        
        if (datatypes[j] == "FA") {
          fa_consc <- abind(fa_consc, datmat.n, rev.along=1)
        } else if (datatypes[j] == "MD") {
          md_consc <- abind(md_consc, datmat.n, rev.along=1)
        } else if (datatypes[j] == "AD") {
          ad_consc <- abind(ad_consc, datmat.n, rev.along=1)
        } else if (datatypes[j] == "RD") {
          rd_consc <- abind(rd_consc, datmat.n, rev.along=1)
        } else if (datatypes[j] == "fibercount") {
          fiber_consc <- abind(fiber_consc, datmat.n, rev.along=1)
        } else {
          stop("somehow, none of the datatypes were chosen in this loop...")
        }
        
      } else if (conscstatMRI == 1 & conscstatDch == 1) {  # Unconscious who didn't recover
        
        if (datatypes[j] == "FA") {
          fa_unconsc <- abind(fa_unconsc, datmat.n, rev.along=1)
        } else if (datatypes[j] == "MD") {
          md_unconsc <- abind(md_unconsc, datmat.n, rev.along=1)
        } else if (datatypes[j] == "AD") {
          ad_unconsc <- abind(ad_unconsc, datmat.n, rev.along=1)
        } else if (datatypes[j] == "RD") {
          rd_unconsc <- abind(rd_unconsc, datmat.n, rev.along=1)
        } else if (datatypes[j] == "fibercount") {
          fiber_unconsc <- abind(fiber_unconsc, datmat.n, rev.along=1)
        } else {
          stop("somehow, none of the datatypes were chosen in this loop...")
        }
        
      }
      
      
    } else {
      
      
      ## Do this when not checking for just the unconscious at MRI patients
      
      if (conscstat == 0) {
        
        if (datatypes[j] == "FA") {
          fa_consc <- abind(fa_consc, datmat.n, rev.along=1)
        } else if (datatypes[j] == "MD") {
          md_consc <- abind(md_consc, datmat.n, rev.along=1)
        } else if (datatypes[j] == "AD") {
          ad_consc <- abind(ad_consc, datmat.n, rev.along=1)
        } else if (datatypes[j] == "RD") {
          rd_consc <- abind(rd_consc, datmat.n, rev.along=1)
        } else if (datatypes[j] == "fibercount") {
          fiber_consc <- abind(fiber_consc, datmat.n, rev.along=1)
        } else {
          stop("somehow, none of the datatypes were chosen in this loop...")
        }
        
      } else if (conscstat == 1) {
        
        if (datatypes[j] == "FA") {
          fa_unconsc <- abind(fa_unconsc, datmat.n, rev.along=1)
        } else if (datatypes[j] == "MD") {
          md_unconsc <- abind(md_unconsc, datmat.n, rev.along=1)
        } else if (datatypes[j] == "AD") {
          ad_unconsc <- abind(ad_unconsc, datmat.n, rev.along=1)
        } else if (datatypes[j] == "RD") {
          rd_unconsc <- abind(rd_unconsc, datmat.n, rev.along=1)
        } else if (datatypes[j] == "fibercount") {
          fiber_unconsc <- abind(fiber_unconsc, datmat.n, rev.along=1)
        } else {
          stop("somehow, none of the datatypes were chosen in this loop...")
        }
        
      } else {
        stop("conscious status isn't defined I guess?")
      }
      
      
    }
    
    
    
    
  }
}




### Average across patients
fa_unconsc_avg <-  apply(fa_unconsc, c(1,2), mean, na.rm=T)
fa_consc_avg <-  apply(fa_consc, c(1,2), mean, na.rm=T)
md_unconsc_avg <-  apply(md_unconsc, c(1,2), mean, na.rm=T)
md_consc_avg <-  apply(md_consc, c(1,2), mean, na.rm=T)
ad_unconsc_avg <-  apply(ad_unconsc, c(1,2), mean, na.rm=T)
ad_consc_avg <-  apply(ad_consc, c(1,2), mean, na.rm=T)
rd_unconsc_avg <-  apply(rd_unconsc, c(1,2), mean, na.rm=T)
rd_consc_avg <-  apply(rd_consc, c(1,2), mean, na.rm=T)
fiber_unconsc_avg <-  apply(fiber_unconsc, c(1,2), mean, na.rm=T)
fiber_consc_avg <-  apply(fiber_consc, c(1,2), mean, na.rm=T)


### Deltas between unconscious and conscious
#fa_delta <- fa_unconsc_avg - fa_consc_avg
#fiber_delta1 <- fiber_unconsc_avg - fiber_consc_avg

fiber_delta <- fiber_consc_avg - fiber_unconsc_avg
# percent change -- doesn't work b/c have some 0 values in fiber_unconsc_avg
#fiber_delta2 <- ((fiber_consc_avg - fiber_unconsc_avg) / fiber_unconsc_avg) * 100

# use this other measure? divides by average?
#fiber_delta2 <- (fiber_consc_avg - fiber_unconsc_avg) / ((fiber_consc_avg + fiber_unconsc_avg)/2)
###### So, ANY change from 0 will always result in a 2 or -2.
###### e.g.: (-0.000000000000000000000000001-0)/((0.000000000000000000000000001-0)/2)
###### ...this equals -2.

# percent change adding 1 to denominator (this one is better)
fiber_delta3 <- ((fiber_consc_avg - fiber_unconsc_avg) / (fiber_unconsc_avg+1)) * 100


# normalize between -1 and 1
# no this doesn't do anything--all it does it change the scale.
#fiber_delta2 <- fiber_consc_avg - fiber_unconsc_avg
#fiber_delta2 <- 2 * ( (fiber_delta2 - min(fiber_delta2, na.rm=T)) / (max(fiber_delta2, na.rm=T) - min(fiber_delta2, na.rm=T)) ) - 1






### Melt matrices in order to plot in ggplot
fiber_unconsc_m <- melt(fiber_unconsc_avg, na.rm=T)
fiber_consc_m <- melt(fiber_consc_avg, na.rm=T)


#fiber_delta_m1 <- melt(fiber_delta1, na.rm=T)
fiber_delta_m <- melt(fiber_delta, na.rm=T)
#fiber_delta2_m <- melt(fiber_delta2, na.rm=T)
fiber_delta3_m <- melt(fiber_delta3, na.rm=T)

# normalize values - feature scaling
#fiber_delta_m2$value <- (fiber_delta_m2$value - min(fiber_delta_m2$value)) / (max(fiber_delta_m2$value) - min(fiber_delta_m2$value))

# standardize values
#fiber_delta_m2$value <- fiber_delta_m2$value - (mean(fiber_delta_m2$value, na.rm=T) / sd(fiber_delta_m2$value, na.rm=T))

# find percent change
#fiber_delta_m2$value <- 







#### Make plots for all patients who were unconscious at MRI ####
if (conscious_time == "unconsious.at.time.of.MRI" & checkUnconscMRI != 1) {

  ## unconscious
  p1 <- ggplot(fiber_unconsc_m, aes(x=Var1, y=Var2, fill=value)) + 
          geom_tile(color="white") +
          scale_fill_gradient2(low="blue", high="red", mid="white",
                               midpoint = 0, limit=c(0, 895), space="Lab", name="num fibers") +   #c(-460, 420)  c(-82, 1000)
          theme_minimal() +
          xlab("") + ylab("") +  # make axis labels blank
          theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "left") +
          scale_y_discrete(position = "right") +
          ggtitle("Fiber Count", subtitle="Unconscious at MRI    n=7") +
          coord_fixed()
  
  ## conscious
  p2 <- ggplot(fiber_consc_m, aes(x=Var1, y=Var2, fill=value)) + 
          geom_tile(color="white") +
          scale_fill_gradient2(low="blue", high="red", mid="white",
                               midpoint = 0, limit=c(0, 895), space="Lab", name="num fibers") +   #c(-460, 420)  c(-82, 1000) 0, 730
          theme_minimal() +
          xlab("") + ylab("") +  # make axis labels blank
          theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "left") +
          scale_y_discrete(position = "right") +
          ggtitle("Fiber Count", subtitle="Conscious at MRI    n=7") +
          coord_fixed()
  
  
  p3 <- ggplot(fiber_delta_m, aes(x=Var1, y=Var2, fill=value)) + 
          geom_tile(color="white") +
          scale_fill_gradient2(low="blue", high="red", mid="white",
                               midpoint = 0, limit=c(-460, 420), space="Lab", name="delta fiber num") +   #c(-460, 420)  c(-82, 1000)
          theme_minimal() +
          xlab("") + ylab("") +  # make axis labels blank
          theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "left") +
          scale_y_discrete(position = "right") +
          ggtitle("Fiber Delta: consc - unconsc", subtitle="at time of MRI") +
          coord_fixed()
  
  
  # p4 <- ggplot(fiber_delta2_m, aes(x=Var1, y=Var2, fill=value)) + 
  #         geom_tile(color="white") +
  #         scale_fill_gradient2(low="blue", high="red", mid="white",
  #                              midpoint = 0, limit=c(-2, 2), space="Lab", name="% change") +   #c(-460, 420)  c(-82, 1000)
  #         theme_minimal() +
  #         xlab("") + ylab("") +  # make axis labels blank
  #         theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "left") +
  #         scale_y_discrete(position = "right") +
  #         ggtitle("Fiber Delta % Difference: (consc - unconsc) / ((consc + unconsc)/2)", subtitle="at time of MRI") +
  #         coord_fixed()
  
  
  p5 <- ggplot(fiber_delta3_m, aes(x=Var1, y=Var2, fill=value)) + 
          geom_tile(color="white") +
          scale_fill_gradient2(low="blue", high="red", mid="white",
                               midpoint = 0, limit=c(-95, 3575), space="Lab", name="% change") +   #c(-460, 420)  c(-82, 1000)
          theme_minimal() +
          xlab("") + ylab("") +  # make axis labels blank
          theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "left") +
          scale_y_discrete(position = "right") +
          ggtitle("Fiber Delta % Difference: (consc - unconsc) / (uncosc + 1) * 100", subtitle="at time of MRI") +
          coord_fixed()
  
  
  print(p1)
  print(p2)
  print(p3)
  #print(p4)
  print(p5)
  
  if (savePlots == 1) {
    ggsave("../new heatmaps/unconscious at MRI.tiff", p1, device="tiff", dpi=600, width=8, height=8, units="in")
    ggsave("../new heatmaps/conscious at MRI.tiff", p2, device="tiff", dpi=600, width=8, height=8, units="in")
    ggsave("../new heatmaps/delta at MRI.tiff", p3, device="tiff", dpi=600, width=8, height=8, units="in")
    #ggsave("../new heatmaps/percent change at MRI.tiff", p4, device="tiff", dpi=600)
    ggsave("../new heatmaps/percent change at MRI.tiff", p5, device="tiff", dpi=600, width=8, height=8, units="in")
  }
  


}






#### Make plots for all patients who were unconscious at discharge ####
if (conscious_time == "unconscious.on.discharge" & checkUnconscMRI != 1) {
  
  ## unconscious
  p1 <- ggplot(fiber_unconsc_m, aes(x=Var1, y=Var2, fill=value)) + 
          geom_tile(color="white") +
          scale_fill_gradient2(low="blue", high="red", mid="white",
                               midpoint = 0, limit=c(0, 895), space="Lab", name="num fibers") +   #c(-460, 420)  c(-82, 1000)
          theme_minimal() +
          xlab("") + ylab("") +  # make axis labels blank
          theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "left") +
          scale_y_discrete(position = "right") +
          ggtitle("Fiber Count", subtitle="Unconscious at Discharge    n=4") +
          coord_fixed()
  
  ## conscious
  p2 <- ggplot(fiber_consc_m, aes(x=Var1, y=Var2, fill=value)) + 
          geom_tile(color="white") +
          scale_fill_gradient2(low="blue", high="red", mid="white",
                               midpoint = 0, limit=c(0, 895), space="Lab", name="num fibers") +   #c(-460, 420)  c(-82, 1000)
          theme_minimal() +
          xlab("") + ylab("") +  # make axis labels blank
          theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "left") +
          scale_y_discrete(position = "right") +
          ggtitle("Fiber Count", subtitle="Conscious at Discharge    n=10") +
          coord_fixed()
  
  
  p3 <- ggplot(fiber_delta_m, aes(x=Var1, y=Var2, fill=value)) + 
          geom_tile(color="white") +
          scale_fill_gradient2(low="blue", high="red", mid="white",
                               midpoint = 0, limit=c(-460, 420), space="Lab", name="delta fiber num") +   #c(-460, 420)  c(-82, 1000)
          theme_minimal() +
          xlab("") + ylab("") +  # make axis labels blank
          theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "left") +
          scale_y_discrete(position = "right") +
          ggtitle("Fiber Delta: consc - unconsc", subtitle="at Discharge") +
          coord_fixed()
  
  
  # p4 <- ggplot(fiber_delta2_m, aes(x=Var1, y=Var2, fill=value)) + 
  #         geom_tile(color="white") +
  #         scale_fill_gradient2(low="blue", high="red", mid="white",
  #                              midpoint = 0, limit=c(-2, 2), space="Lab", name="% change") +   #c(-460, 420)  c(-82, 1000)
  #         theme_minimal() +
  #         xlab("") + ylab("") +  # make axis labels blank
  #         theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "left") +
  #         scale_y_discrete(position = "right") +
  #         ggtitle("Fiber Delta % Difference: consc - unconsc", subtitle="at Discharge") +
  #         coord_fixed()
  
  
  p5 <- ggplot(fiber_delta3_m, aes(x=Var1, y=Var2, fill=value)) + 
          geom_tile(color="white") +
          scale_fill_gradient2(low="blue", high="red", mid="white",
                               midpoint = 0, limit=c(-95, 3575), space="Lab", name="% change") +   #c(-460, 420)  c(-82, 1000)
          theme_minimal() +
          xlab("") + ylab("") +  # make axis labels blank
          theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "left") +
          scale_y_discrete(position = "right") +
          ggtitle("Fiber Delta % Difference: (consc - unconsc) / (uncosc + 1) * 100", subtitle="at Discharge") +
          coord_fixed()
        
  
  print(p1)
  print(p2)
  print(p3)
  #print(p4)
  print(p5)
  
  if (savePlots == 1) {
    ggsave("../new heatmaps/unconscious at Dch.tiff", p1, device="tiff", dpi=600, width=8, height=8, units="in")
    ggsave("../new heatmaps/conscious at Dch.tiff", p2, device="tiff", dpi=600, width=8, height=8, units="in")
    ggsave("../new heatmaps/delta at Dch.tiff", p3, device="tiff", dpi=600, width=8, height=8, units="in")
    #ggsave("../new heatmaps/percent change at Dch.tiff", p4, device="tiff", dpi=600)
    ggsave("../new heatmaps/percent change at Dch.tiff", p5, device="tiff", dpi=600, width=8, height=8, units="in")
  }
  
  
  
}






#### Make plots for all patients who were unconscious at MRI and whether they were conscious at discharge ####
if (checkUnconscMRI == 1) {
  
  ## unconscious
  p1 <- ggplot(fiber_unconsc_m, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile(color="white") +
    scale_fill_gradient2(low="blue", high="red", mid="white",
                         midpoint = 0, limit=c(0, 895), space="Lab", name="num fibers") +   #c(-460, 420)  c(-82, 1000)
    theme_minimal() +
    xlab("") + ylab("") +  # make axis labels blank
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "left") +
    scale_y_discrete(position = "right") +
    ggtitle("Fiber Count", subtitle="Unconscious at MRI and Remain Unconscious at Discharge    n=4") +
    coord_fixed()
  
  ## conscious
  p2 <- ggplot(fiber_consc_m, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile(color="white") +
    scale_fill_gradient2(low="blue", high="red", mid="white",
                         midpoint = 0, limit=c(0, 895), space="Lab", name="num fibers") +   #c(-460, 420)  c(-82, 1000)
    theme_minimal() +
    xlab("") + ylab("") +  # make axis labels blank
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "left") +
    scale_y_discrete(position = "right") +
    ggtitle("Fiber Count", subtitle="Unconscious at MRI and Recover Conscious at Discharge    n=3") +
    coord_fixed()
  
  
  p3 <- ggplot(fiber_delta_m, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile(color="white") +
    scale_fill_gradient2(low="blue", high="red", mid="white",
                         midpoint = 0, limit=c(-530, 530), space="Lab", name="delta fiber num") +   #c(-460, 420)  c(-82, 1000)
    theme_minimal() +
    xlab("") + ylab("") +  # make axis labels blank
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "left") +
    scale_y_discrete(position = "right") +
    ggtitle("Fiber Delta: consc - unconsc", subtitle="Outcome after unconsc at MRI") +
    coord_fixed()
  
  
  # p4 <- ggplot(fiber_delta2_m, aes(x=Var1, y=Var2, fill=value)) + 
  #   geom_tile(color="white") +
  #   scale_fill_gradient2(low="blue", high="red", mid="white",
  #                        midpoint = 0, limit=c(-2, 2), space="Lab", name="% change") +   #c(-460, 420)  c(-82, 1000)
  #   theme_minimal() +
  #   xlab("") + ylab("") +  # make axis labels blank
  #   theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "left") +
  #   scale_y_discrete(position = "right") +
  #   ggtitle("Fiber Delta % Difference: consc - unconsc / average", subtitle="Outcome after unconsc at MRI") +
  #   coord_fixed()
  
  
  p5 <- ggplot(fiber_delta3_m, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile(color="white") +
    scale_fill_gradient2(low="blue", high="red", mid="white",
                         midpoint = 0, limit=c(-1300, 1300), space="Lab", name="% change") +   #c(-460, 420)  c(-82, 1000)
    theme_minimal() +
    xlab("") + ylab("") +  # make axis labels blank
    theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "left") +
    scale_y_discrete(position = "right") +
    ggtitle("Fiber Delta % Difference: (consc - unconsc) / (uncosc + 1) * 100", subtitle="Outcome after unconsc at MRI") +
    coord_fixed()
  
  
  print(p1)
  print(p2)
  print(p3)
  #print(p4)
  print(p5)
  
  if (savePlots == 1) {
    ggsave("../new heatmaps/UNCONSC MRI unconscious at Dch.tiff", p1, device="tiff", dpi=600, width=8, height=8, units="in")
    ggsave("../new heatmaps/UNCONSC MRI conscious at Dch.tiff", p2, device="tiff", dpi=600, width=8, height=8, units="in")
    ggsave("../new heatmaps/UNCONSC MRI delta at Dch.tiff", p3, device="tiff", dpi=600, width=8, height=8, units="in")
    #ggsave("../new heatmaps/UNCONSC MRI percent change at Dch.tiff", p4, device="tiff", dpi=600, width=8, height=8, units="in")
    # This seems to work to produce a decent tiff w/ dpi=600 and still has the color legend bar
    ggsave("../new heatmaps/UNCONSC MRI percent change at Dch.tiff", p5, device="tiff", dpi=600, width=8, height=8, units="in")
  }
  
  
  
}






