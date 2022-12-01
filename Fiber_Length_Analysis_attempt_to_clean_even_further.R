library(openxlsx)
library(ggplot2)
library(dplyr)

# Extract data from the fiber length text file output from
# whatever that program Batool used in the first place.
ParseFiberLenData <- function(txtfiles, lens, file_dir, control=FALSE, conscdat=NULL) {
  library(reshape2)
  library(dplyr)
  # initialize empty table
  mainTable <- data.frame()
  # empty list to put all tables for each data type in
  out <- list()
  # parse text files
  split_txt <- strsplit(txtfiles, "_")
  
  # extract "columns"
  ID <- sapply(split_txt, "[", 1)
  ID.u <- unique(ID)
  
  ### loops to extract data from text files
  # first loops over all unique patient IDs,
  # then loops over each length-type for each patient,
  # then loops over each file for each length-type
  # (there normally should just be one file per length-type anyway)
  for (i in 1:length(ID.u)) {
    print(paste("Processing ID: ", ID.u[i], sep=""))
    
    patfiles <- grep(ID.u[i], txtfiles, value=T)
    
    
    for (j in 1:length(lens)) {
      #print(paste("for length: ", lens[j], sep=""))
      
      lenfiles <- grep(lens[j], patfiles, value=T)
      
      if (length(lenfiles) == 0) {
        print(paste(ID.u[i], "doesn't have file for type:", lens[j]))
        next
      }
      
      
      for (file in lenfiles) {  # this should always be of length==1 anyway
        file_path <- paste(file_dir, "/", file, sep="")
        dat <- try(read.table(file_path, skip=6), silent=TRUE)  # there's also a blank line at the beginning of the file to skip
        if (class(dat) == "try-error") {
          print(paste("something wrong with file:", file, "-- skipping over it"))
          next
        }
        colnames(dat) <- c("param", "calc", "value")
        # I guess I can cast the data
        # also apparently NOW I don't need an aggragate function for casting??
        dat.c <- dcast(dat, param ~ calc)
        # add in a column for the patient ID and length-type
        dat.c$ID <- ID.u[i]
        dat.c$len <- lens[j]
        # then append this table to the main table
        mainTable <- rbind(mainTable, dat.c)
        
        
      }
      
    }
    
  }
  
  # make lengths factors
  mainTable$len <- as.factor(mainTable$len)
  
  # join consciousness data if "control" == FALSE
  # "for now", this used the hard-coded names in the conscdat table.
  # it probably will never change, so be wary if using another data source for consciousness info.
  # lololol and that just happened to me! Now I have to add in the 3 new groups manually here.
  if (control == FALSE) {
    mainTable <- inner_join(mainTable, 
                            select(conscdat, DTI.Number, unconsious.on.admission, unconsious.at.time.of.MRI, unconscious.on.discharge, Recovery.of.consciousness, conscious.throughout, MRI.conscious.recover, MRI.conscious.no.recover, group),
                            by=c("ID" = "DTI.Number"))
    
  } else {
    # actually want to add the same columns to the control data, except make them all NAs
    mainTable <- mainTable %>%
      mutate(unconsious.on.admission=NA, unconsious.at.time.of.MRI=NA, unconscious.on.discharge=NA, Recovery.of.consciousness=NA, conscious.throughout=NA, MRI.conscious.recover=NA, MRI.conscious.no.recover=NA, group="control")
    
  }
  
  
  ### save each parameter type to its own data frame
  
  # first get all unique data type names
  datnames <- unique(mainTable$param)
  # use these later
  newnms <- paste(datnames, "Dat", sep="")
  
  # then assign in a loop
  for (name in datnames) {
    assign(
      paste(name, "Dat", sep=""),
      mainTable %>%
        filter(param == name & len != "L10_L80")  # don't care about the overall huge length
    )
    
  }
  
  
  # calculate ratios for the L10_L20 and L80_L160 lengths and save to output list.
  # I guess I could have just put the entire dplyr part into the list instead of assigning first.
  # oh well.
  for (name in newnms) {
    # calculate ratio
    assign(
      name,
      eval(as.name(name)) %>%  # sometimes this works without as.name(). No idea why.
        group_by(ID) %>%
        mutate(ratio.80_60.10_20 = mean[len=="L80_L160"] / mean[len=="L10_L20"])
    )
    
    # put into output list
      out[[name]] <- eval(as.name(name))
    
  }
  
  # also return entire data set
  out[["MainTable"]] <- mainTable
  
  return(out)
  
}



## Make plots as a function and save them all to a list.
## Can work in a way to skip over the MainTable plots, cause they're useless.
## For now, can just ignore them.
PlotEachDataType <- function(datalist, consc_groups, control=FALSE) {
  out <- list()
  names_to_use <- names(datalist)[names(datalist) != "MainTable"]
  for (name in names_to_use) {
    # make plot names
    if (control == FALSE) {
      #plt_names <- paste(name, "_", consc_groups, sep="")
      plt_names <- name
    } else {
      plt_names <- paste(name, "Control")
    }
    
    # assign plots
    for (i in 1:length(plt_names)) {
      if (control==FALSE) {
        # this uses an obnoxious way to remove rows missing data for the current conscious group
        out[[plt_names[i]]] <- ggplot(datalist[[name]][!is.na(datalist[[name]][consc_groups[i]]),], aes_string(x="len", y="mean", group="ID", color=consc_groups[i])) +
          geom_line(size=1.1) + theme_bw() +
          labs(title=paste(name, consc_groups[i]), x="length", y="mean", color="key")
        
      } else {
        
        out[[plt_names[i]]] <- ggplot(datalist[[name]], aes_string(x="len", y="mean", group="ID")) +
          geom_line(size=1.1) + theme_bw() +
          labs(title=paste("Control", name), x="length", y="mean", color="key")
        
      }
    }
  }
  return(out)
}



removeDiag <- function(df, rm.region) {
  require(reshape2)
  # don't really need to specify value.var in dcast, 
  #as it will just take the average over one value, which just equals the value anyway
  # this defaults to use the "delta" column for the delta tables, which is fortunate, so I'm just going with it.
  mat <- as.matrix(dcast(df, Var1 ~ Var2))
  #dat2 <- dat[-1,-1]
  mat2 <- mat[,-1]
  #rownames(dat2) <- dat[,1][-1]
  rownames(mat2) <- mat[,1]
  
  # CONVERT VALUES BACK TO NUMERIC. They all turn to character due to the rownames starting out as the first column
  mat2 <- apply(mat2, MARGIN=c(1,2), FUN=as.numeric)
  
  if (rm.region == "upper") {
    mat2[upper.tri(mat2)] <- NA
    
  } else if (rm.region == "lower") {
    mat2[lower.tri(mat2)] <- NA
    
  } else {
    stop("Must specify 'upper' or 'lower' for rm.region")
  }
  
  # remelt matrix with triangle removed
  mat.m <- melt(mat2, na.rm=T)
  mat.m$value <- as.numeric(mat.m$value)
  #colnames(mat.m) <- c("Var1", "Var2", "value")
  
  return(mat.m)
  
}



### Need to make a function to do the fibercount heatmaps too
### Take the code from "Just lesion loc Plot_Matrix.R" and turn into a function
CreateConnectivityMatrix <- function(txtfiles, file_dir, datatypes, control=FALSE, checkUnconscMRI=FALSE, conscious_time=NULL, conscdat=NULL, locdat=NULL) {
  library(corrplot)
  
  # make empty dataframes, depending on if analyzing controls or not
  if (control==FALSE) {
    fa_unconsc.df <- data.frame()
    fa_consc.df <- data.frame()
    md_unconsc.df <- data.frame()
    md_consc.df <- data.frame()
    ad_unconsc.df <- data.frame()
    ad_consc.df <- data.frame()
    rd_unconsc.df <- data.frame()
    rd_consc.df <- data.frame()
    fiber_unconsc.df <- data.frame()
    fiber_consc.df <- data.frame()
  } else {
    fa_control.df <- data.frame()
    md_control.df <- data.frame()
    ad_control.df <- data.frame()
    rd_control.df <- data.frame()
    fiber_control.df <- data.frame()
  }
  
  # list to put all dataframes into
  out_main <- list()
 
  
  # get only the "Analysis2" files, which are the connectivity matrices
  analysis2 <- grep("Analysis2", txtfiles, value=T)
  # extract unique patient identifiers
  split_txt <- strsplit(analysis2, "_")
  ID <- sapply(split_txt, "[", 1)
  ID.u <- unique(ID)
  
  # for each patient...
  for (i in 1:length(ID.u)) {
    print(paste("Processing ID: ", ID.u[i], sep=""))
    
    patfiles <- grep(ID.u[i], analysis2, value=T)
    
    # go through all of the connectivity matrix files for them
    for (j in 1:length(datatypes)) {
      
      typefile <- grep(datatypes[j], patfiles, value=T)
      
      # check if patient had file for this datatype
      if (length(typefile) == 0) {
        print(paste(ID.u[i], "doesn't have file for type:", datatypes[j]))
        next
      }
      
      file_path <- paste(file_dir, "/", typefile, sep="")
      dat <- try(read.table(file_path)[-1,-1], silent=TRUE)
      if (class(dat) == "try-error") {
        print(paste("something wrong with file:", file, "-- skipping over it"))
        next
      }
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
      ### we never used this?
      ###var.order <- labelorder#colnames(datmat.n) #VECTOR OF NAMES SHOULD GO HERE
      #View(datmat.n[var.order, var.order])
      
      ### Change variable names to have "ipsilateral" or "contralateral"
      ### based off of the injury location
      
      ### ....We have a problem, seeing as the injury location table
      ### uses MRNs to denote patients, while the consciousness status table
      ### AND the filenames use the DTI (ICH) number to denote patients.
      
      # MRNS ARE IN THE FILE AFTER ALL.
      # find the mrn for this DTI number
      # do this only for the not-controls
      if (control == FALSE) {
        
        curr_MRN <- unique(filter(conscdat, DTI.Number == ID.u[i])$MRN)
        
        # once I can check the MRN to DTI number, though, I could do something like this:
        sub_tent <- filter(locdat, MRN == curr_MRN)$sub_tent
        sus_tent <- filter(locdat, MRN == curr_MRN)$sus_tent
        
        print(paste("Patient", ID.u[i], "has lesion on", sus_tent))
        
        # get col and row name indices
        # the row and column indices should probably always be the same anyway
        colixL <- grep("^L_", colnames(datmat.n), value=F)
        rowixL <- grep("^L_", rownames(datmat.n), value=F)
        colixR <- grep("^R_", colnames(datmat.n), value=F)
        rowixR <- grep("^R_", rownames(datmat.n), value=F)
        
        # check if we have location data for the current patient
        if (!curr_MRN %in% locdat$MRN) {
          print(paste("Missing location data for patient:", curr_MRN, "/", ID.u[i]))
          break  # because I tell it to skip to next patient, ICH16 is STILL not included in the analysis.
        }
        
        if (sus_tent == "L") {
          # (the caret indicates that this should be found only at the beginning of the string)
          # ipsilateral
          colnames(datmat.n)[colixL] <- gsub("^L_", "ipsilateral_", colnames(datmat.n)[colixL])
          rownames(datmat.n)[rowixL] <- gsub("^L_", "ipsilateral_", rownames(datmat.n)[rowixL])
          # contralateral
          colnames(datmat.n)[colixR] <- gsub("^R_", "contralateral_", colnames(datmat.n)[colixR])
          rownames(datmat.n)[rowixR] <- gsub("^R_", "contralateral_", rownames(datmat.n)[rowixR])
          
        } else if (sus_tent == "R") {
          # ipsilateral
          colnames(datmat.n)[colixR] <- gsub("^R_", "ipsilateral_", colnames(datmat.n)[colixR])
          rownames(datmat.n)[rowixR] <- gsub("^R_", "ipsilateral_", rownames(datmat.n)[rowixR])
          # contralateral
          colnames(datmat.n)[colixL] <- gsub("^L_", "contralateral_", colnames(datmat.n)[colixL])
          rownames(datmat.n)[rowixL] <- gsub("^L_", "contralateral_", rownames(datmat.n)[rowixL])
          
        } else {
          next  #### Need to find out what to do when lesion is "Both" or "None". Just go to next loop iteration for now (might be able to make this a break statement to go to next MRN loop)
        }
        
      }
       #### Also need to find out what to do with the sub_tent variable.
      
      ###=====================================###
      ### so then MAYBE I can (instead of storing these matrices as the 3rd dimension in the big arrays)
      ### just melt the matrix *now* and then append it to a big data frame.
      ### Then later, when I take the average across patients, 
      ### I can just group_by DTI/ICH/MRN/whatever AND the "variables" (i.e., the matrix columns and rows)
      ### and then mutate to find the mean of these.
      ###=====================================##
      # DON'T REMOVE THE UPPER HALF FIRST when making the melted temp data frame
      ## Because: one patient might have (e.g.,) "ipsilateral_Amygdala" and "contralateral_Amydala"
      ## as var1 and var2 respectivelly, but if another patient has the ipsi and contra flipped (b/c lesion on other side),
      ## then THIS patient's var1 and var2 for these same names will be removed, because this combo now occurs in the other tri
      ## Can just remove duplicates from the final data frame later or something
      
      # melt into temp df
      tmp <- melt(datmat.n, na.rm=T)
      # add in column for MRN (if not control) and DTI number
      if (control==FALSE) {
        tmp$MRN <- curr_MRN
        tmp$DTI <- ID.u[i]
      } else {
        tmp$DTI <- ID.u[i]
      }
      
      
      
      
      # remove upper half (for heatmap plot. if plot in ggplot, then this removes the lower half, as expected.)
      #datmat.n[lower.tri(datmat.n)] <- NA
      # for ggplot plotting
      #datmat.n[upper.tri(datmat.n)] <- NA
      
      # place into array depending on datatype and conscious status
      ## Do this only for non-controls
      ## (and probably want to change this to use the new group label)
      if (control == FALSE) {
        
        if (checkUnconscMRI == TRUE) {
          conscstatMRI <- filter(conscdat, DTI.Number == ID.u[i])[,"unconsious.at.time.of.MRI"] # Want this to equal 1
          conscstatDch <- filter(conscdat, DTI.Number == ID.u[i])[,"unconscious.on.discharge"] # Then want to look at patients with 0 and 1 for this
        } else {
          conscstat <- filter(conscdat, DTI.Number == ID.u[i])[,conscious_time]#$unconscious.on.discharge
        }
        
        
        
        # Do the below, but only for patients who were unconsciuous at time of MRI
        if (checkUnconscMRI == TRUE) {
          
          if (conscstatMRI == 1 & conscstatDch == 0) {  # Unconscious who recovered
            
            if (datatypes[j] == "FA") {
              
              # also append to data frame
              fa_consc.df <- rbind(fa_consc.df, tmp)
            } else if (datatypes[j] == "MD") {
              
              md_consc.df <- rbind(md_consc.df, tmp)
            } else if (datatypes[j] == "AD") {
              
              ad_consc.df <- rbind(ad_consc.df, tmp)
            } else if (datatypes[j] == "RD") {
              
              rd_consc.df <- rbind(rd_consc.df, tmp)
            } else if (datatypes[j] == "fibercount") {
              
              fiber_consc.df <- rbind(fiber_consc.df, tmp)
            } else {
              stop("somehow, none of the datatypes were chosen in this loop...")
            }
            
          } else if (conscstatMRI == 1 & conscstatDch == 1) {  # Unconscious who didn't recover
            
            if (datatypes[j] == "FA") {
              
              fa_unconsc.df <- rbind(fa_unconsc.df, tmp)
            } else if (datatypes[j] == "MD") {
              
              md_unconsc.df <- rbind(md_unconsc.df, tmp)
            } else if (datatypes[j] == "AD") {
              
              ad_unconsc.df <- rbind(ad_unconsc.df, tmp)
            } else if (datatypes[j] == "RD") {
              
              rd_unconsc.df <- rbind(rd_unconsc.df, tmp)
            } else if (datatypes[j] == "fibercount") {
              
              fiber_unconsc.df <- rbind(fiber_unconsc.df, tmp)
            } else {
              stop("somehow, none of the datatypes were chosen in this loop...")
            }
            
          }
          
          
        } else {
          
          
          ## Do this when not checking for just the unconscious at MRI patients
          
          if (conscstat == 0) {
            
            if (datatypes[j] == "FA") {
              
              fa_consc.df <- rbind(fa_consc.df, tmp)
            } else if (datatypes[j] == "MD") {
              
              md_consc.df <- rbind(md_consc.df, tmp)
            } else if (datatypes[j] == "AD") {
              
              ad_consc.df <- rbind(ad_consc.df, tmp)
            } else if (datatypes[j] == "RD") {
              
              rd_consc.df <- rbind(rd_consc.df, tmp)
            } else if (datatypes[j] == "fibercount") {
              
              fiber_consc.df <- rbind(fiber_consc.df, tmp)
            } else {
              stop("somehow, none of the datatypes were chosen in this loop...")
            }
            
          } else if (conscstat == 1) {
            
            if (datatypes[j] == "FA") {
              
              fa_unconsc.df <- rbind(fa_unconsc.df, tmp)
            } else if (datatypes[j] == "MD") {
              
              md_unconsc.df <- rbind(md_unconsc.df, tmp)
            } else if (datatypes[j] == "AD") {
              
              ad_unconsc.df <- rbind(ad_unconsc.df, tmp)
            } else if (datatypes[j] == "RD") {
              
              rd_unconsc.df <- rbind(rd_unconsc.df, tmp)
            } else if (datatypes[j] == "fibercount") {
              
              fiber_unconsc.df <- rbind(fiber_unconsc.df, tmp)
            } else {
              stop("somehow, none of the datatypes were chosen in this loop...")
            }
            
          } else {
            stop("conscious status isn't defined I guess?")
          }
          
          
        }
        
      } else {
        if (datatypes[j] == "FA") {
          fa_control.df <- rbind(fa_control.df, tmp)
          
        } else if (datatypes[j] == "MD") {
          md_control.df <- rbind(md_control.df, tmp)
          
        } else if (datatypes[j] == "AD") {
          ad_control.df <- rbind(ad_control.df, tmp)
          
        } else if (datatypes[j] == "RD") {
          rd_control.df <- rbind(rd_control.df, tmp)
          
        } else if (datatypes[j] == "fibercount") {
          fiber_control.df <- rbind(fiber_control.df, tmp)
          
        } else {
          stop("somehow, none of the datatypes were chosen in this loop...?")
        }
      }
      
      
      
      
      
      
    }
  }
  
  ## average across patients
  # ...I can't think of a good way to use lapply on a list of these dataframes
  # and perform the two dplyr functions on them..
  if (control==FALSE) {
    fa_unconsc_avg.df <- fa_unconsc.df %>%
      group_by(Var1, Var2) %>%
      summarise(avg = mean(value, na.rm=T))
    
    fa_consc_avg.df <- fa_consc.df %>%
      group_by(Var1, Var2) %>%
      summarise(avg = mean(value, na.rm=T))
    
    md_unconsc_avg.df <- md_unconsc.df %>%
      group_by(Var1, Var2) %>%
      summarise(avg = mean(value, na.rm=T))
    
    md_consc_avg.df <- md_consc.df %>%
      group_by(Var1, Var2) %>%
      summarise(avg = mean(value, na.rm=T))
    
    ad_unconsc_avg.df <- ad_unconsc.df %>%
      group_by(Var1, Var2) %>%
      summarise(avg = mean(value, na.rm=T))
    
    ad_consc_avg.df <- ad_consc.df %>%
      group_by(Var1, Var2) %>%
      summarise(avg = mean(value, na.rm=T))
    
    rd_unconsc_avg.df <- rd_unconsc.df %>%
      group_by(Var1, Var2) %>%
      summarise(avg = mean(value, na.rm=T))
    
    rd_consc_avg.df <- rd_consc.df %>%
      group_by(Var1, Var2) %>%
      summarise(avg = mean(value, na.rm=T))
    
    fiber_unconsc_avg.df <- fiber_unconsc.df %>%
      group_by(Var1, Var2) %>%
      summarise(avg = mean(value, na.rm=T))
    
    fiber_consc_avg.df <- fiber_consc.df %>%
      group_by(Var1, Var2) %>%
      summarise(avg = mean(value, na.rm=T))
    
    # remove diagonals using lapply after saving all dataframes to lists
    out1 <- list(fa_consc=fa_consc_avg.df, md_consc=md_consc_avg.df, ad_consc=ad_consc_avg.df, rd_consc=rd_consc_avg.df, fiber_consc=fiber_consc_avg.df)
    out2 <- list(fa_unconsc=fa_unconsc_avg.df, md_unconsc=md_unconsc_avg.df, ad_unconsc=ad_unconsc_avg.df, rd_unconsc=rd_unconsc_avg.df, fiber_unconsc=fiber_unconsc_avg.df)
    
    out1 <- lapply(out1, removeDiag, rm.region="upper")
    out2 <- lapply(out2, removeDiag, rm.region="upper")
    
    ## Deltas
    fa_delta <- fa_unconsc_avg.df %>%
      inner_join(fa_consc_avg.df, by=c("Var1", "Var2")) %>%
      mutate(delta = ((avg.y - avg.x) / (avg.x + 1)) * 100)
    
    md_delta <- md_unconsc_avg.df %>%
      inner_join(md_consc_avg.df, by=c("Var1", "Var2")) %>%
      mutate(delta = ((avg.y - avg.x) / (avg.x + 1)) * 100)
    
    ad_delta <- ad_unconsc_avg.df %>%
      inner_join(ad_consc_avg.df, by=c("Var1", "Var2")) %>%
      mutate(delta = ((avg.y - avg.x) / (avg.x + 1)) * 100)
    
    rd_delta <- rd_unconsc_avg.df %>%
      inner_join(rd_consc_avg.df, by=c("Var1", "Var2")) %>%
      mutate(delta = ((avg.y - avg.x) / (avg.x + 1)) * 100)
    
    fiber_delta <- fiber_unconsc_avg.df %>%
      inner_join(fiber_consc_avg.df, by=c("Var1", "Var2")) %>%
      mutate(delta = ((avg.y - avg.x) / (avg.x + 1)) * 100)
    
    out3 <- list(fa_delta=fa_delta, md_delta=md_delta, ad_delta=ad_delta, rd_delta=rd_delta, fiber_delta=fiber_delta)
    out3 <- lapply(out3, removeDiag, rm.region="upper")
    
    
  } else {
    fa_control_avg.df <- fa_control.df %>%
      group_by(Var1, Var2) %>%
      summarise(avg = mean(value, na.rm=T))
    
    md_control_avg.df <- md_control.df %>%
      group_by(Var1, Var2) %>%
      summarise(avg = mean(value, na.rm=T))
    
    ad_control_avg.df <- ad_control.df %>%
      group_by(Var1, Var2) %>%
      summarise(avg = mean(value, na.rm=T))
    
    rd_control_avg.df <- rd_control.df %>%
      group_by(Var1, Var2) %>%
      summarise(avg = mean(value, na.rm=T))
    
    fiber_control_avg.df <- fiber_control.df %>%
      group_by(Var1, Var2) %>%
      summarise(avg = mean(value, na.rm=T))
    
    # save and remove diagonals
    out <- list(fa_control=fa_control_avg.df, md_control=md_control_avg.df, ad_control=ad_control_avg.df, rd_control=rd_control_avg.df, fiber_control=fiber_control_avg.df)
    out <- lapply(out, removeDiag, rm.region="upper")
    
    
  }
  
  
  # save dataframe lists to master output list depending on whether controls were analyzed or not
  if (control == FALSE) {
    out_main <- list(conscious=out1, unconscious=out2, deltas=out3)
    return(out_main)
  } else {
    # can only output just the control data here. can't make deltas if there's no conscious state difference.
    out_main <- out
    return(out_main)
  }
  
  
}



PlotEachConnectivityMatrix <- function(datalist, lim_list, names_to_use) {
  out <- list()
  # choose the min/max of the plot based on whichever datatype is currently being iterated over
  for (name in names_to_use) {
    if (grepl("fiber", name, ignore.case=T)) {
      min_lim <- lim_list[["fiber"]][1]
      max_lim <- lim_list[["fiber"]][2]
    } else if (grepl("fa", name, ignore.case=T)) {
      min_lim <- lim_list[["fa"]][1]
      max_lim <- lim_list[["fa"]][2]
    } else if (grepl("md", name, ignore.case=T)) {
      min_lim <- lim_list[["md"]][1]
      max_lim <- lim_list[["md"]][2]
    } else if (grepl("ad", name, ignore.case=T)) {
      min_lim <- lim_list[["ad"]][1]
      max_lim <- lim_list[["ad"]][2]
    } else if (grepl("rd", name, ignore.case=T)) {
      min_lim <- lim_list[["rd"]][1]
      max_lim <- lim_list[["rd"]][2]
    } else {
      stop("Limit list doesn't contain the parameter name (fa, md, ad, rd, fiber)")
    }
    out[[name]] <- ggplot(datalist[[name]], aes(x=Var1, y=Var2, fill=value)) + 
      geom_tile(color="white") +
      scale_fill_gradient2(low="blue", high="red", mid="white",
                           midpoint = 0, limit=c(min_lim, max_lim), space="Lab", name="num fibers") +   #c(-460, 420)  c(-82, 1000)
      theme_minimal() +
      xlab("") + ylab("") +  # make axis labels blank
      theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "left") +
      scale_y_discrete(position = "right") +
      ggtitle(paste(name)) +
      coord_fixed()
    
  }
  return(out)
}



SavePlots <- function(plotlist, path) {
  for (i in 1:length(plotlist)) {
    fname <- gsub("\\.", "_", names(plotlist)[i])
    save_path <- paste(path, "/", fname, ".pdf", sep="")
    ggsave(save_path, plotlist[[i]], device="pdf")
  }
}



#### Load data ####
# path to save plots
# save plots
path <- "E:/Experiments/ICH_MRI/New Length Analysis"
path2 <- "E:/Experiments/ICH_MRI/New Length Analysis/Controls"
path3 <- "E:/Experiments/ICH_MRI/New Length Analysis/New Groups"
# set up files for reading
file_dir <- "E:/Experiments/ICH_MRI/Results_JanClaassen_ICH_DTI"
file_dir_ctrl <- "E:/Experiments/ICH_MRI/Results_Controls_ICH"

txtfiles <- list.files(path=file_dir, pattern=".txt", full.names=FALSE)
txtfiles_ctrl <- list.files(path=file_dir_ctrl, pattern=".txt", full.names=FALSE)

# datatypes to look at for connectivity matrices
datatypes <- c("FA", "MD", "AD", "RD", "fibercount")

# patient info for non-controls
conscdat <- read.xlsx("E:/Experiments/ICH_MRI/MRI DTI set.xlsx", 1)
conscdat$DTI.Number <- gsub(" ", "", gsub("DTI", "ICH", gsub("DTI ", "ICH", gsub("A.*", "", conscdat$DTI.Number))))

# lesion location info
locdat <- read.xlsx("E:/Experiments/ICH_MRI/lesion_location.xlsx")
# add in new grouping categories
# conscdat <- conscdat %>%
#   mutate(conscious.throughout = case_when(unconsious.on.admission==0 & unconsious.at.time.of.MRI==0 & unconscious.on.discharge==0 ~ 1,
#                                           TRUE ~ 0),
#          MRI.conscious.recover = case_when(unconsious.at.time.of.MRI==1 & unconscious.on.discharge==0 ~ 1,
#                                            TRUE ~ 0),
#          MRI.conscious.no.recover = case_when(unconsious.at.time.of.MRI==0 & unconscious.on.discharge==1 ~ 1,
#                                               TRUE ~ 0))


# should this instead be:
conscdat <- conscdat %>%
  mutate(conscious.throughout = case_when(unconsious.on.admission==0 & unconsious.at.time.of.MRI==0 & unconscious.on.discharge==0 ~ 1,
                                          TRUE ~ 0),
         MRI.conscious.recover = case_when((unconsious.on.admission==1 | unconsious.at.time.of.MRI==1) & unconscious.on.discharge==0 ~ 1,
                                           TRUE ~ 0),
         MRI.conscious.no.recover = case_when((unconsious.on.admission==1 | unconsious.at.time.of.MRI==1) & unconscious.on.discharge==1 ~ 1,
                                              TRUE ~ 0))


# make "summary" group column
conscdat <- conscdat %>%
  mutate(group = case_when(conscious.throughout==1 ~ "consc throughout",
                           MRI.conscious.recover==1 ~ "MRI consc recover",
                           MRI.conscious.no.recover==1 ~ "MRI consc no recover"))


# convert conscious columns to factor
conscdat[,c("unconsious.on.admission", "unconsious.at.time.of.MRI", "unconscious.on.discharge", "Recovery.of.consciousness", "conscious.throughout", "MRI.conscious.recover", "MRI.conscious.no.recover", "group")] <- lapply(
  conscdat[,c("unconsious.on.admission", "unconsious.at.time.of.MRI", "unconscious.on.discharge", "Recovery.of.consciousness", "conscious.throughout", "MRI.conscious.recover", "MRI.conscious.no.recover", "group")],
  as.factor
)

# conscious classification group names for plotting (from conscdat table)
consc_groups <- c("unconsious.on.admission", "unconsious.at.time.of.MRI", "unconscious.on.discharge", "Recovery.of.consciousness")
# New groups should be:
consc_groups2 <- c("conscious.throughout", "MRI.conscious.recover", "MRI.conscious.no.recover")
# the one summary group
consc_groups3 <- "group"

### ACTUALLY want to make an overall summary group column
### where each group has its own identifier number
### and THEN plot the data.


# fiber length names (as specified in the text file names)
lens <- c("L10_L20", "L10_L80", "L20_L40", "L40_L80", "L80_L160")






#### Analysis ####

# pull the data
parsedData <- ParseFiberLenData(txtfiles=txtfiles, lens=lens, file_dir=file_dir, control=FALSE, conscdat=conscdat)
parsedDataControls <- ParseFiberLenData(txtfiles=txtfiles_ctrl, lens=lens, file_dir=file_dir_ctrl, control=TRUE)

#### Combine relative data frames between both lists together ####
parsedDataAll <- mapply(rbind, parsedData, parsedDataControls, SIMPLIFY=F)

# pull data with new groups added. I can just use this as the "normal" variable now, since it's the same as the above one
# but now has the 3 new group columns added in.
#parsedData_newgroups <- ParseFiberLenData(txtfiles=txtfiles, lens=lens, file_dir=file_dir, control=FALSE, conscdat=conscdat)


# plots for each data type, saved into their respective list for controls and not controls.
# Each plot is an entry in the list
plots_alldata <- PlotEachDataType(datalist=parsedData, consc_groups=consc_groups3, control=FALSE)
plots_controls <- PlotEachDataType(datalist=parsedDataControls, consc_groups=consc_groups, control=TRUE)

# plots with the 3 conscious groups + controls
plots_All <- PlotEachDataType(datalist=parsedDataAll, consc_groups=consc_groups3)


# plotting using new groups
#plots_alldata_newgroups <- PlotEachDataType(datalist=parsedData_newgroups, consc_groups=consc_groups3, control=FALSE)


### plot showing patient IDs using the directlabels package
plots_All$faDat + directlabels::geom_dl(aes(label=ID), method=list(directlabels::dl.combine("first.points", "last.points"), cex=0.8))
plots_All$adDat + directlabels::geom_dl(aes(label=ID), method=list(directlabels::dl.combine("first.points", "last.points"), cex=0.8))






#### Save DTI plots ####

# SavePlots(plotlist=plots_alldata, path=path)
# SavePlots(plotlist=plots_controls, path=path2)
# # saving plots with new groups
# #SavePlots(plotlist=plots_alldata_newgroups, path=path3)
# 
# # save plots with the new 4 groups on same plot
# SavePlots(plotlist=plots_All, path=path3)









#### Connectivity Matrices ####

PatientConnectivity <- CreateConnectivityMatrix(txtfiles, file_dir, datatypes, control=FALSE, conscious_time="unconscious.on.discharge", conscdat=conscdat, locdat=locdat)


ControlConnectivity <- CreateConnectivityMatrix(txtfiles_ctrl, file_dir_ctrl, datatypes, control=TRUE, conscdat=NULL)




## test to see if get same results.
# yes, they are
# ggplot(PatientConnectivity$conscious$fiber_consc, aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile(color="white") +
#   scale_fill_gradient2(low="blue", high="red", mid="white",
#                        midpoint = 0, limit=c(0, 900), space="Lab", name="num fibers") +   #c(-460, 420)  c(-82, 1000)
#   theme_minimal() +
#   xlab("") + ylab("") +  # make axis labels blank
#   theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.position = "left") +
#   scale_y_discrete(position = "right") +
#   ggtitle("Fiber Count", subtitle="Conscious at MRI") +
#   coord_fixed()





#### Connectivity Matrix Plots ####

## Set parameters for connectivity matrix plots
# md, ad, and rd have the highest max value in the unconscious patients
lim_list <- list(fiber=c(0, 1700), fa=c(0, 0.5), md=c(0, 2), ad=c(0, 2.4), rd=c(0, 1.8))
consc_names_to_use <- names(PatientConnectivity$conscious) #"fiber_consc"
unconsc_names_to_use <- names(PatientConnectivity$unconscious) #"fiber_unconsc"
control_names_to_use <- names(ControlConnectivity) #"fiber_control"



## plots for controls
MatrixPlotsControls <- PlotEachConnectivityMatrix(datalist=ControlConnectivity, lim_list=lim_list, names_to_use=control_names_to_use)
MatrixPlotsControls$fa_control
MatrixPlotsControls$md_control
MatrixPlotsControls$ad_control
MatrixPlotsControls$rd_control
MatrixPlotsControls$fiber_control


## plots for consc patients
MatrixPlotsConsc <- PlotEachConnectivityMatrix(datalist=PatientConnectivity$conscious, lim_list=lim_list, names_to_use=consc_names_to_use)
MatrixPlotsConsc$fa_consc
MatrixPlotsConsc$md_consc
MatrixPlotsConsc$ad_consc
MatrixPlotsConsc$rd_consc
MatrixPlotsConsc$fiber_consc


## plots for unconscious patients
MatrixPlotsUnconsc <- PlotEachConnectivityMatrix(datalist=PatientConnectivity$unconscious, lim_list=lim_list, names_to_use=unconsc_names_to_use)
MatrixPlotsUnconsc$fa_unconsc
MatrixPlotsUnconsc$md_unconsc
MatrixPlotsUnconsc$ad_unconsc
MatrixPlotsUnconsc$rd_unconsc
MatrixPlotsUnconsc$fiber_unconsc







## For fiber count, looks like controls have same patterns, but they're stronger.
## number of fibers increases in this order: unconsc < consc < control





