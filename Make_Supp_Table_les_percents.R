makeLesionCountTable <- function(df, names_use, time="MRI", use_type="ICH", count_type="percent") {
  ## use_type can be either "ICH" or "edema"; time can be either "MRI" or "Dch
  # don't really need to specify names in advance..
  out <- data.frame(name=character(0), cIps=character(0), cCont=character(0), ucIps=character(0), ucCont=character(0),
                    stringsAsFactors = FALSE)
  
  names_use2 <- paste(names_use, use_type)
  
  for (i in 1:length(names_use)) {
    ix <- grepl(names_use[i], df$variable) #& grepl("ipsi", consc_loc_table$variable)
    row_use <- df[ix, ]
    if (time == "MRI") {
      val_cIps <- (row_use$cMRI_pcnt[row_use$type == use_type & (grepl("ipsi", row_use$variable) | grepl("_C", row_use$variable) )])
      #edemaval_cIps <- row_use$cMRI_pcnt[row_use$type == "edema" & grepl("ipsi", row_use$variable)]
      val_cCont <- (row_use$cMRI_pcnt[row_use$type == use_type & (grepl("contro", row_use$variable) | grepl("_C", row_use$variable) )])
      val_ucIps <- (row_use$ucMRI_pcnt[row_use$type == use_type & (grepl("ipsi", row_use$variable) | grepl("_C", row_use$variable) )])
      val_ucCont <- (row_use$ucMRI_pcnt[row_use$type == use_type & (grepl("contro", row_use$variable) | grepl("_C", row_use$variable) )])
      
      if (names_use[i] %in% c("MB", "AntPons", "Teg")) {
        val_cIps <- (row_use$cMRI_pcnt[row_use$type == use_type & grepl("_C", row_use$variable)])
        #edemaval_cIps <- row_use$cMRI_pcnt[row_use$type == "edema" & grepl("ipsi", row_use$variable)]
        val_cCont <- (row_use$cMRI_pcnt[row_use$type == use_type & grepl("_C", row_use$variable)])
        val_ucIps <- (row_use$ucMRI_pcnt[row_use$type == use_type &  grepl("_C", row_use$variable)])
        val_ucCont <- (row_use$ucMRI_pcnt[row_use$type == use_type & grepl("_C", row_use$variable)])
      }
      
      if (names_use[i] %in% c("Vermis", "IVH")) { #(length(val_cIps)==0 & length(val_cCont)==0 & length(val_ucIps)==0 & length(val_ucCont)==0) {
        val_cIps <- (row_use$cMRI_pcnt[row_use$type == use_type])
        #edemaval_cIps <- row_use$cMRI_pcnt[row_use$type == "edema" & grepl("ipsi", row_use$variable)]
        val_cCont <- (row_use$cMRI_pcnt[row_use$type == use_type])
        val_ucIps <- (row_use$ucMRI_pcnt[row_use$type == use_type])
        val_ucCont <- (row_use$ucMRI_pcnt[row_use$type == use_type])
      }
      
      if (count_type == "percent") {
        insert_row <- c(name=names_use2[i], cIps=val_cIps*100, cCont=val_cCont*100, ucIps=val_ucIps*100, ucCont=val_ucCont*100)
      } else if (count_type == "count") {
        insert_row <- c(name=names_use2[i], cIps=val_cIps, cCont=val_cCont, ucIps=val_ucIps, ucCont=val_ucCont)
      }
      
      #out <- rbind(out, insert_row)
      out[nrow(out)+1, ] <- insert_row
      
    } else if (time == "Dch") {
      val_cIps <- (row_use$cDC_pcnt[row_use$type == use_type & (grepl("ipsi", row_use$variable) | grepl("_C", row_use$variable) )])
      #edemaval_cIps <- row_use$cMRI_pcnt[row_use$type == "edema" & grepl("ipsi", row_use$variable)]
      val_cCont <- (row_use$cDC_pcnt[row_use$type == use_type & (grepl("contro", row_use$variable) | grepl("_C", row_use$variable) )])
      val_ucIps <- (row_use$ucDC_pcnt[row_use$type == use_type & (grepl("ipsi", row_use$variable) | grepl("_C", row_use$variable) )])
      val_ucCont <- (row_use$ucDC_pcnt[row_use$type == use_type & (grepl("contro", row_use$variable) | grepl("_C", row_use$variable) )])
      
      if (length(val_cIps)==0 & length(val_cCont)==0 & length(val_ucIps)==0 & length(val_ucCont)==0) {
        val_cIps <- (row_use$cDC_pcnt[row_use$type == use_type])
        #edemaval_cIps <- row_use$cMRI_pcnt[row_use$type == "edema" & grepl("ipsi", row_use$variable)]
        val_cCont <- (row_use$cDC_pcnt[row_use$type == use_type])
        val_ucIps <- (row_use$ucDC_pcnt[row_use$type == use_type])
        val_ucCont <- (row_use$ucDC_pcnt[row_use$type == use_type])
      }
      
      if (count_type == "percent") {
        insert_row <- c(name=names_use2[i], cIps=val_cIps*100, cCont=val_cCont*100, ucIps=val_ucIps*100, ucCont=val_ucCont*100)
      } else if (count_type == "count") {
        insert_row <- c(name=names_use2[i], cIps=val_cIps, cCont=val_cCont, ucIps=val_ucIps, ucCont=val_ucCont)
      }
      
      #out <- rbind(out, insert_row)
      out[nrow(out)+1, ] <- insert_row
      
    }
    
  }
  return(out)
}











# old versions
makeLesionPcntTable_oldv1 <- function(df, time="MRI", use_type="ICH") {
  ## use_type can be either "ICH" or "edema"; time can be either "MRI" or "Dch
  # don't really need to specify names in advance..
  out <- data.frame(name=character(0), cIps=character(0), cCont=character(0), ucIps=character(0), ucCont=character(0),
                    stringsAsFactors = FALSE)
  
  names_use <- c("FCx", "PCx", "TCx", "OCx", "INS", "IC_ant", "IC_post", "Caudate", "PUT", "GP", "TH", "BF", "Hypo", "MB_peduncle", "MB", "AntPons", "Teg", "Cereb", "Vermis", "IVH")
  names_use2 <- paste(names_use, use_type)
  
  for (i in 1:length(names_use)) {
    ix <- grepl(names_use[i], df$variable) #& grepl("ipsi", consc_loc_table$variable)
    row_use <- df[ix, ]
    if (time == "MRI") {
      val_cIps <- (row_use$cMRI_pcnt[row_use$type == use_type & (grepl("ipsi", row_use$variable) | grepl("_C", row_use$variable) )]) * 100
      #edemaval_cIps <- row_use$cMRI_pcnt[row_use$type == "edema" & grepl("ipsi", row_use$variable)]
      val_cCont <- (row_use$cMRI_pcnt[row_use$type == use_type & (grepl("contro", row_use$variable) | grepl("_C", row_use$variable) )]) * 100
      val_ucIps <- (row_use$ucMRI_pcnt[row_use$type == use_type & (grepl("ipsi", row_use$variable) | grepl("_C", row_use$variable) )]) * 100
      val_ucCont <- (row_use$ucMRI_pcnt[row_use$type == use_type & (grepl("contro", row_use$variable) | grepl("_C", row_use$variable) )]) * 100
      
      if (length(val_cIps)==0 & length(val_cCont)==0 & length(val_ucIps)==0 & length(val_ucCont)==0) {
        val_cIps <- (row_use$cMRI_pcnt[row_use$type == use_type]) * 100
        #edemaval_cIps <- row_use$cMRI_pcnt[row_use$type == "edema" & grepl("ipsi", row_use$variable)]
        val_cCont <- (row_use$cMRI_pcnt[row_use$type == use_type]) * 100
        val_ucIps <- (row_use$ucMRI_pcnt[row_use$type == use_type]) * 100
        val_ucCont <- (row_use$ucMRI_pcnt[row_use$type == use_type]) * 100
      }
      
      insert_row <- c(name=names_use2[i], cIps=val_cIps, cCont=val_cCont, ucIps=val_ucIps, ucCont=val_ucCont)
      #out <- rbind(out, insert_row)
      out[nrow(out)+1, ] <- insert_row
      
    } else if (time == "Dch") {
      val_cIps <- (row_use$cDC_pcnt[row_use$type == use_type & (grepl("ipsi", row_use$variable) | grepl("_C", row_use$variable) )]) * 100
      #edemaval_cIps <- row_use$cMRI_pcnt[row_use$type == "edema" & grepl("ipsi", row_use$variable)]
      val_cCont <- (row_use$cDC_pcnt[row_use$type == use_type & (grepl("contro", row_use$variable) | grepl("_C", row_use$variable) )]) * 100
      val_ucIps <- (row_use$ucDC_pcnt[row_use$type == use_type & (grepl("ipsi", row_use$variable) | grepl("_C", row_use$variable) )]) * 100
      val_ucCont <- (row_use$ucDC_pcnt[row_use$type == use_type & (grepl("contro", row_use$variable) | grepl("_C", row_use$variable) )]) * 100
      
      if (length(val_cIps)==0 & length(val_cCont)==0 & length(val_ucIps)==0 & length(val_ucCont)==0) {
        val_cIps <- (row_use$cDC_pcnt[row_use$type == use_type]) * 100
        #edemaval_cIps <- row_use$cMRI_pcnt[row_use$type == "edema" & grepl("ipsi", row_use$variable)]
        val_cCont <- (row_use$cDC_pcnt[row_use$type == use_type]) * 100
        val_ucIps <- (row_use$ucDC_pcnt[row_use$type == use_type]) * 100
        val_ucCont <- (row_use$ucDC_pcnt[row_use$type == use_type]) * 100
      }
      
      insert_row <- c(name=names_use2[i], cIps=val_cIps, cCont=val_cCont, ucIps=val_ucIps, ucCont=val_ucCont)
      #out <- rbind(out, insert_row)
      out[nrow(out)+1, ] <- insert_row

    }
    
  }
  return(out)
}




makeLesionPcntTable_old <- function(df, time="MRI", use_type="ICH") {
  ## use_type can be either "ICH" or "edema"; time can be either "MRI" or "Dch
  # don't really need to specify names in advance..
  out <- data.frame(name=character(0), cIps=character(0), cCont=character(0), ucIps=character(0), ucCont=character(0),
                    stringsAsFactors = FALSE)
  
  names_use <- c("FCx", "PCx", "TCx", "OCx", "INS", "IC_ant", "IC_post", "Caudate", "PUT", "GP", "TH", "BF", "Hypo", "MB_peduncle", "MB", "AntPons", "Teg", "Cereb", "Vermis", "IVH")
  names_use2 <- paste(names_use, use_type)
  
  for (i in 1:length(names_use)) {
    ix <- grepl(names_use[i], df$variable) #& grepl("ipsi", consc_loc_table$variable)
    row_use <- df[ix, ]
    if (time == "MRI") {
      val_cIps <- (row_use$cMRI_pcnt[row_use$type == use_type & (grepl("ipsi", row_use$variable) | grepl("_C", row_use$variable) )]) * 100
      #edemaval_cIps <- row_use$cMRI_pcnt[row_use$type == "edema" & grepl("ipsi", row_use$variable)]
      val_cCont <- (row_use$cMRI_pcnt[row_use$type == use_type & (grepl("contro", row_use$variable) | grepl("_C", row_use$variable) )]) * 100
      val_ucIps <- (row_use$ucMRI_pcnt[row_use$type == use_type & (grepl("ipsi", row_use$variable) | grepl("_C", row_use$variable) )]) * 100
      val_ucCont <- (row_use$ucMRI_pcnt[row_use$type == use_type & (grepl("contro", row_use$variable) | grepl("_C", row_use$variable) )]) * 100
      
      if (names_use[i] %in% c("MB", "AntPons", "Teg")) {
        val_cIps <- (row_use$cMRI_pcnt[row_use$type == use_type & grepl("_C", row_use$variable)]) * 100
        #edemaval_cIps <- row_use$cMRI_pcnt[row_use$type == "edema" & grepl("ipsi", row_use$variable)]
        val_cCont <- (row_use$cMRI_pcnt[row_use$type == use_type & grepl("_C", row_use$variable)]) * 100
        val_ucIps <- (row_use$ucMRI_pcnt[row_use$type == use_type &  grepl("_C", row_use$variable)]) * 100
        val_ucCont <- (row_use$ucMRI_pcnt[row_use$type == use_type & grepl("_C", row_use$variable)]) * 100
      }
      
      if (names_use[i] %in% c("Vermis", "IVH")) { #(length(val_cIps)==0 & length(val_cCont)==0 & length(val_ucIps)==0 & length(val_ucCont)==0) {
        val_cIps <- (row_use$cMRI_pcnt[row_use$type == use_type]) * 100
        #edemaval_cIps <- row_use$cMRI_pcnt[row_use$type == "edema" & grepl("ipsi", row_use$variable)]
        val_cCont <- (row_use$cMRI_pcnt[row_use$type == use_type]) * 100
        val_ucIps <- (row_use$ucMRI_pcnt[row_use$type == use_type]) * 100
        val_ucCont <- (row_use$ucMRI_pcnt[row_use$type == use_type]) * 100
      }
      
      insert_row <- c(name=names_use2[i], cIps=val_cIps, cCont=val_cCont, ucIps=val_ucIps, ucCont=val_ucCont)
      #out <- rbind(out, insert_row)
      out[nrow(out)+1, ] <- insert_row
      
    } else if (time == "Dch") {
      val_cIps <- (row_use$cDC_pcnt[row_use$type == use_type & (grepl("ipsi", row_use$variable) | grepl("_C", row_use$variable) )]) * 100
      #edemaval_cIps <- row_use$cMRI_pcnt[row_use$type == "edema" & grepl("ipsi", row_use$variable)]
      val_cCont <- (row_use$cDC_pcnt[row_use$type == use_type & (grepl("contro", row_use$variable) | grepl("_C", row_use$variable) )]) * 100
      val_ucIps <- (row_use$ucDC_pcnt[row_use$type == use_type & (grepl("ipsi", row_use$variable) | grepl("_C", row_use$variable) )]) * 100
      val_ucCont <- (row_use$ucDC_pcnt[row_use$type == use_type & (grepl("contro", row_use$variable) | grepl("_C", row_use$variable) )]) * 100
      
      if (length(val_cIps)==0 & length(val_cCont)==0 & length(val_ucIps)==0 & length(val_ucCont)==0) {
        val_cIps <- (row_use$cDC_pcnt[row_use$type == use_type]) * 100
        #edemaval_cIps <- row_use$cMRI_pcnt[row_use$type == "edema" & grepl("ipsi", row_use$variable)]
        val_cCont <- (row_use$cDC_pcnt[row_use$type == use_type]) * 100
        val_ucIps <- (row_use$ucDC_pcnt[row_use$type == use_type]) * 100
        val_ucCont <- (row_use$ucDC_pcnt[row_use$type == use_type]) * 100
      }
      
      insert_row <- c(name=names_use2[i], cIps=val_cIps, cCont=val_cCont, ucIps=val_ucIps, ucCont=val_ucCont)
      #out <- rbind(out, insert_row)
      out[nrow(out)+1, ] <- insert_row
      
    }
    
  }
  return(out)
}




makeLesionCountTable_old <- function(df, time="MRI", use_type="ICH") {
  ## use_type can be either "ICH" or "edema"; time can be either "MRI" or "Dch
  # don't really need to specify names in advance..
  out <- data.frame(name=character(0), cIps=character(0), cCont=character(0), ucIps=character(0), ucCont=character(0),
                    stringsAsFactors = FALSE)
  
  names_use <- c("FCx", "PCx", "TCx", "OCx", "INS", "IC_ant", "IC_post", "Caudate", "PUT", "GP", "TH", "BF", "Hypo", "MB_peduncle", "MB", "AntPons", "Teg", "Cereb", "Vermis", "IVH")
  names_use2 <- paste(names_use, use_type)
  
  for (i in 1:length(names_use)) {
    ix <- grepl(names_use[i], df$variable) #& grepl("ipsi", consc_loc_table$variable)
    row_use <- df[ix, ]
    if (time == "MRI") {
      val_cIps <- (row_use$cMRI_pcnt[row_use$type == use_type & (grepl("ipsi", row_use$variable) | grepl("_C", row_use$variable) )])
      #edemaval_cIps <- row_use$cMRI_pcnt[row_use$type == "edema" & grepl("ipsi", row_use$variable)]
      val_cCont <- (row_use$cMRI_pcnt[row_use$type == use_type & (grepl("contro", row_use$variable) | grepl("_C", row_use$variable) )])
      val_ucIps <- (row_use$ucMRI_pcnt[row_use$type == use_type & (grepl("ipsi", row_use$variable) | grepl("_C", row_use$variable) )])
      val_ucCont <- (row_use$ucMRI_pcnt[row_use$type == use_type & (grepl("contro", row_use$variable) | grepl("_C", row_use$variable) )])
      
      if (length(val_cIps)==0 & length(val_cCont)==0 & length(val_ucIps)==0 & length(val_ucCont)==0) {
        val_cIps <- (row_use$cMRI_pcnt[row_use$type == use_type])
        #edemaval_cIps <- row_use$cMRI_pcnt[row_use$type == "edema" & grepl("ipsi", row_use$variable)]
        val_cCont <- (row_use$cMRI_pcnt[row_use$type == use_type])
        val_ucIps <- (row_use$ucMRI_pcnt[row_use$type == use_type])
        val_ucCont <- (row_use$ucMRI_pcnt[row_use$type == use_type])
      }
      
      insert_row <- c(name=names_use2[i], cIps=val_cIps, cCont=val_cCont, ucIps=val_ucIps, ucCont=val_ucCont)
      #out <- rbind(out, insert_row)
      out[nrow(out)+1, ] <- insert_row
      
    } else if (time == "Dch") {
      val_cIps <- (row_use$cDC_pcnt[row_use$type == use_type & (grepl("ipsi", row_use$variable) | grepl("_C", row_use$variable) )])
      #edemaval_cIps <- row_use$cMRI_pcnt[row_use$type == "edema" & grepl("ipsi", row_use$variable)]
      val_cCont <- (row_use$cDC_pcnt[row_use$type == use_type & (grepl("contro", row_use$variable) | grepl("_C", row_use$variable) )])
      val_ucIps <- (row_use$ucDC_pcnt[row_use$type == use_type & (grepl("ipsi", row_use$variable) | grepl("_C", row_use$variable) )])
      val_ucCont <- (row_use$ucDC_pcnt[row_use$type == use_type & (grepl("contro", row_use$variable) | grepl("_C", row_use$variable) )])
      
      if (length(val_cIps)==0 & length(val_cCont)==0 & length(val_ucIps)==0 & length(val_ucCont)==0) {
        val_cIps <- (row_use$cDC_pcnt[row_use$type == use_type])
        #edemaval_cIps <- row_use$cMRI_pcnt[row_use$type == "edema" & grepl("ipsi", row_use$variable)]
        val_cCont <- (row_use$cDC_pcnt[row_use$type == use_type])
        val_ucIps <- (row_use$ucDC_pcnt[row_use$type == use_type])
        val_ucCont <- (row_use$ucDC_pcnt[row_use$type == use_type])
      }
      
      insert_row <- c(name=names_use2[i], cIps=val_cIps, cCont=val_cCont, ucIps=val_ucIps, ucCont=val_ucCont)
      #out <- rbind(out, insert_row)
      out[nrow(out)+1, ] <- insert_row
      
    }
    
  }
  return(out)
}