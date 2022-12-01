
# set empty multi-dimentional matrices for each datatype, with a dimention for each patient?
# test <- array(rep(NA,900) , c(30, 30, 3))  # 3 should be the number of patients?
# # can extract/assign to a given dimension for a patient with:
# test[,,1] <- datmat.n
# 
# 
# ## OR NOT. just go through each ID, then each variabletype,
# ## then generate the matrix
# ## want to bind with rev.along=1 to bind on last dimension (which holds the matrices for each datatype)
# test <- abind(test, datmat.n, rev.along=1)
# 
# # can build an empty array with the correct number of dimensions
# # set the 3rd D to 0. it's technically index-able, but there's no values in it.
# test <- array(data=NA, dim=c(30,30,0))
# test <- abind(test, datmat.n, rev.along=1)




## make real empty arrays
# fa_array <- array(data=NA, dim=c(30,30,0))
# md_array <- array(data=NA, dim=c(30,30,0))
# ad_array <- array(data=NA, dim=c(30,30,0))
# rd_array <- array(data=NA, dim=c(30,30,0))
# fiber_array <- array(data=NA, dim=c(30,30,0))



### Old Heatmaps
# heatmap(fiber_unconsc_avg, scale="none", Rowv=NA, Colv=NA, revC=F,
#         col=rev(heat.colors(15)),
#         main="Fiber Count - Unconscious at MRI",
#         # this max range part just chooses whichever max from the unconsc or consc is higher.
#         breaks=seq(0.0001, 895,#max(max(fiber_unconsc_avg, na.rm=T), max(fiber_consc_avg, na.rm=T)),
#                    length.out=length(rev(heat.colors(15)))+1))
# 
# heatmap(fiber_consc_avg, scale="none", Rowv=NA, Colv=NA, revC=F,
#         col=rev(heat.colors(15)),
#         main="Fiber Count - Conscious at MRI",
#         breaks=seq(0.0001, 895,#max(max(fiber_unconsc_avg, na.rm=T), max(fiber_consc_avg, na.rm=T)),
#                    length.out=length(rev(heat.colors(15)))+1))
# 
# 
# 
# heatmap(fa_delta, scale="none", Rowv=NA, Colv=NA, revC=F,
#         col=rev(heat.colors(15)),
#         main="Fiber Count Delta - at MRI",
#         breaks=seq(-0.4, 0.3,#max(max(fiber_unconsc_avg, na.rm=T), max(fiber_consc_avg, na.rm=T)),
#                    length.out=length(rev(rainbow(15)))+1))




### What if I melt the "correlation matrix" and then can plot with ggplot?
#fiber_consc_l <- fiber_consc_avg[upper.tri(fiber_consc_avg)] <- NA

# fiber_consc_m <- melt(fiber_consc_avg, na.rm=T)
# fiber_consc_m <- melt(fiber_consc_avg, na.rm=T)
# 
# # why are the NAs here all dark grey?
# # it was because I forgot to add in "na.rm=T" to the melt function.
# ggplot(fiber_consc_m, aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile(color="white") +
#   scale_fill_gradient2(low="blue", high="red", mid="white",
#                        midpoint = 0, limit=c(0, 895), space="Lab", name="value") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle=90, vjust=1, hjust=1)) +
#   coord_fixed()

## Ok, the above plot works now.

## Try to plot the deltas w/ blue negative and red positive









# ggplot(fiber_delta_m1, aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile(color="white") +
#   scale_fill_gradient2(low="blue", high="red", mid="white",
#                        midpoint = 0, limit=c(min(fiber_delta1), max(fiber_delta1)), space="Lab", name="delta fiber num") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle=90, vjust=1, hjust=1)) +
#   ggtitle("fiber delta: unconsc - consc") +
#   coord_fixed()













## Want to use this one for delta
# ggplot(fiber_delta_m, aes(x=Var1, y=Var2, fill=value)) + 
#   geom_tile(color="white") +
#   scale_fill_gradient2(low="blue", high="red", mid="white",
#                        midpoint = 0, limit=c(-460, 420), space="Lab", name="delta fiber num") +   #c(-460, 420)  c(-82, 1000)
#   theme_minimal() +
#   xlab("") + ylab("") +  # make axis labels blank
#   theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
#   ggtitle("fiber delta: consc - unconsc", subtitle="Conscious at MRI") +
#   coord_fixed()


#limit=c(min(fiber_delta2, na.rm=T), max(fiber_delta2, na.rm=T))





###







## previous analysis going through all files.
# for (i in 1:length(analysis2)) {
#   
# 
#   dat <- read.table(analysis2[i])[-1,-1]
#   #title <- strsplit(analysis2[1], "_")[[1]][1]
#   title <- analysis2[i]
#   
#   dat2 <- dat[-1,-1]
#   # very complicated way to get the names first row, turn from data.frame to a vector, then turn only these names to factor levels (because this keeps the levels for all the different data types in each column)
#   colnames(dat2) <- as.factor(as.character(unname(unlist(dat[1,][-1]))))
#   rownames(dat2) <- dat[,1][-1]
#   
#   # turn all columns to numeric
#   #apply(dat2, as.numeric)
#   
#   
#   datmat <- as.matrix(dat2)
#   
#   # turn to numeric
#   datmat.n <- apply(datmat, c(1,2), as.numeric)
#   # put back row
#   
#   ### NOTE:
#   ### Can remove the lower/upper "triangle" of the diagonal matrix
#   ### using this function: 
#   ### lower.tri(datmat.n)   or  upper.tri(datmat.n)
#   ### and then just setting the values indexed there to NA:
#   ### datmat.n[lower.tri(datmat.n)] <- NA
#   
#   datmat.n[lower.tri(datmat.n)] <- NA
#   
#   
#   # save image as pdf
#   pdf(file=paste("../heatmaps/", title, ".pdf", sep=""))
#   
#   # Can make the diagonal go the other way if set revC to FALSE
#   heatmap(datmat.n, scale="none", Rowv=NA, Colv=NA, revC=T, 
#           col=rev(heat.colors(15)),
#           main=title,
#           # The range of this seq for the breaks refers to the values in the matrix
#           # that this will set the color-gradient too.
#           # e.g., 0 should be the smallest value, since there can't be a value
#           # smaller than 0 in these matrices.
#           # ...but setting it to a very small number above 0 allows for
#           # the actual 0-diagonal to appear completely white.
#           breaks=seq(0.0001, 3, length.out=length(rev(heat.colors(15)))+1))
#   
#   
#   dev.off()
# 
# 
# 
# }
# 



