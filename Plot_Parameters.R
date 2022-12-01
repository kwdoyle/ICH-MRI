library(openxlsx)
library(reshape2)
library(ggplot2)
library(dplyr)

### Analysis 1


setwd("E:/Experiments/ICH_MRI/Results_JanClaassen_ICH_DTI")


txtfiles <- list.files(pattern=".txt")
conscdat <- read.xlsx("../MRI DTI set.xlsx", 1)

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


split_txt <- strsplit(txtfiles, "_")

# extract "columns"
ID <- sapply(split_txt, "[", 1)
ID.u <- unique(ID)


# for each ID, find the files that correspond to it
# then can take these files and look for each length-type
# maybe set those types as parameters to search for
lens <- c("L10_L20", "L10_L80", "L20_L40", "L40_L80", "L80_L160")
# "L10_L160" is the entire length and the data from these are different,
# so don't include them in this part

# table to store data from all files in
mainTable <- data.frame()

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
      dat <- read.table(file, skip=6)  # there's also a blank line at the beginning of the file to skip
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

# want to make the length-types factors for plotting
mainTable$len <- as.factor(mainTable$len)

# join the consciousness data to mainTable
mainTable <- inner_join(mainTable, select(conscdat, DTI.Number, unconsious.on.admission, unconsious.at.time.of.MRI, unconscious.on.discharge, Recovery.of.consciousness),
           by=c("ID" = "DTI.Number"))



# filter data for 'fa'
faDat <- mainTable %>%
  filter(param == "fa" & len != 'L10_L80')

# filter data for 'md'
mdDat <- mainTable %>%
  filter(param == "md" & len != 'L10_L80')

# filter data for 'ad'
adDat <- mainTable %>%
  filter(param == "ad" & len != 'L10_L80')

# filter data for 'rd'
rdDat <- mainTable %>%
  filter(param == "rd" & len != 'L10_L80')

# filter data for 'txx'
txxDat <- mainTable %>%
  filter(param == "txx" & len != 'L10_L80')

# filter data for 'txy'
txyDat <- mainTable %>%
  filter(param == "txy" & len != 'L10_L80')

# filter data for 'txz'
txzDat <- mainTable %>%
  filter(param == "txz" & len != 'L10_L80')

# filter data for 'tyy'
tyyDat <- mainTable %>%
  filter(param == "tyy" & len != 'L10_L80')

# filter data for 'tyz'
tyzDat <- mainTable %>%
  filter(param == "tyz" & len != 'L10_L80')

# filter data for 'tzz'
tzzDat <- mainTable %>%
  filter(param == "tzz" & len != 'L10_L80')




# take ratio of shortest length to longest length
faDatRatio <- faDat %>%
  filter(len=="L10_L20" | len=="L80_L160") %>%
  group_by(ID) %>%
  mutate(ratio = mean[len=="L80_L160"] / mean[len=="L10_L20"])

# then do t test
t.test(ratio ~ unconscious.on.discharge, data=faDatRatio)
t.test(ratio ~ unconsious.at.time.of.MRI, data=faDatRatio)


  

### FA
p1 <- ggplot(faDat, aes(x=len, y=mean, group=ID, color=unconsious.at.time.of.MRI)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="FA - unconsious at time of MRI", x="length", y="mean", color="key")

p2 <- ggplot(faDat, aes(x=len, y=mean, group=ID, color=unconsious.on.admission)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="FA - unconsious on admission", x="length", y="mean", color="key")

p3 <- ggplot(faDat, aes(x=len, y=mean, group=ID, color=unconscious.on.discharge)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="FA - unconscious. on discharge", x="length", y="mean", color="key")

p4 <- ggplot(filter(faDat, !is.na(Recovery.of.consciousness)), aes(x=len, y=mean, group=ID, color=Recovery.of.consciousness)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="FA - Recovery of consciousness", x="length", y="mean", color="key")

ggsave("../lengths/FA unconscious at MRI.pdf", p1, device="pdf")
ggsave("../lengths/FA unconscious at Admission.pdf", p2, device="pdf")
ggsave("../lengths/FA unconscious at Dch.pdf", p3, device="pdf")
ggsave("../lengths/FA recovery of consciousness.pdf", p4, device="pdf")



### MD
p1 <- ggplot(mdDat, aes(x=len, y=mean, group=ID, color=unconsious.at.time.of.MRI)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="MD - unconsious at time of MRI", x="length", y="mean", color="key")

p2 <- ggplot(mdDat, aes(x=len, y=mean, group=ID, color=unconsious.on.admission)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="MD - unconsious on admission", x="length", y="mean", color="key")

p3 <- ggplot(mdDat, aes(x=len, y=mean, group=ID, color=unconscious.on.discharge)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="MD - unconscious on discharge", x="length", y="mean", color="key")

p4 <- ggplot(filter(mdDat, !is.na(Recovery.of.consciousness)), aes(x=len, y=mean, group=ID, color=Recovery.of.consciousness)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="MD - Recovery of consciousness", x="length", y="mean", color="key")

ggsave("../lengths/MD unconscious at MRI.pdf", p1, device="pdf")
ggsave("../lengths/MD unconscious at Admission.pdf", p2, device="pdf")
ggsave("../lengths/MD unconscious at Dch.pdf", p3, device="pdf")
ggsave("../lengths/MD recovery of consciousness.pdf", p4, device="pdf")



### AD
p1 <- ggplot(adDat, aes(x=len, y=mean, group=ID, color=unconsious.at.time.of.MRI)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="AD - unconsious at time of MRI", x="length", y="mean", color="key")

p2 <- ggplot(adDat, aes(x=len, y=mean, group=ID, color=unconsious.on.admission)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="AD - unconsious on admission", x="length", y="mean", color="key")

p3 <- ggplot(adDat, aes(x=len, y=mean, group=ID, color=unconscious.on.discharge)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="AD - unconscious on discharge", x="length", y="mean", color="key")

p4 <- ggplot(filter(adDat, !is.na(Recovery.of.consciousness)), aes(x=len, y=mean, group=ID, color=Recovery.of.consciousness)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="AD - Recovery of consciousness", x="length", y="mean", color="key")

ggsave("../lengths/AD unconscious at MRI.pdf", p1, device="pdf")
ggsave("../lengths/AD unconscious at Admission.pdf", p2, device="pdf")
ggsave("../lengths/AD unconscious at Dch.pdf", p3, device="pdf")
ggsave("../lengths/AD recovery of consciousness.pdf", p4, device="pdf")



### RD
p1 <- ggplot(rdDat, aes(x=len, y=mean, group=ID, color=unconsious.at.time.of.MRI)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="RD - unconsious at time of MRI", x="length", y="mean", color="key")

p2 <- ggplot(rdDat, aes(x=len, y=mean, group=ID, color=unconsious.on.admission)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="RD - unconsious on admission", x="length", y="mean", color="key")

p3 <- ggplot(rdDat, aes(x=len, y=mean, group=ID, color=unconscious.on.discharge)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="RD - unconscious on discharge", x="length", y="mean", color="key")

p4 <- ggplot(filter(rdDat, !is.na(Recovery.of.consciousness)), aes(x=len, y=mean, group=ID, color=Recovery.of.consciousness)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="RD - Recovery of consciousness", x="length", y="mean", color="key")

ggsave("../lengths/RD unconscious at MRI.pdf", p1, device="pdf")
ggsave("../lengths/RD unconscious at Admission.pdf", p2, device="pdf")
ggsave("../lengths/RD unconscious at Dch.pdf", p3, device="pdf")
ggsave("../lengths/RD recovery of consciousness.pdf", p4, device="pdf")



### TXX
p1 <- ggplot(txxDat, aes(x=len, y=mean, group=ID, color=unconsious.at.time.of.MRI)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TXX - unconsious at time of MRI", x="length", y="mean", color="key")

p2 <- ggplot(txxDat, aes(x=len, y=mean, group=ID, color=unconsious.on.admission)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TXX - unconsious on admission", x="length", y="mean", color="key")

p3 <- ggplot(txxDat, aes(x=len, y=mean, group=ID, color=unconscious.on.discharge)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TXX - unconscious on discharge", x="length", y="mean", color="key")

p4 <- ggplot(filter(txxDat, !is.na(Recovery.of.consciousness)), aes(x=len, y=mean, group=ID, color=Recovery.of.consciousness)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TXX - Recovery of consciousness", x="length", y="mean", color="key")

ggsave("../lengths/TXX unconscious at MRI.pdf", p1, device="pdf")
ggsave("../lengths/TXX unconscious at Admission.pdf", p2, device="pdf")
ggsave("../lengths/TXX unconscious at Dch.pdf", p3, device="pdf")
ggsave("../lengths/TXX recovery of consciousness.pdf", p4, device="pdf")



### TXY
p1 <- ggplot(txyDat, aes(x=len, y=mean, group=ID, color=unconsious.at.time.of.MRI)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TXY - unconsious at time of MRI", x="length", y="mean", color="key")

p2 <- ggplot(txyDat, aes(x=len, y=mean, group=ID, color=unconsious.on.admission)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TXY - unconsious on admission", x="length", y="mean", color="key")

p3 <- ggplot(txyDat, aes(x=len, y=mean, group=ID, color=unconscious.on.discharge)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TXY - unconscious on discharge", x="length", y="mean", color="key")

p4 <- ggplot(filter(txyDat, !is.na(Recovery.of.consciousness)), aes(x=len, y=mean, group=ID, color=Recovery.of.consciousness)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TXY - Recovery of consciousness", x="length", y="mean", color="key")

ggsave("../lengths/TXY unconscious at MRI.pdf", p1, device="pdf")
ggsave("../lengths/TXY unconscious at Admission.pdf", p2, device="pdf")
ggsave("../lengths/TXY unconscious at Dch.pdf", p3, device="pdf")
ggsave("../lengths/TXY recovery of consciousness.pdf", p4, device="pdf")



### TXZ
p1 <- ggplot(txzDat, aes(x=len, y=mean, group=ID, color=unconsious.at.time.of.MRI)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TXZ - unconsious at time of MRI", x="length", y="mean", color="key")

p2 <- ggplot(txzDat, aes(x=len, y=mean, group=ID, color=unconsious.on.admission)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TXZ - unconsious on admission", x="length", y="mean", color="key")

p3 <- ggplot(txzDat, aes(x=len, y=mean, group=ID, color=unconscious.on.discharge)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TXZ - unconscious on discharge", x="length", y="mean", color="key")

p4 <- ggplot(filter(txzDat, !is.na(Recovery.of.consciousness)), aes(x=len, y=mean, group=ID, color=Recovery.of.consciousness)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TXZ - Recovery of consciousness", x="length", y="mean", color="key")

ggsave("../lengths/TXZ unconscious at MRI.pdf", p1, device="pdf")
ggsave("../lengths/TXZ unconscious at Admission.pdf", p2, device="pdf")
ggsave("../lengths/TXZ unconscious at Dch.pdf", p3, device="pdf")
ggsave("../lengths/TXZ recovery of consciousness.pdf", p4, device="pdf")



### TYY
p1 <- ggplot(tyyDat, aes(x=len, y=mean, group=ID, color=unconsious.at.time.of.MRI)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TYY - unconsious at time of MRI", x="length", y="mean", color="key")

p2 <- ggplot(tyyDat, aes(x=len, y=mean, group=ID, color=unconsious.on.admission)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TYY - unconsious on admission", x="length", y="mean", color="key")

p3 <- ggplot(tyyDat, aes(x=len, y=mean, group=ID, color=unconscious.on.discharge)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TYY - unconscious on discharge", x="length", y="mean", color="key")

p4 <- ggplot(filter(tyyDat, !is.na(Recovery.of.consciousness)), aes(x=len, y=mean, group=ID, color=Recovery.of.consciousness)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TYY - Recovery of consciousness", x="length", y="mean", color="key")

ggsave("../lengths/TYY unconscious at MRI.pdf", p1, device="pdf")
ggsave("../lengths/TYY unconscious at Admission.pdf", p2, device="pdf")
ggsave("../lengths/TYY unconscious at Dch.pdf", p3, device="pdf")
ggsave("../lengths/TYY recovery of consciousness.pdf", p4, device="pdf")



### TYZ
p1 <- ggplot(tyzDat, aes(x=len, y=mean, group=ID, color=unconsious.at.time.of.MRI)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TYZ - unconsious at time of MRI", x="length", y="mean", color="key")

p2 <- ggplot(tyzDat, aes(x=len, y=mean, group=ID, color=unconsious.on.admission)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TYZ - unconsious on admission", x="length", y="mean", color="key")

p3 <- ggplot(tyzDat, aes(x=len, y=mean, group=ID, color=unconscious.on.discharge)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TYZ - unconscious on discharge", x="length", y="mean", color="key")

p4 <- ggplot(filter(tyzDat, !is.na(Recovery.of.consciousness)), aes(x=len, y=mean, group=ID, color=Recovery.of.consciousness)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TYZ - Recovery of consciousness", x="length", y="mean", color="key")

ggsave("../lengths/TYZ unconscious at MRI.pdf", p1, device="pdf")
ggsave("../lengths/TYZ unconscious at Admission.pdf", p2, device="pdf")
ggsave("../lengths/TYZ unconscious at Dch.pdf", p3, device="pdf")
ggsave("../lengths/TYZ recovery of consciousness.pdf", p4, device="pdf")



### TZZ
p1 <- ggplot(tzzDat, aes(x=len, y=mean, group=ID, color=unconsious.at.time.of.MRI)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TZZ - unconsious at time of MRI", x="length", y="mean", color="key")

p2 <- ggplot(tzzDat, aes(x=len, y=mean, group=ID, color=unconsious.on.admission)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TZZ - unconsious on admission", x="length", y="mean", color="key")

p3 <- ggplot(tzzDat, aes(x=len, y=mean, group=ID, color=unconscious.on.discharge)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TZZ - unconscious on discharge", x="length", y="mean", color="key")

p4 <- ggplot(filter(tzzDat, !is.na(Recovery.of.consciousness)), aes(x=len, y=mean, group=ID, color=Recovery.of.consciousness)) +
  geom_line(size=1.1) + theme_bw() +
  labs(title="TZZ - Recovery of consciousness", x="length", y="mean", color="key")

ggsave("../lengths/TZZ unconscious at MRI.pdf", p1, device="pdf")
ggsave("../lengths/TZZ unconscious at Admission.pdf", p2, device="pdf")
ggsave("../lengths/TZZ unconscious at Dch.pdf", p3, device="pdf")
ggsave("../lengths/TZZ recovery of consciousness.pdf", p4, device="pdf")

