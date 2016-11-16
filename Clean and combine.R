setwd("data")

##############################################################        ##############################################################
##########               LA                   ################        ##########               LA                   ################ 
##############################################################        ##############################################################

Lab="LA" # Enter name of fatty acid here
Y=read.csv(file= "linoleic.csv",header=T)
rownames(Y)=Y$Name


Y=Y[grepl(pattern="blank",Y$Type,ignore.case = T)==F,] ## removes blanks from dataframe
LOD=rep(T,time=nrow(Y)); LOQ=rep(T,time=nrow(Y)) ## creates variables for LOD and LOQ

####  Standards  ####
stnd=Y[grepl(pattern="standard",Y$Type,ignore.case = T)==T,] ## making a dataframe of the standards
stnd=stnd[grepl(pattern="bbX",stnd$Primary.Flags,ignore.case = F)==F,] ## removing excluded standards from dataframe 
stnd[is.na(stnd)]<-0
max.stnd=max(stnd$Conc.) # creates a value for the maximum standard curve value
max.val=rep(T,time=nrow(Y)); Conc=numeric()

# This loop will create: 
for(i in 1:nrow(Y)) { 
  Y$LOD[i]= ifelse(Y$S.N[i]<3, F,T) # limit of detection (LOD) flag
  Y$LOQ[i]= ifelse(Y$S.N[i]<10, F, T) #limit of quantification (LOQ) flag
  Conc[i] = ifelse(Y$S.N[i] < 10, NA, Y$Conc.[i]) #removing concentration values if below LOQ
}

Y$Conc=Conc   #adding new concentration variable to dataframe

### determining values above the max standard curve ###
Y$Above.curve=max.val ## adds a vector with T/F for if concentration is above max.stnd value
table(Y$Conc.>max.stnd) ## number of samples above the max.stnd
which(Y$Conc.>max.stnd)   # which sample rows are above the max.stnd
plot(x=Y$RT,y= Y$Conc.,xlab = 'Retention Time',ylab = 'Concentration',main = Lab); abline(h=max.stnd, lty=2); points(y=stnd$Conc.,x=stnd$RT ,col='red', pch=13)


###### Normalizing Concentrations #######
Y$Conc= Y$Conc * dilfac * FA

setwd("..");setwd("objects")
######  Saving objects #####
LA=Y #Enter the name of fatty acid here for naming dataframe
save(LA, file = 'LA.rda') 

#### splitting off the analytes
ana=Y[grepl(pattern="analyte",Y$Type,ignore.case = T)==T,] ## making a dataframe of the analytes
rownames(ana)=ana$Name

# Compiling fatty acid concentration vectors into a single data sheet
Ox.dat=data.frame("File"=ana$Name,"LA"=ana$Conc) # Paste the fatty acid object name here.
rownames(Ox.dat)=Ox.dat$File
try(if((all.equal(rownames(Ox.dat),as.character(ana$Name), as.character(Ox.dat$File)))==F) stop(print(rep("Error", time=200)))) #If names do not match, will report back with "Error"

tmpa=ana;tmpa[is.na(tmpa)]<-0 #converting NAs to 0, so a table can be generated
tmp=table(tmpa$Conc.>max.stnd)
tmp=ifelse(tmp[1]==nrow(ana), NA, nrow(ana)-tmp[1])
AboveMax=data.frame("LA"=tmp[1]) #Reports the number of samples with values above the highest standard curve value
rownames(AboveMax)="Number Above Max Curve Concentration"
tmp2=table(ana$S.N.LOQ)
tmp2=ifelse(tmp2==nrow(ana), NA, tmp2)
BelowLOQ=data.frame("LA"=tmp2[1]) #Reports the number of samples with values above the LOQ (sn>10)
rownames(BelowLOQ)="Number Below LOQ"


######  Archive the min and max detectable Concentration  ######

min.val=min(tmpa$Conc[tmpa$Conc>0])
Min.val=data.frame("LA"= min.val) #index the min value in your samples
row.names(Min.val)="Min Detectable Value"
Max.val=data.frame("LA"=max.stnd) #index the highest standard curve value 
rownames(Max.val)="Max Curve Value"


save(AboveMax,BelowLOQ, file = "MaxLOQ.rda")
save(Min.val,Max.val, file="MinMax.rda")
save(Ox.dat, file = "Ox.dat.rda")


rm(Y, i, LOD, LOQ, max.val,min.val,Lab, tmpa, tmp, tmp2,Conc, AboveMax,ana,BelowLOQ,Max.val,Ox.dat,stnd, LA, max.stnd, Min.val)


setwd("..");setwd("data")
##############################################################        ##############################################################
##########               13_hode                   ##########        ##########               13_hode                   ########## 
##############################################################        ##############################################################

Lab="13_hode" # Enter name of fatty acid here
Y=read.csv(file= "~/dropbox/Austin/Research/Targeted Oxylipid Analysis/2016 data extracion/13_hode.csv",header=T)
rownames(Y)=Y$Name

Y=Y[grepl(pattern="blank",Y$Type,ignore.case = T)==F,] ## removes blanks from dataframe
LOD=rep(T,time=nrow(Y)); LOQ=rep(T,time=nrow(Y)) ## creates variables for LOD and LOQ

####  Standards  ####
stnd=Y[grepl(pattern="standard",Y$Type,ignore.case = T)==T,] ## making a dataframe of the standards
stnd=stnd[grepl(pattern="bbX",stnd$Primary.Flags,ignore.case = F)==F,] ## removing excluded standards from dataframe 
stnd[is.na(stnd)]<-0
max.stnd=max(stnd$Conc.) # creates a value for the maximum standard curve value
max.val=rep(T,time=nrow(Y)); Conc=numeric()

# This loop will create: 
for(i in 1:nrow(Y)) { 
  Y$LOD[i]= ifelse(Y$S.N[i]<3, F,T) # limit of detection (LOD) flag
  Y$LOQ[i]= ifelse(Y$S.N[i]<10, F, T) #limit of quantification (LOQ) flag
  Conc[i] = ifelse(Y$S.N[i] < 10, NA, Y$Conc.[i]) #removing concentration values if below LOQ
}

Y$Conc=Conc   #adding new concentration variable to dataframe

### determining values above the max standard curve ###
Y$Above.curve=max.val ## adds a vector with T/F for if concentration is above max.stnd value
table(Y$Conc.>max.stnd) ## number of samples above the max.stnd
which(Y$Conc.>max.stnd)   # which sample rows are above the max.stnd
plot(x=Y$RT,y= Y$Conc.,xlab = 'Retention Time',ylab = 'Concentration',main = Lab); abline(h=max.stnd, lty=2); points(y=stnd$Conc.,x=stnd$RT ,col='red', pch=13)


###### Normalizing Concentrations ########
Y$Conc= Y$Conc * dilfac * Hodes

setwd("..");setwd("objects")
######  Saving objects #####
X13_hode=Y #Enter the name of fatty acid here for naming dataframe
save(X13_hode, file = '13_hode.rda') 

#### splitting off the analytes
ana=Y[grepl(pattern="analyte",Y$Type,ignore.case = T)==T,] ## making a dataframe of the analytes
rownames(ana)=ana$Name

# Compiling fatty acid concentration vectors into a single data sheet
load("Ox.dat.rda")
try(if((all.equal.character(rownames(Ox.dat), rownames(ana)))==F) 
  stop(print(rep("Error", time=200)))) #If names do not match, will report back with "Error"
Ox.dat=data.frame(Ox.dat,"13_hode"=ana$Conc) # Paste the fatty acid object name here.

####  Indexing number above max standard and number below limit of detection  ####
load("MaxLOQ.rda")
tmpa=ana;tmpa[is.na(tmpa)]<-0 #converting NAs to 0, so a table can be generated
tmp=table(tmpa$Conc.>max.stnd)
tmp=ifelse(tmp[1]==nrow(ana), NA, nrow(ana)-tmp[1])
AboveMax=data.frame(AboveMax,"13_hode"=tmp[1]) #Reports the number of samples with values above the highest standard curve value
tmp2=table(ana$S.N.LOQ)
tmp2=ifelse(tmp2==nrow(ana), NA, tmp2)
BelowLOQ=data.frame(BelowLOQ,"13_hode"=tmp2[1]) #Reports the number of samples with values above the LOQ (sn>10)


######  Archive the min and max detectable concentration  ######
load("MinMax.rda")
min.val=min(tmpa$Conc[tmpa$Conc>0])
Min.val=data.frame(Min.val,"X13_hode"= min.val) #index the min value in your samples
Max.val=data.frame(Max.val,"X13_hode"=max.stnd) #index the highest standard curve value 


save(AboveMax,BelowLOQ, file = "MaxLOQ.rda")
save(Min.val,Max.val, file="MinMax.rda")
save(Ox.dat, file = "Ox.dat.rda")


rm(Y, i, LOD, LOQ, max.val,min.val,Lab, tmpa, tmp, tmp2,Conc, AboveMax,ana,BelowLOQ,Max.val,Ox.dat,stnd, X13_hode, max.stnd, Min.val)

