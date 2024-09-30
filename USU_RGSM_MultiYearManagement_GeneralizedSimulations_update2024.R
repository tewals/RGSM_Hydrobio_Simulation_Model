#==========================================================
# Model for simulating the performance of alternative
#   multi-year water management strategies in terms of 
#   expected RGSM densities across a range of water-year
#   types (wet/dry)
#
# Coded by: Timothy Walsworth
# e-mail: timothy.walsworth@usu.edu
#==========================================================

###########################################
# This code will simulate the predicted
#   relative performance of different
#   broadly-defined MULTI_YEAR water management
#   strategies for meeting RGSM management targets.
#   The broadly-defined strategies are applied
#   to five year sequences of randomized hydrographs
#   to assess their performance under high,
#   low and medium water years.
#
#   Users upload alternative management strategies
#   depicting alternative water storage and
#   release timing strategies
#
# The code will then calculate the hydrologic
#   index (from the PCA in Walsworth
#   and Budy 2020, updated with hydrologic
#   data through 2023) for each management
#   strategy applied to each randomized hydrograph
#   and predict RGSM CPUE in each
#   reach and across the whole MRG with
#   uncertainty. 
#
# The relative performance of each alternative
#   strategy is then compared to a baseline
#   "no action" strategy, as presented in
#   Walsworth and Budy (2022).
##########################################


########################
# 0. Set directories
#    a. Data In
#    b. parameters in
#    c. Results Out
#    d. Plots out 
########################
data.wd<-".../InputData"
parm.wd<-".../ModelParameters"
out.wd<-".../ModelOutput/MultiYear"
fig.wd<-".../ModelOutput/MultiYearFigures"



################################
# 1. Load required Libraries
#################################
packages<-c("quantreg","plotrix","doParallel","doSNOW","foreach","parallel")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      
    }
    library(x, character.only = TRUE)
  }
)

################################
# 2. Source required functions
################################
setwd(parm.wd)
source("RGSM_ManagementComparisonFunctions_update2024.r") # load model functions

##################################
# 3. Load required datasets
#    a. ABQ Central hydrographs
#    b. Model parameters
#    c. PCA parameters
#    d. Mgmt strategy list
#################################
setwd(parm.wd)
load("RGSM_PCAResults_thru23.rdata") # PCA parameters
load("1_RGSM_PCMod_update2024.rdata") # MCMC posterior distributions of model params

sims.out<-mod.out$BUGSoutput$sims.list # rename object containing model parameter estimates

setwd(data.wd) # Change directory
dis<-read.csv("ABQ_CentralGage_DischargeByDay_92-23.csv",header=T) # Read historical ABQ discharge

# Read in mgmt strategies and set column classes
mgmt.list<-read.csv("ManagementStrategies_MultiYear.csv",
                    colClasses=c("factor","character","character","character","integer","integer",
                                 "integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer","integer"))

mgmt.strats<-unique(mgmt.list[,2]) # extract strategy names

# The following line extracts the column numbers for different pieces of information
mgmt.cols<-which(!names(mgmt.list)%in%c("Strategy","Description","Forecast",
                                        "Discretionary","Storage","Type","WYtype"))

type.col<-which(names(mgmt.list)=="Type") # Which column contains mgmt type info?
storage.col<-which(names(mgmt.list)=="Storage") # Which column contains storage info?
wytype.col<-which(names(mgmt.list)=="WYtype") # Which column indicates water year type?

#######################################################
# The following lines create the required year, month, day, and Julian day columns 
#  to be combined with managed hydrographs within simulations.
#########################################################
sample.dat<-dis[dis$year==2023,]
sample.month<-c(rep(1,31),rep(2,28),rep(3,31),rep(4,30),
                rep(5,31),rep(6,30),rep(7,31),rep(8,31),
                rep(9,30),rep(10,31),rep(11,30),rep(12,31))
sample.day<-c(seq(1,31),seq(1,28),seq(1,31),seq(1,30),
              seq(1,31),seq(1,30),seq(1,31),seq(1,31),
              seq(1,30),seq(1,31),seq(1,30),seq(1,31))
sample.cols<-cbind(sample.dat[,1],sample.month,sample.day,sample.dat[,2])
colnames(sample.cols)<-c("year","month","day","doy")

#############################################
# 4. Set global parameters
###############################################
n.hyd.sims<-10000 # number of randomized hydrographs to generate
yrs<-unique(dis$year) #Years with observed ABQ hydrographs
nyr<-length(yrs) # number of years of observed ABQ hydrographs
nreach<-3 # number of reaches
hyd.brks<-c(round(n.hyd.sims/3),round(2*n.hyd.sims/3)) # break point indices between low, medium, and high water years
mcsims<-100# number of mcmc posterior samples to use for model parameters
nsims<-100 # number of sets of hydrographs to simulate RGSM CPUE for
nsamps<-10 # numner of sample sites in each reach
ndays<-365 # number of days in a year
meltstart<-60 # Start of March DOY
meltend<-210  # End of June DOY
nmgmt<-length(mgmt.strats) # number of mgmt strategies
wy.types<-c("Low","Mid","High") # Water year types
ntype<-length(wy.types) # Number of water year types

##############################################
# Set up management strategy parameters
# - scls: the proportion of water stored
#         during months specified for storage
# - damt.acft: Amount of discretionary water
#              in acre-feet
# - sl.acft: Storage capacity in acre feet
# - thresh: CPUE thresholds to examine
##############################################
scls<-c(.1,.25,.5) # Proportions of flow to be stored in different mgmt scenarios
discamt.acft<-c(10000,30000,100000)  # Discretionary water in ACRE FEET
discamt<-discamt.acft*0.5042    # converts acre feet to cubic feet per second (for 24 hours)
sl.acft<-c(25000,50000,100000)  # Storage capacity in acft
storage.limits<-sl.acft*0.5042  # converts acre feet to cubic feet per second (for 24 hours)
thresh<-seq(0.1,5,by=.1) # CPUE thresholds
nthresh<-length(thresh) # number of cpue thresholds
nyr.sim<-5 # Number of years in each simulation

#########################################################
# 5. Generate matrix of simulated, randomized hydrographs
#########################################################
# The following line calculates the total volume (in cfs per day) observed each year
histdist<-aggregate(dis$discharge,list("Year"=dis$year),sum)

# The following section scales randomly selected hydrographs to the total flow
#   volume of another randomly selected year with stochastic variation
#   added (by multiplying the total flow volume by a random number with mean 0,
#   and sd 0.05)

hydsimsout<-matrix(NA,nrow=n.hyd.sims,ncol=ndays)
for(i in 1:n.hyd.sims){
  yruse<-sample(yrs,1)
  yrdis<-dis$discharge[dis$year==yruse]*exp(rnorm(1,-.02^2/2,.02))
  hydsimsout[i,]<-yrdis
}

# The following line orders all daily discharge values in decreasing order.
# The top row of the resulting matrix includes the lowest discharge values for
#  each day of the year, while the last line includes the highest discharge
#  value for each day
sims.ord<-apply(hydsimsout,2,sort,decreasing=F)

# The following section generates randomized hydrographs by:
#  - randomly selecting a starting rank value for the January 1 discharge.
#  - Selecting the index of January 2 discharge by applying a random
#    walk to the January 1 discharge rank.
#  - the random walk step is repeated for each day of the year

rw.dat<-matrix(NA,nrow=ndays,ncol=n.hyd.sims)
for(j in 1:n.hyd.sims){
  rw<-rep(0,ndays) # create empty storage for rank flow indices
  rw[1]<-sample(1:n.hyd.sims,1)# Select January 1 rank flow index
  for(i in 2:ndays){
    rw[i]<-round(rw[i-1]+rnorm(1,0,.01*n.hyd.sims)) # Random walk
    if(rw[i]>n.hyd.sims) rw[i]<-n.hyd.sims # Check for invalid indices (out of bounds)
    if(rw[i]<1) rw[i]<-1# Check for invalid indices (out of bounds)
  }
  
  rw.dat[,j]<-sims.ord[cbind(rw,1:ncol(hydsimsout))]# Create matrix of randomized hydrograph discharge values
  }

# The following line calculates the total flow volume from March through June for
#  each randomized hydrograph so they can be ranked for determining low, medium, 
#  and high water years.
totflows<-colSums(rw.dat[meltstart:meltend,],)

############################################
# 6. Set up parallel processing and run model
#  across all scenario combinations
############################################
paramcombos<-expand.grid(1:length(storage.limits),1:length(scls)) # Create matrix of parameter index combinations


cl<-makeCluster(min(detectCores()-1,length(paramcombos[,1]))) # Set up parallel computing
registerDoParallel(cl) # Begin parallel computing
start.time<-Sys.time() # Start model timer (not necessary but can be helpful)

# The following lines run the model across all parameter combinations. It will
#  run in parallel (running multiple scenarios simultaneously on different
#  processor cores). If you want to run the scenarios in sequence, change the
#  %dopar% to %do% on the following line.
mod.outall<-foreach(qqq=1:nrow(paramcombos),.packages=c("quantreg","plotrix")) %dopar% {
  
  print(qqq) # Print counter
  storagelim<-storage.limits[paramcombos[qqq,1]] # Set temporary storage limit
  mgmt.scale<-scls[paramcombos[qqq,2]] # Set proportional managmement scale (what proportion of water is stored?)
  
  for(d in 1:length(discamt.acft)){ # for each discretionary water volume

    ###################################
    # 7. Create storage matrices/arrays
    ###################################
    cpue.pred<-array(0,dim=c(nmgmt,nsims*mcsims,nreach,3,nyr.sim)) # Create storage for reach-specific CPUE
    annests<-array(0,dim=c(nmgmt,nsims*mcsims,3,nyr.sim)) # Create storage for MRG wide CPUE 
    annflows<-array(0,dim=c(nsims,nyr.sim,ndays)) # create storage for annual flows
    pcs.save<-array(0,dim=c(3,nmgmt,nsims,nyr.sim)) # create storage for PC indices
    wys<-matrix(NA,nrow=nsims,ncol=nyr.sim) # create storage for water year indices
    types<-wys # create storage for water year types
    dbg<-F # Run in debugging mode? T/F
    
    #################################
    # 8. Run model across randomized hydrograph combinations
    #################################
    for(i in 1:nsims){ # For each combination of hydrographs
      
      # The following line will print a counter (will not print when run in parallel)
      if(i%%5 == 0) print(paste("s = ",paramcombos[qqq,1]," of ",length(storage.limits)," p = ",paramcombos[qqq,2]," of ",length(scls),"d = ",d," of ",length(discamt.acft),", i = ",i," of ",nsims,sep=""))
      
      wys[i,]<-sample(1:n.hyd.sims,nyr.sim,replace=T) # sample water year indices
      types[i,as.integer(wys[i,])<(n.hyd.sims/4)]<-"low" # store water year types for low water years
      types[i,as.integer(wys[i,])>(n.hyd.sims/4) & as.integer(wys[i,])<(3*n.hyd.sims/4)]<-"mid" # store water year type for medium water years
      types[i,as.integer(wys[i,])>(3*n.hyd.sims/4)]<-"high" # store water year types for high water years
      
      if(dbg) print("A")# Debugging checkpoint
      
      for(k in 2:nyr.sim){ # for each year
        
        # Set lower water years following higher water years as specific 
        #   year types (e.g., a low water year following a high water year
        #   is designated "highlow")
        if(types[i,k-1]=="high"&types[i,k]=="low") types[i,k]<-"highlow" 
        if(types[i,k-1]=="high"&types[i,k]=="mid") types[i,k]<-"highmid"
        if(types[i,k-1]=="mid"&types[i,k]=="low") types[i,k]<-"midlow"
      }
      
      if(dbg) print("B")# Debugging checkpoint
      
      rn.seeds<-sample(seq(1,1000000),mcsims) # Sample indices for rng seeds
      
      for(y in 1:nyr.sim){
        baselineflows<-rw.dat[,order(totflows)[wys[i,y]]] # Extract randomized hydrograph
        annflows[i,y,]<-baselineflows # Store annual hydrograph
        if(dbg) print("C")# Debugging checkpoint
      }
      
      if(dbg) print("D")# Debugging checkpoint
      
      for(m in 1:nmgmt){  # For each management strategy
        
        if(dbg) print("E")# Debugging checkpoint
        
        stored.avail<-0 # Set initial available stored water
        
        # Following line will check if there is any use of across year storage in this strategy
        storage<-sum(abs(mgmt.list[mgmt.list[,2]==mgmt.strats[m],storage.col]))>0 
        
        ############################
        # The following section will loop through each year being simulated,
        #  convert the randomized hydrograph based on the specified mgmt 
        #  strategy, calculate the integrated hydrologic metrics (PC values)
        #  for the managed hydrograph, and stochastically simulate CPUE values
        ############################
        for(y in 1:nyr.sim){
          
          mgmt.seq<-mgmt.list[mgmt.list[,2]==mgmt.strats[m],] # Extracts mgmt scenario information
          mgmt<-mgmt.seq[mgmt.seq[,wytype.col]==types[i,y],] # Extracts mgmt strategy for current water year type
          storage<-abs(mgmt[,storage.col])>0 # Determines if across-year storage occurs in this strategy
          
          # The following line applies the supplied management strategy to the annual hydrograph
          mgmt.hyd<-convert.hydro.mgmt(hydro=annflows[i,y,],strat = as.numeric(as.character(mgmt[,mgmt.cols])),scale = mgmt.scale,
                                       discretionary = mgmt[,which(names(mgmt.list)=="Discretionary")],discamt = discamt[d],
                                       type=mgmt[,type.col],storage=storage,storage.limit = storagelim,
                                       stored.in = stored.avail,dbg=F)
          if(dbg) print("F")# Debugging checkpoint
          
          # The following line creates the hydrograph matrix needed for calculating annual PC score,
          #  and sets the available temporarily stored water
          if(storage){
            mgmt.hyd.table<-as.data.frame(cbind(sample.cols,mgmt.hyd[[1]]))
            stored.avail<-mgmt.hyd[[2]]
          }
          
          # The following line creates the hydrograph matrix needed for calculating annual PC score
          if(!storage) mgmt.hyd.table<-as.data.frame(cbind(sample.cols,mgmt.hyd))
          
          colnames(mgmt.hyd.table)[5]<-"discharge" # Set discharge column name
          
          ##############################################
          # 2. Calculate PC score for managed hydrograph
          ##############################################
          pcs<-calc.pcscores(daily.avs=mgmt.hyd.table,drying=-1,pca1=pca1)
          pcs.save[,m,i,y]<-pcs[c(1,3,5),1] # Save 10th, 50th, 90th percentile drying PC scores
          
          if(dbg) print("G") # Debugging checkpoint
          
          # The following section checks for an error that was fixed in earlier version of code.
          #  It should not be necessary but has been left in for now.
          if(length(stored.avail)>1){
            print(m)
            stop()
          }
          
          for(k in 1:mcsims){ #Across all MCMC parameter samples being examined
            set.seed(rn.seeds[k]) # set seed
            
            mcuse<-sample(1:length(sims.out$alpha.c[,1]),1) # randomly select a mcmc posterior draw to use for this simulation
            
            
            whichts<-sample(1:31,1) # Sample index for random wlak component
            tsuse<-sims.out$ts[mcuse,whichts]+rnorm(1) # Random walk value for starting year
            for(y in 1:nyr.sim){ # For each year being simulated
              tsuse<-c(tsuse,tsuse[length(tsuse)]+rnorm(1)) # Generate random walk component for hurdle model
            }                                  
            for(y in 1:nyr.sim){ #For each year
              fldsim<-pcs.save[,m,i,y] # PC scores for predicting RGSM CPUE
              for(j in 1:3){ # for each reach
                pred.ca<-sims.out$alpha.c[mcuse,j]*exp(-sims.out$beta.off[mcuse]*exp(-sims.out$beta.denspred[mcuse]*fldsim)) # Density component of model
                pred.ca[pred.ca==0]<-.00000000001 # Add small value to zeros
                predobs<-sims.out$alpha.p[mcuse,j]+sims.out$beta.prespred[mcuse]*fldsim+tsuse[y] # hurdle component of model
                predobs<-exp(predobs)/(1+exp(predobs)) # Anti-logit transform
                prate<-(pred.ca+sqrt(pred.ca^2+4*(sims.out$cv[mcuse]*pred.ca)^2))/(2*(sims.out$cv[mcuse]*pred.ca)^2) # calculate gamma rate
                psh<-1+pred.ca*prate # calculate gamma shape
                predfish<-c() # Create empty storage object
                predfish[1]<-mean(rbinom(10,1,prob=predobs[1])*rgamma(10,shape=psh[1],rate=prate[1])) # Predict low drying CPUE
                predfish[2]<-mean(rbinom(10,1,prob=predobs[2])*rgamma(10,shape=psh[2],rate=prate[2])) # Predict median drying cpue
                predfish[3]<-mean(rbinom(10,1,prob=predobs[3])*rgamma(10,shape=psh[3],rate=prate[3])) # Predict high drying cpue
                
                cpue.pred[m,((i-1)*mcsims+k),j,,y]<-predfish # Store for CPUE for reach j
                
              }
              annests[m,((i-1)*mcsims+k),,y]<-colMeans(cpue.pred[m,((i-1)*mcsims+k),,,y]) # Store MRG-wide CPUE
            }
          }
        }  
      }
    }
    setwd(out.wd) # Changes working directory
    
    # The following line saves necessary simulation output to file
    save(annests,cpue.pred,annflows,pcs.save,wys,types,file=paste("RGSM_MultiYr_s=",paramcombos[qqq,1],"_p=",paramcombos[qqq,2],"_d=",d,".rdata",sep=""),
         compress=T)
    rm(annests,cpue.pred,annflows,pcs.save,wys,types) # removes storage objects from workspace
  }
}

registerDoSEQ() # Closes parallel computing cluster
stopCluster(cl)
Sys.time()-start.time # Prints total elapsed time during model runs (not required)

#########################################
# Model Output Figures
# Read in model output for each scenario
#########################################
setwd(out.wd)
s1files<-list.files()[as.numeric(substr(list.files(),16,16))==1] # Storage scenario 1 files
s2files<-list.files()[as.numeric(substr(list.files(),16,16))==2] # Storage scenario 2 files
s3files<-list.files()[as.numeric(substr(list.files(),16,16))==3] # Storage scenario 3 files

fillist<-list(s1files,s2files,s3files) # List model output files

# Following line creates storage matrices for summary output
my.out.1<-my.out.5<-matrix(NA,nrow=length(mgmt.strats)*3*3*3,ncol=4+nyr.sim+2)

#############################################################
# The following lines create a pdf with model output figures
#############################################################
setwd(data.wd)
# Read in color indicators for figures (can make similar strategies the same color)
mgcols<-read.csv("ManagementStrategies_MultiYear_ColorLineIndicators.csv",header=T)

setwd(fig.wd) # Set working directory for figure outputs
options(scipen=10000000) # Change settings so figure axes do not use Scientific notation

pdf("RGSM_MultiYear_Simulation_Results_3x3.pdf",height=8,width=9)
layout(matrix(c(7,8,9,
                4,5,6,
                1,2,3),nrow=3,ncol=3,byrow=T))
par(oma=c(6,6,2,0),mar=c(.1,.1,.1,.1))

cols<-c("dodgerblue","darkorange","darkgreen","black")

# The following section produces three figures, each with 9 panels (3x3)
# demonstrating the expected number of years each strategy achieves the
# management target of 1 RGSM per 100m2.
ctr<-0
for(i in 1:length(fillist)){
  for(p in 1:length(scls)){
    for(d in 1:length(discamt.acft)){
      setwd(out.wd)
      load(file=fillist[[i]][as.numeric(substr(fillist[[i]],20,20))==p & 
                               as.numeric(substr(fillist[[i]],24,24))==d])

      mgmt.perf1<-mgmt.perf5<-matrix(NA,nrow=dim(annests)[1],ncol=dim(annests)[4]+1)
      plot(NA,NA,ylim=c(0,1),xlim=c(0,dim(annests)[4]),ylab="",
           xlab="",xaxt="n",yaxt="n")
      axis(side=4,padj=1,labels=F)
      if(p==1) axis(side=1)
      if(p==1) mtext(side=1,line=3,"Number of Years Meeting\nCPUE Target of 1",cex=.75)
      if(d==1) mtext(side=2,line=2.5,"Proportion of Simulations",cex=.75)
      if(d==1) axis(side=2,las=2)
      for(j in 1:dim(annests)[1]){
        numsuccess1<-rowSums(annests[j,,3,]>1)
        numsuccess5<-rowSums(annests[j,,3,]>5)
        for(y in 0:dim(annests)[4]){
          mgmt.perf1[j,y+1]<-sum(numsuccess1>=y)/dim(annests)[2]
          mgmt.perf5[j,y+1]<-sum(numsuccess5>=y)/dim(annests)[2]
          
        }
        my.out.1[ctr*dim(annests)[1]+j,1]<-sl.acft[i]
        my.out.1[ctr*dim(annests)[1]+j,2]<-scls[p]
        my.out.1[ctr*dim(annests)[1]+j,3]<-discamt.acft[d]
        my.out.1[ctr*dim(annests)[1]+j,4]<-mgmt.strats[j]
        my.out.1[ctr*dim(annests)[1]+j,5:(ncol(my.out.1)-1)]<-mgmt.perf1[j,]
        my.out.5[ctr*dim(annests)[1]+j,1]<-sl.acft[i]
        my.out.5[ctr*dim(annests)[1]+j,2]<-scls[p]
        my.out.5[ctr*dim(annests)[1]+j,3]<-discamt.acft[d]
        my.out.5[ctr*dim(annests)[1]+j,4]<-mgmt.strats[j]
        my.out.5[ctr*dim(annests)[1]+j,5:(ncol(my.out.5)-1)]<-mgmt.perf5[j,]
        if(is.na(mgcols[j,2])) mgcols[j,2]<-ifelse(j%%4==0,4,j%%4)
        lines(seq(0,dim(annests)[4]),mgmt.perf1[j,],col=cols[mgcols[j,2]])
        
      }
      text(0,.2,labels=paste("Best = ",mgmt.strats[which.max(rowSums(mgmt.perf1))]),
           pos=4)
      text(0,.1,labels=paste("Worst = ",mgmt.strats[which.min(rowSums(mgmt.perf1))]),
           pos=4)
      ranks1<-rank(rowSums(mgmt.perf1),ties.method = "average")
      ranks5<-rank(rowSums(mgmt.perf5),ties.method = "average")
      my.out.1[((ctr*dim(annests)[1])+1):((ctr+1)*dim(annests)[1]),ncol(my.out.1)]<-ranks1
      my.out.5[((ctr*dim(annests)[1])+1):((ctr+1)*dim(annests)[1]),ncol(my.out.1)]<-ranks5
      ctr<-ctr+1
    }
  }
  
  mtext(outer=T,side=1,line=4.5,at=c(.15,.5,.85),c("10K Discretionary","30K Discretionary","100K Discretionary"))
  mtext(outer=T,side=2,line=4,at=c(.16,.5,.84),c("Shift 10%","Shift 25%","Shift 50%"))
  mtext(outer=T,side=3,line=.5,paste("Storage = ",sl.acft[i]," acft",sep=""))
}


# The following section produces three figures, each with 9 panels (3x3)
# demonstrating the expected difference in number of years achieving the
# management target of 1 RGSM per 100m2 between each
# alternative multi-year strategy and the No action strategy.

layout(matrix(c(7,8,9,
                4,5,6,
                1,2,3),nrow=3,ncol=3,byrow=T))
par(oma=c(6,7,2,0),mar=c(.1,.1,.1,.1))

for(i in 1:length(fillist)){
  for(p in 1:length(scls)){
    for(d in 1:length(discamt.acft)){
      
      load(file=fillist[[i]][as.numeric(substr(fillist[[i]],20,20))==p & 
                               as.numeric(substr(fillist[[i]],24,24))==d])

      mgmt.perf1<-mgmt.perf5<-matrix(NA,nrow=dim(annests)[1],ncol=dim(annests)[4]+1)

      plot(NA,NA,ylim=c(-.1,.3),xlim=c(0,dim(annests)[4]),ylab="",
           xlab="",xaxt="n",yaxt="n")
      axis(side=4,padj=1,labels=F)
      if(p==1) axis(side=1)
      if(p==1) mtext(side=1,line=3,"Number of Years Meeting\nCPUE Target of 1",cex=.75)
      if(d==1) mtext(side=2,line=3,"Increased Proportion\nover No Action",cex=.75)
      if(d==1) axis(side=2,las=2)
      for(j in 1:dim(annests)[1]){
        numsuccess1<-rowSums(annests[j,,3,]>1)
        numsuccess5<-rowSums(annests[j,,3,]>5)
        for(y in 0:dim(annests)[4]){
          mgmt.perf1[j,y+1]<-sum(numsuccess1>=y)/dim(annests)[2]
          mgmt.perf5[j,y+1]<-sum(numsuccess5>=y)/dim(annests)[2]
        }
        
        
        lines(seq(0,dim(annests)[4]),mgmt.perf1[j,]-mgmt.perf1[1,],col=cols[mgcols[j,2]])
        
      }
      text(0,.29,labels=paste("Best = ",mgmt.strats[which.max(rowSums(mgmt.perf1))]),
           pos=4)
      text(0,.24,labels=paste("Worst = ",mgmt.strats[which.min(rowSums(mgmt.perf1))]),
           pos=4)
      
      abline(h=0,lwd=2)  
    }
    
  }
  
  mtext(outer=T,side=1,line=4.5,at=c(.15,.5,.85),c("10K Discretionary","30K Discretionary","100K Discretionary"))
  mtext(outer=T,side=2,line=5.5,at=c(.16,.5,.84),c("Shift 10%","Shift 25%","Shift 50%"))
  mtext(outer=T,side=3,line=.5,paste("Storage = ",sl.acft[i]," acft",sep=""))
}
# The following section produces three figures, each with 9 panels (3x3)
# demonstrating the expected number of years each strategy achieves the
# downlisting target of 5 RGSM per 100m2.
for(i in 1:length(fillist)){
  for(p in 1:length(scls)){
    for(d in 1:length(discamt.acft)){
      setwd(out.wd)
      load(file=fillist[[i]][as.numeric(substr(fillist[[i]],20,20))==p & 
                               as.numeric(substr(fillist[[i]],24,24))==d])
      
      mgmt.perf1<-mgmt.perf5<-matrix(NA,nrow=dim(annests)[1],ncol=dim(annests)[4]+1)

      plot(NA,NA,ylim=c(0,1),xlim=c(0,dim(annests)[4]),ylab="",
           xlab="",xaxt="n",yaxt="n")
      axis(side=4,padj=1,labels=F)
      if(p==1) axis(side=1)
      if(p==1) mtext(side=1,line=3,"Number of Years Meeting\nCPUE Target of 5",cex=.75)
      if(d==1) mtext(side=2,line=2.5,"Proportion of Simulations",cex=.75)
      if(d==1) axis(side=2,las=2)
      for(j in 1:dim(annests)[1]){

        numsuccess5<-rowSums(annests[j,,3,]>5)
        for(y in 0:dim(annests)[4]){

          mgmt.perf5[j,y+1]<-sum(numsuccess5>=y)/dim(annests)[2]
        }
        
        lines(seq(0,dim(annests)[4]),mgmt.perf5[j,],col=cols[mgcols[j,2]])
        
      }
      text(0,.2,labels=paste("Best = ",mgmt.strats[which.max(rowSums(mgmt.perf5))]),
           pos=4)
      text(0,.1,labels=paste("Worst = ",mgmt.strats[which.min(rowSums(mgmt.perf5))]),
           pos=4)
      
      
    }
  }
  
  mtext(outer=T,side=1,line=4.5,at=c(.15,.5,.85),c("10K Discretionary","30K Discretionary","100K Discretionary"))
  mtext(outer=T,side=2,line=4,at=c(.16,.5,.84),c("Shift 10%","Shift 25%","Shift 50%"))
  mtext(outer=T,side=3,line=.5,paste("Storage = ",sl.acft[i]," acft",sep=""))
}


# The following section produces three figures, each with 9 panels (3x3)
# demonstrating the expected difference in number of years achieving the
# downlisting target of 5 RGSM per 100m2 between each
# alternative multi-year strategy and the No action strategy.
layout(matrix(c(7,8,9,
                4,5,6,
                1,2,3),nrow=3,ncol=3,byrow=T))
par(oma=c(6,7,2,0),mar=c(.1,.1,.1,.1))

for(i in 1:length(fillist)){
  for(p in 1:length(scls)){
    for(d in 1:length(discamt.acft)){
      
      load(file=fillist[[i]][as.numeric(substr(fillist[[i]],20,20))==p & 
                               as.numeric(substr(fillist[[i]],24,24))==d])
      
      mgmt.perf1<-mgmt.perf5<-matrix(NA,nrow=dim(annests)[1],ncol=dim(annests)[4]+1)

      plot(NA,NA,ylim=c(-.1,.3),xlim=c(0,dim(annests)[4]),ylab="",
           xlab="",xaxt="n",yaxt="n")
      axis(side=4,padj=1,labels=F)
      if(p==1) axis(side=1)
      if(p==1) mtext(side=1,line=3,"Number of Years Meeting\nCPUE Target of 5",cex=.75)
      if(d==1) mtext(side=2,line=3,"Increased Proportion\nover No Action",cex=.75)
      if(d==1) axis(side=2,las=2)
      for(j in 1:dim(annests)[1]){
        numsuccess5<-rowSums(annests[j,,3,]>5)
        for(y in 0:dim(annests)[4]){
          mgmt.perf5[j,y+1]<-sum(numsuccess5>=y)/dim(annests)[2]
        }
        
        
        lines(seq(0,dim(annests)[4]),mgmt.perf5[j,]-mgmt.perf5[1,],col=cols[mgcols[j,2]])
        
      }
      
      text(0,.29,labels=paste("Best = ",mgmt.strats[which.max(rowSums(mgmt.perf5))]),
           pos=4)
      text(0,.24,labels=paste("Worst = ",mgmt.strats[which.min(rowSums(mgmt.perf5))]),
           pos=4)
      
      abline(h=0,lwd=2)  
    }
    
  }
  
  mtext(outer=T,side=1,line=4.5,at=c(.15,.5,.85),c("10K Discretionary","30K Discretionary","100K Discretionary"))
  mtext(outer=T,side=2,line=5.5,at=c(.16,.5,.84),c("Shift 10%","Shift 25%","Shift 50%"))
  mtext(outer=T,side=3,line=.5,paste("Storage = ",sl.acft[i]," acft",sep=""))
}

dev.off() # Closes PDF Plotting Environment




# The following lines set column names for summary output matrices
colnames(my.out.1)<-c("Storage","ProportionMoved","Discretionary",
                      "Strategy",0:5,"Rank")
colnames(my.out.5)<-colnames(my.out.1)

my.out.1<-as.data.frame(my.out.1)# Convert to data frame

# The following line calculates the median and 90% simulation interval for rank
# performance
rankquants<-aggregate(as.numeric(my.out.1$Rank),list(my.out.1$Strategy),quantile,probs=c(.05,.5,.95))
rankquants<-rankquants[order(rankquants[,2][,1]),]# Rank by the worst case performance

###################################
# The following section creates a pdf with the rank performance of 
# alternative multi-year management strategies (as in Figure 10 of Report)
###################################
setwd(fig.wd)
pdf("RGSM_MultiYear_RankStrategies.pdf",height=5,width=6)
par(mar=c(9,5,1,1))
plot(NA,NA,ylim=c(0,nmgmt+.5),xlim=c(0,length(mgmt.strats)+1),ylab="Rank (Higher = Better)",
     xlab="",xaxt="n",yaxs="i",las=2,xaxs="i")

segments(x0=seq(1,length(mgmt.strats)),y0=rankquants[,2][,1],
         x1=seq(1,length(mgmt.strats)),y1=rankquants[,2][,3],lwd=2)
points(seq(1,length(mgmt.strats)),rankquants[,2][,2],pch=16)
segments(x0=seq(1,length(mgmt.strats)),y0=rankquants[,2][,1],x1=seq(1,length(mgmt.strats)),y1=0,
         col="grey80",lwd=0.7)

axis(side=1,at=seq(1,length(mgmt.strats)),labels=rankquants[,1],las=2,
     cex.axis=.75)

dev.off()

