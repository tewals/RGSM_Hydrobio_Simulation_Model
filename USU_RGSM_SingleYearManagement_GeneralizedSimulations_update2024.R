###########################################
# RGSM Generalized Hydrologic Management Strategy
#   Simulations
#
# Coded by: Tim Walsworth
# e-mail: timothy.walsworth@usu.edu
#
###########################################

###########################################
# This code will simulate the predicted
#   relative performance of different
#   broadly defined water management strategies
#   for meeting RGSM management targets.
#   The broadly defined strategies are applied
#   to randomized hydrographs to assess their
#   performance under high, low and medium water
#   years.
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
out.wd<-".../ModelOutput/SingleYearGeneralized"
fig.wd<-".../ModelOutput/SingleYearGeneralizedFigures"


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
mgmt.list<-read.csv("ManagementStrategies_SingleYear.csv",
                    colClasses=c("factor","character","character","integer",
                                 "integer","integer","integer","integer",
                                 "integer","integer","integer","integer",
                                 "integer","integer","integer","integer"))

mgmt.strats<-unique(mgmt.list[,2]) # extract strategy names

# The following line extracts the column numbers for different pieces of information
mgmt.cols<-which(!names(mgmt.list)%in%c("Strategy","Description","Forecast",
                                        "Discretionary","Storage","Type"))

mgmt.type<-which(names(mgmt.list)=="Type")# Which column contains mgmt type info?

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

####################################
# 4. Set global parameters
####################################
n.hyd.sims<-10000 #Number of randomized hydrographs to generate
yrs<-unique(dis$year) # observed discharge years
nyr<-length(yrs) # number of observed discharge years
nreach<-3 # number of reaches

# Set breaks for water year type. Currently set to lowest 40% are low water, 
#  middle 20% are medium, and highest 40% are high water years. Reflects 
#  pattern of generally high or low water years, with fewer intermediate
#  years.
hyd.brks<-c(round(.4*n.hyd.sims),round(.6*n.hyd.sims)) 
mcsims<-100 #Number of parameter values from mcmc posterior distributions to use
nsims<-100  # Number of hydrograph combinations to simulate across
nsamps<-10  # Number of fish sampling events in each reach
ndays<-365  # Number of days in a year
meltstart<-60 # Start of spring snowmelt (for ranking spring flow volumes)
meltend<-210  # End of spring snowmelt (for ranking spring flow volumes)
nmgmt<-length(mgmt.strats) # Number of alternative management strategies
wy.types<-c("Low","Mid","High") # Water year types (need to match inputs)
ntype<-length(wy.types) # How many water year types?

##############################################
# Set up management strategy parameters
# - scls: the proportion of water stored
#         during months specified for storage
# - damt.acft: Amount of discretionary water
#              in acre-feet
# - sl.acft: Storage capacity in acre feet
# - thresh: CPUE thresholds to examine
##############################################
scls<-c(.1,.25,.5) # Proportions of water to be stored during storage periods
damt.acft<-c(10000,30000,100000)  # Discretionary water in ACRE FEET
sl.acft<-c(25000,50000,100000)  # Storage capacity in acft
storage.limits<-sl.acft*0.5042 # converts acre feet to cubic feet per second (for 24 hours)
thresh<-seq(0.1,5,by=.1) # CPUE thresholds to examine
nthresh<-length(thresh) # Number of CPUE thresholds being considered

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
mod.outall<-foreach(qqq=1:length(paramcombos[,1]),.packages=c("quantreg","plotrix")) %dopar% {
  set.seed(1)
  print(qqq) # Print counter
  
  storagelim<-storage.limits[paramcombos[qqq,1]] # Set temporary storage limit

  mgmt.scale<-scls[paramcombos[qqq,2]]# Set proportional managmement scale (what proportion of water is stored?)
  
  for(d in 1:length(damt.acft)){ # for each discretionary water volume
    discamt<-damt.acft[d]*0.5042 # converts acre feet to cubic feet per second (for 24 hours)

    cpue.pred<-array(0,dim=c(nmgmt,nsims*mcsims,nreach,3,ntype))# Create storage for reach-specific CPUE
    annests<-array(0,dim=c(nmgmt,nsims*mcsims,3,ntype)) # Create storage for MRG wide CPUE 
    annflows<-array(0,dim=c(nsims,ntype,ndays)) # create storage for annual flows
    fl.perc.save<-matrix(0,nrow=nsims,ncol=ntype) # create storage for annual flows percentile
    pcs.save<-array(0,dim=c(3,nmgmt,nsims,ntype)) # create storage for PC indices
    
    ###############################
    # 7. Simulate across all water year types
    ################################
    for(wy in 1:ntype){

      yrtype<-wy.types[wy]    # Set water year type to low, mid, or high
      
      ####################
      # Run model across selected randomized hydrographs
      ####################
      for(j in 1:nsims){

        if(j%%5 == 0) print(paste("p = ",paramcombos[qqq,i]," of ",length(scls),"d = ",d," of ",length(damt.acft)," wy = ",wy," of ",ntype,", j = ",j," of ",nsims,sep=""))

        # The following line selects a randomized hydrograph based on water year type
        flnum<-switch(yrtype,
                      Low = sample(seq(1,hyd.brks[1]),1),
                      Mid =  sample(seq(hyd.brks[1]+1,hyd.brks[2]),1),
                      High = sample(seq(hyd.brks[2]+1,n.hyd.sims),1))
        
        fl.perc.save[j,wy]<-flnum/n.hyd.sims # store flow relative to all randomized hydrographs
        
        baselineflows<-rw.dat[,order(totflows)[flnum]] # Extract randomly selected hydrograph
        annflows[j,wy,]<-baselineflows # Store extracted hydrograph
        
        ######################################
        #   Apply each Management strategy
        ######################################
        rn.seeds<-sample(seq(1,1000000),mcsims)
        for(m in 1:nmgmt){

          mgmt<-mgmt.list[mgmt.list[,2]==mgmt.strats[m],] # Extracts mgmt strategy 
          
          # The following line applies the supplied management strategy to the annual hydrograph
          mgmt.hyd<-convert.hydro.mgmt(hydro=baselineflows,strat = as.numeric(mgmt[,mgmt.cols]),scale = mgmt.scale,
                                       discretionary = mgmt[,which(names(mgmt)=="Discretionary")],discamt = discamt,
                                       type=mgmt[,mgmt.type],dbg=F,storage=F,storage.limit = storagelim,
                                       stored.in = 0)
          # The following line creates the hydrograph matrix needed for calculating annual PC score
          mgmt.hyd.table<-as.data.frame(cbind(sample.cols,mgmt.hyd))
          colnames(mgmt.hyd.table)[5]<-"discharge" # Set discharge column name
          
          ##############################################
          # 2. Calculate PC score for managed hydrograph
          ##############################################
          pcs<-calc.pcscores(daily.avs=mgmt.hyd.table,drying=-1,pca1=pca1)
          pcs.save[,m,j,wy]<-pcs[c(1,3,5),1]
          
          for(k in 1:mcsims){ # Predict CPUE Across all MCMC parameter samples being examined
         
            set.seed(rn.seeds[k])
            
            mcuse<-sample(1:length(sims.out$alpha.c[,1]),1) # randomly select a mcmc posterior draw to use for this simulation
  
            fldsim<-pcs[,1] # PC scores for predicting RGSM CPUE
            whichts<-sample(1:31,1) # Sample index for random walk component
            tsuse<-sims.out$ts[mcuse,whichts]+rnorm(1) # Random walk value for starting year
            
            for(i in 1:3){ # for each reach
              pred.ca<-sims.out$alpha.c[mcuse,i]*exp(-sims.out$beta.off[mcuse]*exp(-sims.out$beta.denspred[mcuse]*fldsim))# Density component of model
              pred.ca[pred.ca==0]<-.00000000001 # Add small value to zeros
              predobs<-sims.out$alpha.p[mcuse,i]+sims.out$beta.prespred[mcuse]*fldsim+tsuse# hurdle component of model
              predobs<-exp(predobs)/(1+exp(predobs)) # Anti-logit transform
              prate<-(pred.ca+sqrt(pred.ca^2+4*(sims.out$cv[mcuse]*pred.ca)^2))/(2*(sims.out$cv[mcuse]*pred.ca)^2) # calculate gamma rate
              psh<-1+pred.ca*prate# calculate gamma shape
              predfish<-c() # Create empty storage object
              predfish[1]<-mean(rbinom(10,1,prob=predobs[1])*rgamma(10,shape=psh[1],rate=prate[1])) # Predict low drying CPUE
              predfish[2]<-mean(rbinom(10,1,prob=predobs[2])*rgamma(10,shape=psh[2],rate=prate[2])) # Predict median drying cpue
              predfish[3]<-mean(rbinom(10,1,prob=predobs[3])*rgamma(10,shape=psh[3],rate=prate[3])) # Predict high drying cpue
              cpue.pred[m,((k-1)*nsims+j),i,,wy]<-predfish # Store for CPUE for reach j
              
            }
            annests[m,((k-1)*nsims+j),,wy]<-colMeans(cpue.pred[m,((k-1)*nsims+j),,,wy])   # Store MRG-wide CPUE
             
          }
          
        }
        
      }
      
    }
    setwd(out.wd) # Changes working directory
    
    # The following line saves necessary simulation output to file
    save(annests,cpue.pred, annflows,fl.perc.save,
         pcs.save,file=paste("SingleYear_GeneralizedMgmt_s=",paramcombos[qqq,1],"_p=",paramcombos[qqq,2],"_d=",d,".rdata",sep=""),compress=T)
    
    rm("annests","cpue.pred","annflows","fl.perc.save","pcs.save")# removes storage objects from workspace
    
  }
  
  # }
}

registerDoSEQ() # Closes parallel computing cluster
Sys.time()-start.time # Prints total elapsed time during model runs (not required)

#########################
# Model Figures
#########################
sy.out.1<-sy.out.5<-matrix(NA,nrow=nmgmt*3*3*3*3,ncol=7) # Create storage objects for output summary data
ctr<-0 # start counter object
setwd(out.wd) # set directory

################
# The following section will loop through simulation model output files
#  to calculate summary results
###############
for(vv in 1:length(list.files())){
  
  s<-as.numeric(substr(list.files()[vv],30,30)) #extract storage index
  p<-as.numeric(substr(list.files()[vv],34,34)) # extract proportional mgmt index
  d<-as.numeric(substr(list.files()[vv],38,38)) # extract discretionary water index
  
  dacft<-damt.acft[d] # Set discretionary water acft
  scl<-scls[p] # set proportional mgmt value
  slim<-sl.acft[s] # Set storage limit
  
  load(list.files()[vv]) # load model output
  
  success1<-success5<-array(0, dim=c(nmgmt,5,ntype)) # create storage objects
  
  for(i in 1:nmgmt){ # for each mgmt strategy
    for(j in 1:3){ # for each drying scenario
      for(k in 1:ntype){ # for each water year type
        success1[i,j,k]<-sum(annests[i,,j,k]>1)/length(annests[i,,j,k]) # proportion of simulations meeting managment target
        success5[i,j,k]<-sum(annests[i,,j,k]>5)/length(annests[i,,j,k]) # proportion of simulations meeting downlisting target
        
        # the following section stores summary output
        if(j == 3){
          sy.out.5[ctr*ntype+k,1]<-mgmt.strats[i]
          sy.out.5[ctr*ntype+k,2]<-scls[p]
          sy.out.5[ctr*ntype+k,3]<-damt.acft[d]
          sy.out.5[ctr*ntype+k,4]<-sl.acft[s]
          sy.out.5[ctr*ntype+k,5]<-wy.types[k]
          sy.out.5[ctr*ntype+k,6]<-success5[i,3,k]
          
          sy.out.1[ctr*ntype+k,1]<-mgmt.strats[i]
          sy.out.1[ctr*ntype+k,2]<-scls[p]
          sy.out.1[ctr*ntype+k,3]<-damt.acft[d]
          sy.out.1[ctr*ntype+k,4]<-sl.acft[s]
          sy.out.1[ctr*ntype+k,5]<-wy.types[k]
          sy.out.1[ctr*ntype+k,6]<-success1[i,3,k]
          
        }
        
      }
      
    }
    ctr<-ctr+1 # increment counter up
  }
  
  # The following section calculates the relative performance of each strategy
  #  compared to each other
  rel.perf<-array(0,dim=c(nmgmt,nreach+1,nsims*mcsims*3,ntype))
  rel.rank.perf<-rel.perf
  
  # For each reach
  for(r in 1:(nreach+1)){
    
    # For each water year type
    for(m in 1:ntype){
      
      # For each drying scenario
      for(d in 1:3){
        
        # For Angostura (r=1), Isleta (r=2), and SanAcacia (r=3)
        if(r!=4){
          rel.p<-apply(cpue.pred[,,r,d,m],MARGIN=2,FUN=function(x) (x-mean(x))/ifelse(sd(x)==0,1,sd(x)) ) # Z-score cpue across all mgmt strategies
          rank.p<-apply(cpue.pred[,,r,d,m],MARGIN=2,rank,ties.method="average" ) # rank Z-scored performances
          if(d==1){
            sr.perf<-rel.p# Store standardized performance
            rank.perf<-rank.p# Store rank performance
          } 
          else{
            sr.perf<-cbind(sr.perf,rel.p)# Store standardized performance
            rank.perf<-cbind(rank.perf,rank.p)# Store rank performance
          } 
        }
        # For MRG wide values
        if(r==4){
          rel.p<-apply(annests[,,d,m],MARGIN=2,FUN=function(x) (x-mean(x))/ifelse(sd(x)==0,1,sd(x)) )# Z-score cpue across all mgmt strategies
          rank.p<-apply(annests[,,d,m],MARGIN=2,rank,ties.method="average") # rank Z-scored performances
          
          if(d==1){
            sr.perf<-rel.p # Store standardized performance
            rank.perf<-rank.p # Store rank performance
          } 
          else{
            sr.perf<-cbind(sr.perf,rel.p) # Store standardized performance
            rank.perf<-cbind(rank.perf,rank.p) # Store rank performance
          } 
          
        }
      }
      rel.perf[,r,,m]<-sr.perf # store standardized performance
      rel.rank.perf[,r,,m]<-rank.perf # Store ranked performance
    }
    
    
  }
  
}

####################################
# The following section creates two figures
#  1. Proportional success by management strategy and drying scenario
#     for management target 1 RGSM per 100m2
#  2. Proportional success by management strategy and drying scenario
#     for downlisting target 5 RGSM per 100m2
####################################
setwd(fig.wd)
pdf("RGSM_SingleYear_PropSuccessRanks.pdf",height=6.4,width=7.3)
a<-aggregate(as.numeric(sy.out.1[,6]),list("Strat"=sy.out.1[,1],"WY"=sy.out.1[,5]),quantile,probs=c(.05,.5,.95))

a<-a[order(a$Strat),]
par(mar=c(6,5,2,2),las=2)
cols<-c("dodgerblue","darkorange","darkgreen")
plot(NA,NA,ylim=c(0.0,1),xlim=c(1,nmgmt),xaxt="n",xlab="",ylab="Probability of CPUE >1")
abline(v=seq(1,nmgmt),col="grey70",lwd=.5)
for(i in 1:3){
  subdatL<-a[a$WY==c("Low","Mid","High")[1],]
  subdat<-a[a$WY==c("Low","Mid","High")[i],]
  subdat<-subdat[order(subdatL[,3][,1]),]
  
  points(subdat[,3][,2],pch=16,col=c("dodgerblue","darkorange","darkgreen")[i])
  segments(x0=seq(1,nmgmt),y0=subdat[,3][,1],
           x1=seq(1,nmgmt),y1=subdat[,3][,3],col=c("dodgerblue","darkorange","darkgreen")[i],
           lwd=2)
}
axis(side=1,at=seq(1,nmgmt),labels=subdatL$Strat[order(subdatL[,3][,1])],las=2,cex.axis=.75)


a<-aggregate(as.numeric(sy.out.5[,6]),list("Strat"=sy.out.5[,1],"WY"=sy.out.5[,5]),quantile,probs=c(.05,.5,.95))

a<-a[order(a$Strat),]
par(mar=c(6,5,2,2),las=2)
plot(NA,NA,ylim=c(0,1),xlim=c(1,nmgmt),xaxt="n",xlab="",ylab="Probability of CPUE >5")
abline(v=seq(1,nmgmt),col="grey70",lwd=.5)
for(i in 1:3){
  subdatL<-a[a$WY==c("Low","Mid","High")[1],]
  subdat<-a[a$WY==c("Low","Mid","High")[i],]
  subdat<-subdat[order(subdatL[,3][,1]),]
  
  points(subdat[,3][,2],pch=16,col=c("dodgerblue","darkorange","darkgreen")[i])
  segments(x0=seq(1,nmgmt),y0=subdat[,3][,1],
           x1=seq(1,nmgmt),y1=subdat[,3][,3],col=c("dodgerblue","darkorange","darkgreen")[i],
           lwd=2)
}
axis(side=1,at=seq(1,nmgmt),labels=subdatL$Strat[order(subdatL[,3][,1])],las=2,cex.axis=.75)

dev.off()

#########################################
# Read in model output for each scenario
#########################################
setwd(out.wd)
s1files<-list.files()[as.numeric(substr(list.files(),30,30))==1] # List files for first storage index scenarios
s2files<-list.files()[as.numeric(substr(list.files(),30,30))==2] # List files for second storage index scenarios
s3files<-list.files()[as.numeric(substr(list.files(),30,30))==3] # List files for third storage index scenarios

fillist<-list(s1files,s2files,s3files) # Create list of all model output files

options(scipen=10000000) # Turn of scientific notation on plot axes
cols<-c("dodgerblue","darkorange","darkgreen") # Create llist of colors for figures

# The following lines identify which strategies supplement water in different seasons
#   and is used to organize data plotted in figure below
bs<-which(substr(mgmt.strats,start=nchar(mgmt.strats)-3,stop=nchar(mgmt.strats))=="Both")
sps<-which(substr(mgmt.strats,start=nchar(mgmt.strats)-5,stop=nchar(mgmt.strats))=="Spring")
sums<-which(substr(mgmt.strats,start=nchar(mgmt.strats)-5,stop=nchar(mgmt.strats))=="Summer")
sus<-which(substr(mgmt.strats,start=nchar(mgmt.strats)-5,stop=nchar(mgmt.strats))!="Spring" & 
             substr(mgmt.strats,start=nchar(mgmt.strats)-3,stop=nchar(mgmt.strats))!="Both" &
             substr(mgmt.strats,start=nchar(mgmt.strats)-5,stop=nchar(mgmt.strats))!="Summer")
seasord<-c(1,sps,bs,sums,sus[-1]) # Specifies order of strategies to be plotted


####################################
# The following section creates three figures, each 3x3 panels for
#   different combinations of storage, discretionary water, and 
#   proportional flow storage scenarios
#
# Each panel demonstrates the performance of each mgmt strategy
#  relative to No Action across all drying scenarios
####################################
setwd(fig.wd)
pdf("RGSM_Generalized_SingleYear_3x3RelativePerf.pdf",height=9,width=13)
layout(matrix(c(7,8,9,
                4,5,6,
                1,2,3),nrow=3,ncol=3,byrow=T))
par(oma=c(9,9,2,0),mar=c(.1,.1,.1,.1),las=0)

for(i in 1:length(fillist)){
  for(p in 1:length(scls)){
    for(d in 1:length(damt.acft)){
      setwd(out.wd)
      load(file=fillist[[i]][as.numeric(substr(fillist[[i]],34,34))==p & 
                               as.numeric(substr(fillist[[i]],38,38))==d])
      
      success5<-array(0, dim=c(nmgmt,5,ntype))
      
      for(m in 1:nmgmt){
        for(j in 1:3){
          for(k in 1:ntype){
            success5[m,j,k]<-sum(annests[m,,j,k]>1)/length(annests[m,,j,k])

          }
          
        }
        
      }
      
      plot(success5[seasord,3,1]/success5[1,1,1],ylim=c(.7,2.2),yaxs="i",
           xlim=c(0,nmgmt),ylab="Success Rate Relative to No Action (CPUE > 5)",
           xlab="",xaxt="n",las=2,pch=16,col=cols[1],type="n",yaxt="n")
      if(d==1) mtext(side=2,"Success Rate Relative\nto No Action",cex=.75,line=3)
      if(p==1) axis(side=1, at=seq(1,nmgmt),labels=mgmt.strats[seasord],las=2)
      if(d==1 & p!=1) axis(side=2,at=seq(.7,2.2,by=.1),
                           labels=c("",seq(.8,2.2,by=.1)),las=2)
      if(d==1 & p==1) axis(side=2,at=seq(.7,2.2,by=.1),las=2)
      
      if(d==1 & p==3) text(x=c(3,6,9.5,15),y=rep(2,4),c("Spring","Both","Summer","Partial\nSummer"),srt=90)

      abline(h=1)
      abline(v=c(1.5,4.5,7.5,11.5),col="grey")
      
      points(success5[seasord,3,1]/success5[1,3,1],col=cols[1],pch=16,cex=1.5)
      points(success5[seasord,1,1]/success5[1,1,1],col=cols[3],pch=16,cex=1.5)
      points(success5[seasord,2,1]/success5[1,2,1],col=cols[2],pch=16,cex=1.5)
      
    }
    
  }
  
  mtext(outer=T,side=1,line=6,at=c(.15,.5,.85),c("10K Discretionary","30K Discretionary","100K Discretionary"))
  mtext(outer=T,side=2,line=6,at=c(.16,.5,.84),c("Shift 10%","Shift 25%","Shift 50%"))
  mtext(outer=T,side=3,line=.5,paste("Storage = ",sl.acft[i]," acft",sep=""))
}

dev.off()

