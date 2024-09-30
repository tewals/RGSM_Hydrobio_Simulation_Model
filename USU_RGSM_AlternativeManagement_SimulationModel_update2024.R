###########################################
# RGSM Hydrologic Management Strategy
#   Simulations
#
# Coded by: Tim Walsworth
# e-mail: timothy.walsworth@usu.edu
#
###########################################

###########################################
# This code will simulate the predicted
#   response of RGSM CPUE to alternative
#   water management strategies. Users
#   upload alternative hydrographs 
#   depicting alternative timing and 
#   magnitude of flows at the Central Gage
#   in a single csv file.
#
# The code will then calculate the PC1 and
#   PC2 scores (from the PCA in Walsworth
#   and Budy 2020) for each management
#   strategy and predict RGSM CPUE in each
#   reach and across the whole MRG with
#   uncertainty. 
#
# RGSM density will be predicted from 
#   five different models presented in
#   Walsworth and Budy (2020) and its
#   Appendix B.
##########################################

##########################################
# 1. SET UP DIRECTORIES
#   
#  THESE WILL NEED TO BE SPECIFIED TO 
#   YOUR COMPUTER'S SPECIFIC FILING
#   SYSTEM!
# 
#  Replace the "..." with the directory in 
#  which you saved the downloaded files
##########################################

wd.params<-".../ModelParameters"
wd.strategies<-".../InputData"
wd.output<-".../ModelOutput"
 

# The following lines check whether required packages are installed
#  and installs them if they are not found. Packages quantreg is 
#  needed in the file containing the code for custom functions
#  used in the simulation model (RGSM_ManagementComparisonFunctions_update2024.r)
#  and rjags is used when reading parameter estimates from the
#  empirical model. If a warning or error regarding rjags is 
#  thrown while loading the parameter output, do not worry. This
#  error does not impact the ability to extract the parameter values.
#  The plotrix package is loaded at the end of this section of the code.

if(!require(quantreg)) install.packages("quantreg")
if(!require(rjags)) install.packages("rjags")
if(!require(plotrix)) install.packages("plotrix")
library(plotrix)

##########################################
# 2. Read in Baseline and Alternative 
#    Hydrographs 
##########################################
setwd(wd.strategies)                                              # Set working directory holding alternative management scenarios
strats<-read.csv("RGSM_AlternativeHydrographInput_AllowKnownDrying.csv",header=T)  # Load alternative management strategies

nstrats<-ncol(strats)-4    # How many strategies are being compared? (Includes baseline)

###############################################################
# Specify the format that you have provided the alternative hydrographs
#  being compared to the baseline hydrograph. This should be one of 
#  three possible values:
#  - "cfs" if the alternative hydrographs are provided as daily cfs values.
#  
#  - "differences" if the alternative hydrographs are provided as the
#                  difference in discharge in cfs from the baseline.
#  
#  - "proportion" if the alternative hydrographs are provided as proportions
#                 of the baseline hydrograph (e.g., an alternative hydrograph
#                 with daily discharge 10% higher than the baseline would have
#                 a value of 1.1).
###############################################################
alternative.format<-"differences"

###############################################################
# If alternative strategies are supplied as hydrographs with
#  proposed discharges by day, then no conversions are needed.
#  However, if alternative strategies are supplied as proposed 
#  changes (as differences from the baseline or proportions of
#  the baseline), they need to be converted into daily discharge
#  values.
#################################################################

# If alternative strategies are supplied as proportions of the
#  baseline, 
if(alternative.format=="proportion"){
  for(i in 1:(nstrats-1)){
    strats[-c(1,2),5+i]<-strats[-c(1,2),5]*strats[-c(1,2),5+i]
  }
}



# If alternative strategies are supplied as differences from the
#  baseline in cfs, the following lines of code will convert them to 
#  the appropriate format for the model.
if(alternative.format=="differences"){
  for(i in 1:(nstrats-1)){
    strats[-c(1,2),5+i]<-strats[-c(1,2),5]+strats[-c(1,2),5+i]
  }
}


# Check for negative discharge values and set to a specified
#  minimum value.
flowminset<-5 # Specify the minimum flow value allowed
for(i in 1:nstrats){
  for(j in 3:367){
    if(strats[j,4+i]<0) strats[j,4+i]<-flowminset # Set any negative values to
                                                  # the minimum value
  }
}

#################################################
# 3. Read in PCA Loadings from Walsworth and Budy 2020 report
#    to calculate integrated flow metrics, and load required 
#    custom functions.
#################################################

setwd(wd.params)    # Set working directory to location of parameter estimates
                    #  from Walsworth and Budy 2020

load("RGSM_PCAResults_thru23.rdata") # Read in PCA parameter estimates
source("RGSM_ManagementComparisonFunctions_update2024.r")  # Load custom functions

########################################################
# 5. Set additional simulation parameters
########################################################

# Set required parameters
yearpred<-2024 # For which year are you predicting CPUE?
lastyear<-2023 # What is the last year of data observed in the empirical model?
ndry<-5        # How many levels of drying extent are you examining (5 is default)
nmods<-1       # How many empirical models are you simulating from?

setwd(wd.params) # Change working directory
mod.files<-list.files()[1:nmods] # Save a list of model parameter file names. Model parameter files
                                 # begin with a number, and will be the first listed in the wd.params folder.
                                 # Thus, by selecting the top nmods files, you will only be selecitng the model files.

nreach<-3      # Set number of sampling reaches (should be 3 - Angostura, Isleta, San Acacia)
nsim<-100      # How many stochastic simulations should be run per set of parameter values?
threshs<-seq(0.05,5,by=.05) # Create a vector of CPUE thresholds for calculating 
                            #   probability of meeting different management targets.
nt<-length(threshs) # How many CPUE thresholds are being examined?

#################################################
# The following lines create storage arrays for predicted October cpue in each
# reach (cpue.pred), October CPUE at the MRG wide scale (annests), and the 
# probability of exceeding different CPUE thresholds (probexc).
#
# The dimensions of the cpue.pred array are used to index:
# 1: Management strategy
# 2: Stochastic simulation
# 3: Reach
# 4. Drying severity scenario
# 5. Model (currently only 1 implemented in simulation)
cpue.pred<-array(0,dim=c(nstrats,nsim*3000,nreach,ndry,nmods))
#
# The dimensions of the annests array are used to index:
# 1: Management strategy
# 2: Stochastic simulation
# 3. Drying severity scenario
# 4. Model (currently only 1 implemented in simulation) 
annests<-array(0,dim=c(nstrats,nsim*3000,ndry,nmods))
#
# The dimensions of the annests array are used to index:
# 1: Management strategy
# 2: Drying severity scenario
# 3: Model (currently only 1 implemented in simulation) 
# 4: Reach [indices 1-3] or MRG wide [index 4]
# 5: CPUE Threshold
probexc<-array(0,dim=c(nstrats,ndry,nmods,nreach+1,nt))
#################################################

###############################################
# 5. Run Simulations. The following code will
#    loop through different models (currently
#    using only 1 model), different management
#    strategies, different parameter values from
#    model posterior distribution estimated in
#    Walsworth and Budy (2021),
#    different stochastic simulations, and 
#    different reaches within the MRG. For each
#    iteration, it will sample 10 sites within
#    each reach given the estimated PC scores
#    for the management strategy's hydrologic
#    index (PC score), calculate the reach-
#    specific CPUE and the MRG mean CPUE.
#    Outputs are stored in the arrays designated
#    above for later analysis and storage.
################################################

start.time<-Sys.time() # Start a timer. Not necessary to run the model, but can
                       #  be useful for timing model runs.


for(m in 1:nmods){   # For each empirical model (Only 1 currently)...

  setwd(wd.params)   # Set working directory to load parameter 
                     #   values from model posterior distributions

  load(mod.files[m]) # Load model parameter values from the working directory
 
  sims.out<-mod.out$BUGSoutput$sims.list  # Specify shorter identifier name for parameter values
  
  nmc<-nrow(sims.out$alpha.c) # Store number of mcmc samples for each parameter
  
  for(s in 1:nstrats){    # For each alternative hydrograph (management strategy)
                          #   under consideration...
  
    print(paste("Running Strategy ",s," out of ",nstrats,sep=""))  # Print which management strategy (by number) is currently being simulated
    
    strat.data<-strats[-c(1,2),c(1:4,4+s)]   # Subset hydrograph data for only the management strategy being simulated
    strat.dry<-strats[c(1,2),c(4+s)]     # Extract drying extent indicators (-1 for unknown)
    
    colnames(strat.data)[5]<-"discharge" # Change column name of management strategy to "discharge" 
                                         #   (necessary for the function calculating PC scores)
    
    pcs<-calc.pcscores(daily.avs=strat.data,drying=strat.dry,pca1=pca1) #Calculate PC scores for hydrograph
    
    set.seed(100)      # Set random number generator seed - ensures that all management strategies are
                       #   compared against the same random walk conditions
    for(k in 1:length(sims.out$alpha.c[,1])){ # For each parameter value in the mcmc posterior distribution
      
      for(j in 1:nsim){   # For each stochastic simulation iteration
        
        fldsim<-pcs[,1]   # Set hydrologic index calculated from the pca
        
        tsuse<-sims.out$ts[k,25]+sum(rnorm(yearpred-lastyear)) # generate random walk values for the latent trend in the hurdle component of the simulation model
        for(i in 1:3){  # For each reach...
          
          pred.ca<-sims.out$alpha.c[k,i]*exp(-sims.out$beta.off[k]*exp(-sims.out$beta.denspred[k]*fldsim)) # Calculate expected site-level CPUE
          pred.ca[pred.ca==0]<-.00000000001 # Values cannot be zero, so add small value if zero
          predobs<-sims.out$alpha.p[k,i]+sims.out$beta.prespred[k]*fldsim+tsuse # Calculate the (logit) probability of RGSM being present at a sample site
          predobs<-exp(predobs)/(1+exp(predobs)) # Transform logit probability of presence to probability of presence
          prate<-(pred.ca+sqrt(pred.ca^2+4*(sims.out$cv[k]*pred.ca)^2))/(2*(sims.out$cv[k]*pred.ca)^2) # Calculate gamma distribution rate parameter
          psh<-1+pred.ca*prate              # Calculate gamma distribution shape parameter
          predfish<-c()                     # Create storage vector for predicted CPUE under different drying extent conditions
          predfish[1]<-mean(rbinom(10,1,prob=predobs[1])*rgamma(10,shape=psh[1],rate=prate[1])) # Predict CPUE under low drying extent
          predfish[2]<-mean(rbinom(10,1,prob=predobs[2])*rgamma(10,shape=psh[2],rate=prate[2])) # Predict CPUE under 25th percentile drying extent
          predfish[3]<-mean(rbinom(10,1,prob=predobs[3])*rgamma(10,shape=psh[3],rate=prate[3])) # Predict CPUE under median drying extent
          predfish[4]<-mean(rbinom(10,1,prob=predobs[4])*rgamma(10,shape=psh[4],rate=prate[4])) # Predict CPUE under 75th percentile drying extent
          predfish[5]<-mean(rbinom(10,1,prob=predobs[5])*rgamma(10,shape=psh[5],rate=prate[5])) # Predict CPUE under high drying extent
          cpue.pred[s,((k-1)*nsim+j),i,,m]<-predfish # Store predicted CPUE output
          
        } # end "reach" loop
        annests[s,((k-1)*nsim+j),,m]<-colMeans(cpue.pred[s,((k-1)*nsim+j),,,m]) # Calculate MRG mean CPUE
        
      } # end stochastic simulation loop
      
    } # end posterior distribution parameter value loop
    
    for(z in 1:nt){  # for each CPUE threshold
      probexc[s,,m,4,z]<-colSums(annests[s,,,m]>threshs[z])/length(annests[s,,1,m]) # Calculate the proportion of simulations in which the MRG mean CPUE exceeded the threshold
      probexc[s,,m,1,z]<-colSums(cpue.pred[s,,1,,m]>threshs[z])/length(annests[s,,1,m]) # Calculate the proportion of simulations in which the Angostura CPUE exceeded the threshold
      probexc[s,,m,2,z]<-colSums(cpue.pred[s,,2,,m]>threshs[z])/length(annests[s,,1,m]) # Calculate the proportion of simulations in which the Isleta CPUE exceeded the threshold
      probexc[s,,m,3,z]<-colSums(cpue.pred[s,,3,,m]>threshs[z])/length(annests[s,,1,m]) # Calculate the proportion of simulations in which the San Acacia CPUE exceeded the threshold
    } # end threshold loop
    
    
  } # end management strategy loop
  
  rm(mod.out) # remove model parameter values from R storage
  rm(sims.out) # remove model parameter values from R storage
  
} # end model loop

print(Sys.time()-start.time) # Print out how long full set of simulations took to run

##########################################################
# 6. Calculate relative performance of different strategies
#     within each stochastic iteration. 
##########################################################

rel.perf<-array(0,dim=c(nstrats,nreach+1,nsim*nmc*nmods*ndry)) # Create storage for relative performance output
                                                                # This will hold the z-scored CPUE for each management 
                                                                # strategy relative to all other management strategies
                                                                # for each stochastic simulation. At the end, you will
                                                                # have a vector of relative performance values for each
                                                                # strategy. The first dimension of the array indicates
                                                                # which strategy has its performance stored there. The
                                                                # second dimension indicates which reach the performance
                                                                # is being measured for (1 = Angostura, 2 = Isleta, 3 =
                                                                # San Acacia, 4 = MRG-mean). The third dimension is for 
                                                                # Each stocahstic iteration across all underlying models
                                                                # parameter values from the posterior distribution, and
                                                                # drying extent values.

for(r in 1:(nreach+1)){  # for each reach AND the MRG as a whole

      for(m in 1:nmods){  # for each underlying model (only 1 as default)
        for(d in 1:ndry){ # for each drying scenario
          if(r!=4){  # If r does not equal 4, run the following (r=1,2,3 is for individual reaches, r=4 is the MRG as a whole)
            rel.p<-apply(cpue.pred[,,r,d,m],MARGIN=2,FUN=function(x) (x-mean(x))/ifelse(sd(x)==0,1,sd(x)) ) # Calculate the performance of each management strategy relative to the others for each stochastic simulation
            if(m==1 & d==1) sr.perf<-rel.p  # If this is the first model and the first drying scenario, create storage vector
            else sr.perf<-cbind(sr.perf,rel.p) # Otherwise, bind columns together
          }
            if(r==4){ # If r equals 4, run this section (for MRG mean values)
              rel.p<-apply(annests[,,d,m],MARGIN=2,FUN=function(x) (x-mean(x))/ifelse(sd(x)==0,1,sd(x)) ) # Calculate the performance of each management strategy relative to the others
              if(m==1 & d==1) sr.perf<-rel.p # If this is the first model and the first drying scenario, create storage vector
              else sr.perf<-cbind(sr.perf,rel.p) # Otherwise, bind columns together
            }
        }
      }
  
    rel.perf[,r,]<-sr.perf # Store relative performance in array. 
}

setwd(wd.output)
save.image("SimulationOutput.RData",compress=T) # Save simulation model output as .rdata file.
############################################
# 7. Summary Figures
#  - NOTE: If you are comparing more than 14
#    hydrographs, you will need to set up a
#    new, custom layout to get all figures
#    on a single page. Currently can handle
#    predictions for 14 hydrographs on one
#    page.
############################################

# Set global plotting parameters
ybrks<-seq(0,1,by=.1) # set y-axis breaks to be used when determining the y-axis limits to use
ylimuse<-which.min((min(probexc)-ybrks)[(min(probexc)-ybrks)>=0]) # Determine y-axis minimum limit
threshprobs<-signif(probexc[,3,1,4,c(6,20,100)],3) # Specify which CPUE threshold values you want to print the exceedance probabilities for
tabspot<-ifelse(min(signif(probexc[,3,1,4,c(70)],3))>.4,"bottomleft","topright") # Determine where to place exceedance probability table within plots
rownames(threshprobs)<-colnames(strats[c(5:(4+nstrats))]) # Set the strategy names to be printed in tables
colnames(threshprobs)<-paste("CPE >",threshs[c(6,20,100)]) # Set the CPUE threshold values to be printed as column headers in tables

# The following lines will generate a PDF containing all of the summary figures
setwd(wd.output)  # set output directory
strats.droprows<-strats[-c(1,2),]

#################################################
# If you want to produce plot within R or Rstudio
#   Do not tun the next line. It will open a PDF
#   as the plotting window.
#################################################
pdf("Simulation_Figures.pdf",height=6,width=9) # create and begin editing figure PDF

layout(plotlayouts.list[[nstrats]]) # Set plotting layout
par(mar=c(3,5,1,1),las=1,cex.text=1.5,cex.axis=1) # Set plot margins, axis orientation, text size, and axis text size
cols<-c("black","darkorange","dodgerblue","darkgreen") # Set colors. Use this if only examining 4 strategies, otherwise, run next line
#cols<-c("black","grey",brewer.pal(12,"Paired")) # Use this if examining between 5 and 14 strategies
mos<-format(ISOdatetime(2000,1:12,1,0,0,0),"%b") # Create vector of Month names for axis labels

for(i in 1:nstrats){ #For each management strategy (hydrograph)
  if(i ==1){ # If it is the first hydrograph, run the following lines
    plot(strats.droprows[,4+i],type="l",col=cols[i],ylab="Discharge",
         xlab="",xaxt="n")                                    # Plot Discharge by day
    axis(side=1,at=which(strats.droprows$day==5),
         labels=unique(mos[strats.droprows$month]),las=2) # Add x-axis labels
  } # end code specific to first hydrograph
  if(i >1){ # If it is not the first hydrograph, run the following lines
    plot(strats.droprows[,4+i]-strats.droprows[,4+1],col=cols[i],type="l",
         ylab="Discharge\nDifference",xaxt="n",xlab="",
         ylim=c(1.05*min(strats.droprows[,6:ncol(strats)]-strats.droprows[,5]),1.05*max(strats.droprows[,6:ncol(strats.droprows)]-strats.droprows[,5]))) # Plot the daily difference in discharge between the alternative hydrograph and the baseline hydrograph
    axis(side=1,at=which(strats.droprows$day==5),
         labels=unique(mos[strats.droprows$month]),las=2)  # Add x-axis labels
    abline(h=0,lty=2) # Add horizontal line at y=0
    polygon(c(strats.droprows$doy,rev(strats.droprows$doy)),
            c(rep(0,length(strats.droprows$doy)),rev(strats.droprows[,4+i]-strats.droprows[,4+1])),
            col=cols[i],border=cols[i])  # Plot a filled in polygon demonstrating the discharge difference between the alternative and baseline hydrographs
  } # end code specific to alternative hydrograph
  legend("topleft",bty="n",colnames(strats.droprows[4+i])) # Add text indicating alternative management strategy name
} # End discharge plots section


par(mar=c(5,5,1,1),las=1)  # Set margins and axis orientation
for(r in 4:(nreach+1)){ # For the MRG-mean values only (to run for each reach, set the "4" equal to 1)
  plot(NA,NA,ylim=c(ybrks[ylimuse],1),xlim=c(0,max(threshs)),yaxs="i",
       las=1,
       ylab="Exceedance Probability",xlab="CPUE Threshold") # Create blank plotting space for probability of exceeding CPUE thresholds
  
  for(m in 1:nmods){ # For each underlying model (only 1 as default)
    for(s in 1:nstrats){ # For each strategy...
      for(d in 3:3){ # For median drying scenario [ set to 1 for low drying, 5 for high drying]
        
        lines(threshs,probexc[s,d,m,r,],lwd=2,lty=1,
              col=cols[s]) # Plot lines of probability of exceeding CPUE threshold
      }
    }
  }
}
addtable2plot(tabspot,table=as.data.frame(threshprobs),
              display.rownames = T,bg="white",vlines=T,xpad=.2,cex=.95) # Add table with exceedance values for specified thresholds to plotting region

box() # Add box around plot
par(mar=c(12,5,1,1),las=1) # set figure margins and axis orientation
for(i in 4:(nreach+1)){ # For MRG mean values... (to run for each reach, set the "4" equal to 1)
  plot(NA,NA,ylim=c(-1.05*max(rel.perf,na.rm=T),1.05*max(rel.perf,na.rm=T)),xlim=c(.5,nstrats+.5),ylab="Relative\nPerformance",
       xlab="",xaxt="n",yaxt="n") # Create blank plotting space for relative performance figure
  axis(side=1,at=seq(1,nstrats),cex.axis=1,
       labels=colnames(strats[5:(4+nstrats)]),las=2) # Add x-axis with strategy names
  abline(h=0,lty=2) # Add horizontal line at y=0
  
  for(j in 1:nstrats){ # For each strategy
    ranks<-quantile(rel.perf[j,i,],probs=c(.025,.05,.1,.25,.45,.5,
                                           .55,.75,.9,.95,.975),na.rm=T) # Calculate relative performance quantiles
    segments(x0=j,y0=ranks[1],x1=j,y1=ranks[11],col=cols[j]) # draw line for 95% CI
    segments(x0=j,y0=ranks[2],x1=j,y1=ranks[10],lwd=3,col=cols[j]) # Draw line for 90% CI
    segments(x0=j,y0=ranks[3],x1=j,y1=ranks[9],lwd=5,col=cols[j]) # Draw line for 80% CI
    segments(x0=j,y0=ranks[4],x1=j,y1=ranks[8],lwd=7,col=cols[j]) # Draw line for 50% CI
    segments(x0=j,y0=ranks[5],x1=j,y1=ranks[7],lwd=13,col=cols[j]) # Draw line for 10% CI
  }
  
}

text((par()$usr[2]-par()$usr[1])*.025+(par()$usr[1]),
     (par()$usr[4]-par()$usr[3])*.9+(par()$usr[3]),"Better",pos=4) # Add text to plot indicating region of better performing strategies
text((par()$usr[2]-par()$usr[1])*.025+(par()$usr[1]),
     (par()$usr[4]-par()$usr[3])*.1+(par()$usr[3]),"Worse",pos=4) # Add text to plot indicating region of more poorly performing strategies

dev.off() # Close and save PDF of figures. You need to run this line before you can open and view your figures.


###################################
# The following section will make a figure demonstrating the predicted probability
#  of meeting a specified CPUE threshold (default is 1) in each individual reach
#  across all management strategies examined (see Figure 17 in report for example)
###################################
cpuethresh<-1 # What CPUE threshold are you interested in?
thresh1<-which(threshs==cpuethresh)


setwd(wd.output)
#################################################
# If you want to produce plot within R or Rstudio
#   Do not tun the next line. It will open a PDF
#   as the plotting window.
#################################################
# create and begin editing figure PDF
pdf("SimulationModel_Figures_ReachSpecificProbabilities.pdf",height=8,width=9)

# The following line sets up the plot layout as a function of the number of
#  alternative hydrographs being examined
layout(matrix(seq(1:(ceiling(sqrt(nstrats+1))^2)),nrow=ceiling(sqrt(nstrats+1)),
              ncol=ceiling(sqrt(nstrats+1)),byrow=F),
       heights=c(rep(1,9),1.75),widths=c(1,rep(1,6)))

rowsplt<-ceiling(sqrt(nstrats+1)) # Calculation for plotting axes

par(mar=c(.5,.5,.5,.5),oma=c(6,6,0,1),xpd=F) # Set plotting parameters

# The following section plots the expected performance for the baseline hydrograph
plot(NA,ylim=c(0,1),xlim=c(0.5,4.5),xaxt="n",xlab="Reach",
     ylab=paste("P(CPUE > ",cpuethresh,")",sep=""),
     las=1,yaxs="i",yaxt="n")
abline(h=c(.25,.5,.75),lty=2,col="grey") # add horizontal lines at .25, .5, and .75
axis(side=2,at=c(0,1),las=1,line=.25,hadj=0)
#axis(side=1,at=seq(1,4),labels=c("Angostura","Isleta","San Acacia","MRG"),
#    las=0,cex.axis=1.5) # Add reach names to x-axis
text(.75,.9,labels="1",col="blue") # Add strategy ID number to panel
text(.07,.5,paste("P(CPUE > ",cpuethresh,")",sep=""),srt=90,xpd=NA)
# The following section loops through each reach and the MRG wide results
# and plots the probability of exceeding the threshold specified above
for(j in 1:4){
  for(k in 1:3){
    points(j,probexc[1,c(1,3,5)[k],1,j,thresh1],col=cols[k],pch=16,cex=1.5)
  }
  
}

# The following section loops through each alternative management strategy
#   and plots the probability of meeting the management target 
for(i in 2:nstrats){
  
  par(mar=c(.5,.5,.5,.5)) #Set plotting parameters
  if(i%in%c(rowsplt*seq(1,rowsplt))) par(mar=c(3,.5,.5,.5)) # set plotting parameters
  # The following section plots the expected performance for hydrograph i
  plot(NA,ylim=c(0,1),xlim=c(.5,4.5),xaxt="n",yaxt="n",yaxs="i")
  abline(h=c(.25,.5,.75),lty=2,col="grey") # add horizontal lines
  # The following section loops through each reach and the MRG wide results
  # and plots the probability of exceeding the threshold specified above
  for(j in 1:4){
    for(k in 1:3){
      points(j,probexc[i,c(1,3,5)[k],1,j,thresh1],col=cols[k],pch=16,cex=1.5)
    }
    
  }
  
  text(.75,.9,labels=i,col="blue") # Add strategy identifier
  
  # The following lines add x and y-axes to the left and bottom panels
  if(i%in%c(rowsplt*seq(1,rowsplt))) axis(side=1,at=seq(1,4),labels=c("Angostura","Isleta","San Acacia","MRG"),las=2,
                                          cex.axis=1.25)
  oddev<-rowsplt%%2
  if(i%in% c(seq(1,ifelse(oddev,rowsplt,rowsplt-1)))){
    axis(side=2,at=c(0,1),las=1,hadj=0,line=1.55)
    text(0.07,.5,paste("P(CPUE > ",cpuethresh,")",sep=""),srt=90,xpd=NA)
  } 
  if(i%in% c(seq(2,ifelse(!oddev,rowsplt,rowsplt-1)))){
    axis(side=2,at=c(0,1),las=1,line=.25,hadj=0)
    text(0.07,.5,paste("P(CPUE > ",cpuethresh,")",sep=""),srt=90,xpd=NA)
  } 
  
}

# The following lines add a legend
par(mar=c(.5,.5,.5,.5))
frame()
legend("center",pch=16,col=cols,bty="n",
       legend=c("Low drying","Median drying","High drying"))

#mtext(side=2,outer=T,line=3.5,paste("P(CPUE > ",cpuethresh,")",sep="")) # Add y-axis label

dev.off() # Close and save PDF of figures. You need to run this line before you can open and view your figures.
