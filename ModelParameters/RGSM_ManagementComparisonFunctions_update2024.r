#########################################
# Functions to be used within the script
#  comparing the performance of
#  alternative water management
#  strategies. This file is called from
#  within the main script.
#
# Coded by: Tim Walsworth
# e-mail: timothy.walsworth@usu.edu
#########################################

#########################################
# Load required libraries and data sets
#########################################
library(quantreg) # Loads library used for quantile regressions predicting
                  # summer dying extent

pcadat<-read.csv("MRG_HydroPCA_thru23.csv",header=T)# Loads hydrologic data
                                                    # required to fit quantile
                                                    # regressions for predicting
                                                    # summer channel drying

#########################################################
# Function: calc.pcscores
#
# - This function reads in a hydrograph and PCA
#   output to convert the hydrograph into PC scores
#   used by the simulation model to predict RGSM
#   densities.
#
# - Inputs:
#   - daily.avs: Matrix or data frame including
#                the following columns:
#                - year = Year
#                - doy = Julian day of year
#                - day = day of month
#                - month = integer indicating month
#                - discharge = daily average discharge
#                              at the ABQ Central Gage
#
#   - pca1: PCA Model Object. For the RGSM simulation
#           model, this should be set to equal pca1, 
#           which is read into the workspace in the 
#           main simulation file.
#
# - Output: a 5x2 matrix of PC scores for the 
#            hydrograph provided. The first column
#            is PC1 (used in the main simulation
#            model). The second column is PC2 (used
#            in some alternative model structures
#            reported in the Walsworth and Budy 2020
#            report to USBR, but which is not
#            yet integrated into the simulation model).
#            The different rows represent the PC scores
#            under different drying severity scenarios.
#            Row 1 represents the PC scores under 10th 
#            percentile drying, row 2 is 25th percentile,
#            row 3 is 50th percentile, row 4 is 75th
#            percentile, and row 5 is 90th percentile.
#
############################################################
calc.pcscores<-function(daily.avs,drying,pca1){

  abqmin<-min(daily.avs$discharge[daily.avs$doy>150 & daily.avs$doy<300 ]) # Finds the minimum daily 
                                                      # discharge after May
  
  abqspr<-daily.avs[daily.avs$month%in%c(5,6),] # Subsets out May and June, which
                                                #  are used to calculate spring
                                                #  high flow conditions.
  
  # The following line calculates mean discharge for each day during May and June
  bd<-aggregate(abqspr$discharge,list(abqspr$year,abqspr$month,abqspr$day),mean)
  
  by<-aggregate(bd$x,list(bd$Group.1),mean) # Calculates mean daily discharge across
                                            # May and June
  
  by$vol<-by$x*60*60*24*2.29569e-5*61       # Calculates total volume of flow
                                            # in May and June (cfs x 60s x 60min
                                            # x 24hr x 61days to acft)
  
  # The following four (4) lines calculate the number of days of spring high
  #  flows above the following thresholds: 1500, 2000, 2500, and 3000 cfs
  bt15<-aggregate(bd$x,list(bd$Group.1),function(x) sum(x>1500))
  bt20<-aggregate(bd$x,list(bd$Group.1),function(x) sum(x>2000))
  bt25<-aggregate(bd$x,list(bd$Group.1),function(x) sum(x>2500))
  bt30<-aggregate(bd$x,list(bd$Group.1),function(x) sum(x>3000))
  
  # The following code calculates the timing of the spring high flow peak
  #   by finding the peak 5-day moving average peak.
  mov.5<-rep(0,length(daily.avs[,1])) # Creates empty storage vector
  for(i in 5:length(daily.avs[,1])){ # Loops through each day of the year
    mov.5[i]<-mean(daily.avs$discharge[((i-4):i)]) # Calculates 5-day avg. discharge
  } # close loop through days
  
  doy<-seq(1,length(daily.avs$discharge)) # Create a vector of day of year
  
  # The following line identifies the day of year with the max 5-day avg. discharge
  pk<-doy[which(mov.5[1:220]==max(mov.5[1:220],na.rm = T))] 
  if(length(pk)>1) pk<-pk[round(length(pk)/2)] # If two identical peaks, set peak to midpoint
  
  if(drying[1]<0){
    # The following lines calculate the quantile regression parameters for predicting
    #  summer channel drying extent in San Acacia and Isleta from observed summer
    #  minimum discharge values. The 10th, 25th, 50th, 75th, and 90th
    #  percentiles are predicted
    cSA<-coef(rq(log(SanAca_MileDayDry+100)~CentralMin,data=pcadat,tau=c(.1,.25,.5,.75,.9)))
    cIS<-coef(rq(log(Isleta_MileDayDry+100)~CentralMin,data=pcadat,tau=c(.1,.25,.5,.75,.9)))
    
    
    # The following lines take the parameter estimates from the quantile regression
    #   above and use them to estimate expected summer channel drying for the
    #   year of interest across range of drying intensities (percentiles)
    SAdry<-exp(cSA[1,]+abqmin*cSA[2,])-100
    ISdry<-exp(cIS[1,]+abqmin*cIS[2,])-100
    
    SAdry[SAdry<0]<-0 # Checks for negative drying extent and sets to 0
    ISdry[ISdry<0]<-0 # Checks for negative drying extent and sets to 0
    
  }
  
  if(drying[1]>=0)
  {
    SAdry<-rep(drying[2],5)
    ISdry<-rep(drying[1],5)
  }
  
  # The following lines compile the calculated hydrologic metrics for the hydrograph
  #  of interest to be converted  to PC scores, and centers and scales all values
  #  relative to the data used to fit the PCA in Walsworth and Budy 2020.
  # This is repeated for each of the drying severity scenarios examined.
  pcadat1<-c(by[1,3],abqmin,by[1,2],bt15[1,2],bt20[1,2],bt25[1,2],bt30[1,2],SAdry[1],ISdry[1],pk)
  pcadat1.cs<-(pcadat1-pca1$center)/pca1$scale # Center and scale all values
  
  pcadat2<-c(by[1,3],abqmin,by[1,2],bt15[1,2],bt20[1,2],bt25[1,2],bt30[1,2],SAdry[2],ISdry[2],pk)
  pcadat2.cs<-(pcadat2-pca1$center)/pca1$scale # Center and scale all values
  
  pcadat3<-c(by[1,3],abqmin,by[1,2],bt15[1,2],bt20[1,2],bt25[1,2],bt30[1,2],SAdry[3],ISdry[3],pk)
  pcadat3.cs<-(pcadat3-pca1$center)/pca1$scale # Center and scale all values
  
  pcadat4<-c(by[1,3],abqmin,by[1,2],bt15[1,2],bt20[1,2],bt25[1,2],bt30[1,2],SAdry[4],ISdry[4],pk)
  pcadat4.cs<-(pcadat4-pca1$center)/pca1$scale # Center and scale all values
  
  pcadat5<-c(by[1,3],abqmin,by[1,2],bt15[1,2],bt20[1,2],bt25[1,2],bt30[1,2],SAdry[5],ISdry[5],pk)
  pcadat5.cs<-(pcadat5-pca1$center)/pca1$scale # Center and scale all values
  
  # The following lines multiply the centered and scaled hydrologic metrics by
  #  the PCA loadings to calculate PC1 scores
  pc1.1<-sum(pcadat1.cs*pca1$rotation[,1])
  pc1.2<-sum(pcadat2.cs*pca1$rotation[,1])
  pc1.3<-sum(pcadat3.cs*pca1$rotation[,1])
  pc1.4<-sum(pcadat4.cs*pca1$rotation[,1])
  pc1.5<-sum(pcadat5.cs*pca1$rotation[,1])
  
  # The following lines multiply the centered and scaled hydrologic metrics by
  #  the PCA loadings to calculate PC2 scores
  pc2.1<-sum(pcadat1.cs*pca1$rotation[,2])
  pc2.2<-sum(pcadat2.cs*pca1$rotation[,2])
  pc2.3<-sum(pcadat3.cs*pca1$rotation[,2])
  pc2.4<-sum(pcadat4.cs*pca1$rotation[,2])
  pc2.5<-sum(pcadat5.cs*pca1$rotation[,2])
  
  # Compiles the calculated PC scores into a 5x2 matrix
  pcs.out<-cbind(c(pc1.1,pc1.2,pc1.3,pc1.4,pc1.5),
                 c(pc2.1,pc2.2,pc2.3,pc2.4,pc2.5))
  colnames(pcs.out)<-c("PC1","PC2") # Names columns of output matrix
  
  return(pcs.out) # Returns output matrix containing PC scores
  
}


######################################################################
# Function: convert.hydro.mgmt
#
# This function applies generalized flow management strategies
#  to a hydrograph, returning the managed hydrograph.
#
# - Inputs:
#   - hydro - starting hydrograph to which management is
#             applied
#   - strat - generalized management strategy being applied.
#             Management strategy should be formatted as 
#             a vector of length 12, each value representing 
#             the action in each month. Values should be one 
#             of (-1, 0, 1), with -1 meaning flow is stored
#             during this month, 0 meaning no change to that
#             month, and 1 meaning flow is added to the month.
#             Flow is added evenly across all days to which
#             mgmt specifies for flow additions.
#   - scale - Specifies the proportion of daily flows to be
#             stored during months specified for storage.
#   - discretionary - specifies whether discretionary water is
#                     available (0 means no, any other value means yes)
#   - discamt - specifies the volume of discretionary water
#               available
#   - type - specifies the type of management strategy being
#            applied. Currently only supports proportional
#            shifts in flow (type = "prop")
#   - storage - Logical (T/F) - should end of year water 
#               volume stored be returned.
#   - storage.limit - limit to total storage.
#   - stored.in - volume of stored water transferred from previous year
#   - dbg - Logical (T/F) - run in debug mode? Debug mode returns 
#           text checkpoints throughout the script so errors can
#           be located.
#
# - Output:
#   - mod.hyd - modified hydrograph
#   - total.stored - remaining water in storage.
#
########################################################################
convert.hydro.mgmt<-function(hydro=test.h, strat=mgmt, scale=scl,
                             discretionary=0, discamt = 0, 
                             type=c("prop"),storage=F,storage.limit=NULL,
                             stored.in=0,dbg=F){

  if(dbg) print("Debug Checkpoint A") # Debugging checkpoint
  
  ###########################################################################
  # The following section checks how many management periods are specified
  #   and specifies which days of year should be modified by each
  #   management period.
  
  ncells<-length(strat) # how many management periods?
  if(ncells==3)  mod.cells<-list(c(60:120),c(121:181),c(213:273))
  if(ncells==5) mod.cells<-list(c(60:90),c(91:120),c(121:151),c(152:181),
                                c(182:212),c(213:242),c(243:273))
  if(ncells==12) mod.cells<-list(c(1,31),c(32:59),c(60:90),c(91:120),c(121:151),
                                 c(152:181),c(182:212),c(213:243),c(244:273),
                                 c(274:304),c(305:334),c(335:365))
  ###########################################################################
  
  total.stored<-stored.in # Add carry-over storage to total stored water volume
  mod.hyd<-hydro # Set mod.hyd equal to the hydrograph of daily mean flows being
                 #  managed
  
  ################################################
  # The following section checks whether any management actions are specified.
  #  If not, then it returns the unmodified hydrograph and (if storage = T) the
  #  volume of stored water to be carried over to the next year.
  if(sum(abs(strat))==0){
    if(storage) return(list(mod.hyd,total.stored)) 
    if(!storage) return(mod.hyd)
  } 
  ################################################
  
  if(dbg) print("Debug Checkpoint B") # Debugging checkpoint
  
  ################################################
  # The following sections sets the days in which water is being stored
  #  and the days in which water is being supplemented by management actions.
  store.cells<-unlist(mod.cells[which(strat<0)]) # Set storage days
  add.cells<-unlist(mod.cells[which(strat>0)])   # Set supplementation days
  ################################################
  
  if(dbg) print("Debug Checkpoint C") # Debugging checkpoint
  
  ###############################################################
  # The following section of code is only run if type = "prop"
  ###############################################################
  if(type=="prop"){
    if(dbg) print("Debug Checkpoint D") # Debugging checkpoint
    
    # The following section of code will only be run if there are storage
    #   periods in the specified management strategy
    if(length(store.cells)>0){
      for(s in 1:length(store.cells)){ # For each day in which water is to be stored
        
        # The following section of code will only be run if there is storage 
        #   capacity remaining.
        if(total.stored<storage.limit){ 
          # The following line calculates how much water is stored from the current day.
          #  If the specified proportion of flow is greater than the amount of available
          #  storage, only enough water to fill the storage is stored.
          store.amt<-min(mod.hyd[store.cells[s]]*scale,storage.limit-total.stored)
          
          # The following line removes the stored amount of flow from the current day.
          mod.hyd[store.cells[s]]<-mod.hyd[store.cells[s]]-store.amt
          total.stored<- total.stored+store.amt # Add the volume stored to total
        } 
      }
    }

    if(dbg) print("Debug Checkpoint E") # Debugging checkpoint
    
    # The following line adds discretionary water to the stored volume
    if(discretionary!=0) total.stored<-total.stored+discamt

    # The following section of code will only be run if there are supplementation
    #   periods in the specified management strategy
    if(length(add.cells)>0){
      add.amnt<-total.stored/length(add.cells) # Divide total stored water by
                                               #  number of days in supplementation
                                               #  period
      
      mod.hyd[add.cells]<-hydro[add.cells]+add.amnt # Add supplementation flows
                                                    #  to each day of supplementation
      
      # Remove water from storage
      total.stored<-total.stored-add.amnt*length(add.cells)
    } 
    
    if(dbg) print("Debug Checkpoint F") # Debugging checkpoint
    if(dbg) print(length(add.cells)) # Debugging check
    if(dbg) print(total.stored) # Debugging check
    
    # The following section returns the modified hydrograph and (if storage = T) 
    #  the volume of stored water to be carried over to the next year.
    if(storage) return(list(mod.hyd,total.stored))
    if(!storage) return(mod.hyd)
  }
  
}

#######################################################
# The following section generates figure layouts
#  for a specified number of alternative strategies
#  being explored.
#######################################################
plotlayouts.list<-list("one"= matrix(c(1,2,
                                       3,3),nrow=2,ncol=2,byrow=T),
                       "two"= matrix(c(1,3,
                                       2,4),nrow=2,ncol=2,byrow=T),
                       "three" = matrix(c(1,4,4,
                                          1,4,4,
                                          2,5,5,
                                          2,5,5,
                                          3,5,5,
                                          3,5,5),nrow=6,ncol=3,byrow=T),
                       "four"=matrix(c(1,5,5,
                                       1,5,5,
                                       2,5,5,
                                       2,6,6,
                                       3,6,6,
                                       3,6,6,
                                       4,6,6,
                                       4,6,6),nrow=8,ncol=3,byrow=T),
                       "five"= matrix(c(1,6,6,6,
                                        2,6,6,6,
                                        3,7,7,7,
                                        4,7,7,7,
                                        5,7,7,7),nrow=5,ncol=4,byrow=T),
                       "six" = matrix(c(1,1,7,7,7,
                                        2,5,7,7,7,
                                        3,6,8,8,8,
                                        4,0,8,8,8),nrow=4,ncol=5,byrow=T),
                       "seven" = matrix(c(1,1,8,8,8,
                                          2,5,8,8,8,
                                          3,6,9,9,9,
                                          4,7,9,9,9),nrow=4,ncol=5,byrow=T),
                       "eight" = matrix(c(1,1,9,9,9,
                                          2,6,9,9,9,
                                          3,7,10,10,10,
                                          4,8,10,10,10,
                                          5,0,10,10,10),nrow=5,ncol=5,byrow=T),
                       "nine"= matrix(c(1,1,10,10,10,
                                        2,6,10,10,10,
                                        3,7,11,11,11,
                                        4,8,11,11,11,
                                        5,9,11,11,11),nrow=5,ncol=5,byrow=T),
                       "ten"=matrix(c(1,1,11,11,11,
                                      2,7,11,11,11,
                                      3,8,11,11,11,
                                      4,9,12,12,12,
                                      5,10,12,12,12,
                                      6,0,12,12,12),nrow=6,ncol=5,byrow=T),
                       "eleven"=matrix(c(1,1,12,12,12,
                                         2,7,12,12,12,
                                         3,8,12,12,12,
                                         4,9,13,13,13,
                                         5,10,13,13,13,
                                         6,11,13,13,13),nrow=6,ncol=5,byrow=T),
                       "twelve"=matrix(c(1,1,12,12,12,
                                         2,7,12,12,12,
                                         3,8,12,12,12,
                                         4,9,13,13,13,
                                         5,10,13,13,13,
                                         6,11,13,13,13),nrow=6,ncol=5,byrow=T),
                       "thirteen"=matrix(c(1, 1,12,12,12,
                                           2, 7,12,12,12,
                                           3, 8,12,12,12,
                                           4, 9,13,13,13,
                                           5,10,13,13,13,
                                           6,11,13,13,13),nrow=6,ncol=5,byrow=T),
                       "fourteen"=matrix(c(1,1, 15,15,15,
                                           2,9, 15,15,15,
                                           3,10,15,15,15,
                                           4,11,16,16,16,
                                           5,12,16,16,16,
                                           6,13,16,16,16,
                                           7,14,16,16,16,
                                           8,0,16,16,16),nrow=8,ncol=5,byrow=T))


