

## This contains the compartement function 'doSim':
source("200624_doSim.R")

## Remember where we are at the start.
root <- getwd()

## This array holds the observed data.
## The data here is collected from the ONS. It only covers the starting condition and desired ending condition.
theData <- read.csv("inputData.csv",stringsAsFactors = FALSE)
dateVector <- as.Date(theData$Date,format = "%d/%m/%Y", origin = "1900-01-01")
obs_Em <- theData$Employed
names(obs_Em) <- dateVector
obs_Fu <- theData$Furloughed
names(obs_Fu) <- dateVector
obs_Un <- theData$Unemployed
names(obs_Un) <-dateVector

##Initialise the parameters.
endFur <- 101
peri_furMax <- 0.29132
peri_dEmFu <- 0.1
peri_dEmUn <- 0.000008753
peri_dFuEm <- 0.006760946
peri_UnMin <- 0.039
peri_dUnEm <- 0.069
post_dEmUn <- 0.904792147
post_dFuEm <- 0.366032341
post_dFuUn <- 0.487674979
post_UnMin <- 0.039
post_dUnEm <- 0.069

## Run time parameters
theStart <- as.Date("2020-06-22")
xx <- which(theStart == as.Date(dateVector, origin="01-01-1900"))
rowNum <-length(obs_Em)
obs_Em <- obs_Em[xx:rowNum]
obs_Fu <- obs_Fu[xx:rowNum]
obs_Un <- obs_Un[xx:rowNum]

## This is the number of people in the population being simulated.
## In this case it represents the economically active population of the UK, aged 16+
pop <- 34326606

##Estimate the simluation starting parameters from the observed data.
Em<-obs_Em[1]/pop
Fu<-obs_Fu[1]/pop
Un<-obs_Un[1]/pop
## check the values add to 1.
#checkTotal <- Em+Fu+Un

##Create the init array for passing parameter values to the simulation.
inits <- c(
        Em=Em,
        Fu=Fu,
        Un=Un
)
names(inits)<- c("Em","Fu","Un")
##Format the values and display.
format(round(inits*pop,0),scientific=15)
format(sum(inits[1:3])*pop,scientific=15)

##define the simTime, i.e. how many days to simulate.
simTime <- length(obs_Em)

##define the number of simulations.

numSim <- 2

setwd(root)

## Set up tables to receive results.
        ## Separate results table for each compartment.
theResult_Em <- matrix(rep(0,simTime),nrow=simTime,ncol=numSim)
theResult_Fu <- matrix(rep(0,simTime),nrow=simTime,ncol=numSim)
theResult_Un <- matrix(rep(0,simTime),nrow=simTime,ncol=numSim)
        ##and the parameters used.
theParams <- data.frame(
        "propEmFu"=double(),
        "propFuUn"=double()
)

##Define the fidelity. This is the number of slices to divide each day into for the numerical solving of the differential equations.
fidelity <- 100000
##This step generates randomly parameterised simulations that can subsequently be sampled.
startTime <- Sys.time()

## Specify the distributions of the parameters.
## propEmFu - the proportion of employed who are furloughed each day.
propEmFu_mean <- 0.001    #mean estimated value
propEmFu_sd <- 0.0001     #standard deviation 
propEmFu_lb <- 0        #lower bound
propEmFu_ub <- 0.002      #upper bound
propEmFu_mean2 <- (propEmFu_mean-propEmFu_lb)/(propEmFu_ub-propEmFu_lb) # This is the mean value for the beta function which is the mean proportion of the interval between the lower and upper limits.
propEmFu_sd2 <- propEmFu_sd*propEmFu_mean2/propEmFu_mean
propEmFu_a <- abs(((1 - propEmFu_mean2) / propEmFu_sd2 ^ 2 - 1 / propEmFu_mean2) * propEmFu_mean2 ^ 2)  # the first (alpha) R0 shape parameter for a beta distribution. The larger the number the narrower the range.
propEmFu_b <- abs(propEmFu_a * (1 / propEmFu_mean2 - 1))    # the second (beta) parameter.
## propFuUn - the proportion of furloughed who become unemployed each day.
propFuUn_mean <- 0.001    #mean estimated value
propFuUn_sd <- 0.0001     #standard deviation 
propFuUn_lb <- 0        #lower bound
propFuUn_ub <- 0.002      #upper bound
propFuUn_mean2 <- (propFuUn_mean-propFuUn_lb)/(propFuUn_ub-propFuUn_lb) # This is the mean value for the beta function which is the mean proportion of the interval between the lower and upper limits.
propFuUn_sd2 <- propFuUn_sd*propFuUn_mean2/propFuUn_mean
propFuUn_a <- abs(((1 - propFuUn_mean2) / propFuUn_sd2 ^ 2 - 1 / propFuUn_mean2) * propFuUn_mean2 ^ 2)  # the first (alpha) R0 shape parameter for a beta distribution. The larger the number the narrower the range.
propFuUn_b <- abs(propFuUn_a * (1 / propFuUn_mean2 - 1))    # the second (beta) parameter.


## This loop generates a large number of scenarios simulations over the time interval of the observed
## data, so that we can select the best fitting examples.
for (i in 1:numSim){
       ## Print out some user feeback to help them keep track of progress.
        if (round(i/100,0)==i/100){
                cat(i,"   ")
        }
        ##Here we use the parameter distributions to pick a set of parameter values.
        # This beta distribution will generate random values.
        # See Imai et al Report 3. R0 2.6 (1.5-3.5); Zhao 2.24 (1.96-2.55) to 2.58 (2.89-4.39) (different model estimates).
        #propEmFu<- propEmFu_lb + qbeta(runif(1),propEmFu_a,propEmFu_b) * (propEmFu_ub-propEmFu_lb)
        #propFuUn<- propFuUn_lb + qbeta(runif(1),propFuUn_a,propFuUn_b) * (propFuUn_ub-propFuUn_lb)
        #print(paste("propEmFu:", propEmFu))
        #print(paste("propFuUn:", propFuUn))
        ##Now run the simulation.
        theResult <-
                doSim(inits,
                      simTime,
                      pop,
                      fidelity,
                      theStart,
                      endFur,
                      peri_furMax,
                      peri_dEmFu,
                      peri_dEmUn,
                      peri_dFuEm,
                      peri_UnMin,
                      peri_dUnEm,
                      post_dEmUn,
                      post_dFuEm,
                      post_dFuUn,
                      post_UnMin,
                      post_dUnEm
                )
        ## Take the results table and add them to to each of the compartment tables.
        theResult_Em[,i] <-as.matrix(theResult[,"Em"])
        theResult_Fu[,i] <-as.matrix(theResult[,"Fu"])
        theResult_Un[,i] <-as.matrix(theResult[,"Un"])
        
        ## theParams is a convenient array holding the parameter values of this simulation.
        ## The array is then added to the table holding the parameter values for each scenario.
        theParams <-
                rbind(theParams,
                        c(simTime,
                                pop,
                                fidelity,
                                theStart,
                                endFur,
                                peri_furMax,
                                peri_dEmFu,
                                peri_dEmUn,
                                peri_dFuEm,
                                peri_UnMin,
                                peri_dUnEm,
                                post_dEmUn,
                                post_dFuEm,
                                post_dFuUn,
                                post_UnMin,
                                post_dUnEm)
                      )
        colnames(theParams)<- c("simTime",
                                "pop",
                                "fidelity",
                                "theStart",
                                "endFur",
                                "peri_furMax",
                                "peri_dEmFu",
                                "peri_dEmUn",
                                "peri_dFuEm",
                                "peri_UnMin",
                                "peri_dUnEm",
                                "post_dEmUn",
                                "post_dFuEm",
                                "post_dFuUn",
                                "post_UnMin",
                                "post_dUnEm")
}

##Developer output for timings.
endTime <- Sys.time()
taken <- endTime - startTime
print(taken)

## Make sure we are back in the correct directory.
setwd(root)
write.csv(theResult_Em,paste("theResult_Em.csv",sep=""))
write.csv(theResult_Fu,paste("theResult_Fu.csv",sep=""))
write.csv(theResult_Un,paste("theResult_Un.csv",sep=""))

## Calculate the fitting score to the inputData
## 'weighting' is an array of weighting values that is applied to the mean squared difference
weighting <- rep(1, length(obs_Em))
## Calculate a score form the difference between the simulated number and the observed number.
## Square it to penalise distant ones more harshly.
## Base the score on Unemployed numbers only.
##Create a matrix of target values.
theTarget_Un<- as.matrix(theData$Unemployed)
##calculate the squared differences relative to the absolute target values.
gap_Un<-
        (theResult_Un[1:length(obs_Un),] - theTarget_Un[1:length(obs_Un)])^2 / theTarget_Un[1:length(obs_Un)]
##Sum the scores for each day of the simulation to give a single score for each whole scenario
scores <- apply(gap_Un,2,sum,na.rm=TRUE)
##Link the scores to the scenario parameters.
theParams <- cbind(theParams,scores)
##Now put the paramters values into order with the best fits first.
theParamsOrdered <- theParams[order(scores),]
write.csv(theParamsOrdered,"BestFittingParams.csv")

##Now lets do a longer simulation of the best fitting scenario.
## Find the parameter values of the best fitting scenario and remember them.
oendFur<- theParamsOrdered[1,"endFur"]
operi_furMax<- theParamsOrdered[1,"peri_furMax"]
operi_dEmFu<- theParamsOrdered[1,"peri_dEmFu"]
operi_dEmUn<- theParamsOrdered[1,"peri_dEmUn"]
operi_dFuEm<- theParamsOrdered[1,"peri_dFuEm"]
operi_UnMin<- theParamsOrdered[1,"peri_UnMin"]
operi_dUnEm<- theParamsOrdered[1,"peri_dUnEm"]
opost_dEmUn<- theParamsOrdered[1,"post_dEmUn"]
opost_dFuEm<- theParamsOrdered[1,"post_dFuEm"]
opost_dFuUn<- theParamsOrdered[1,"post_dFuUn"]
opost_UnMin<- theParamsOrdered[1,"post_UnMin"]
opost_dUnEm<- theParamsOrdered[1,"post_dUnEm"]


# Choose how long to project the simulations.
simTime <- 500
bestFitResult <- doSim(inits,
                       simTime,
                       pop,
                       fidelity,
                       theStart,
                       oendFur,
                       operi_furMax,
                       operi_dEmFu,
                       operi_dEmUn,
                       operi_dFuEm,
                       operi_UnMin,
                       operi_dUnEm,
                       opost_dEmUn,
                       opost_dFuEm,
                       opost_dFuUn,
                       opost_UnMin,
                       opost_dUnEm
)
write.csv(bestFitResult, "bestFitResult.csv")

## Plot the best fitting scenario and save it to a png file.
##identify which is the target day.
vertLine <- length(obs_Un)
##identify which is the target number
horizontalLine <- as.numeric(obs_Un[length(obs_Un)])
plotName<-"bestFit"
png(paste(plotName,".png",sep=""))
matplot(
        x = 0:(simTime-1),
        y = bestFitResult[,1:3],
        type = "l",
        xlab = "Time (days)",
        ylab = "Economically active individuals (Age 16+)",
        lwd = 2,
        lty = 1,
        bty = "l",
        col = c("dark green","orange","red"),
        ylim = c(0,1*pop),
        font.lab=2
)
abline(h=horizontalLine,col="grey")
abline(v=vertLine,col="grey")
## Add a legend
legend(
        "topright",c("Employed",
                "Furloughed",
                "Unemployed"
        ),
        pch = 19,
        col = c("dark green","orange","red"),
        ncol=1,
        bty = "n"
)
dev.off()
