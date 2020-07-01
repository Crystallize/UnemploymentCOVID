###################################################################

## The simulation function
## inits is the list of starting values for the compartments,
        ## in the order Em,Fu,Un.
## simTime is the number of days the simulation should run for.
## pop is the number of people in the population.
## fidelity is the number of slices a day is divided into for numerical estimation of teh integrals
## theStart is the date the simulation begins. [It needs to match the date of the first day of observed data.]
## endFur is the day into the simulation that the furlough scheme ends, and represents the switching date between the two sets of partial equations.
## peri params are those that relate to the period during furlough.
## post params are those that relate to the period after the furlough scheme ends.
## peri_furMax is the maximum number of people that can be furloughed, estimated from the timeseries.
## peri_dEmFu is the rate of transfer from employed to furloughed.
## peri_dEmUn is the rate of transfer from employed to unemployed.
## peri_UnMin is the minimum proportion of the economically active that can be unemployed, based on historical data.
## peri_dFuEm is the rate of transfer from furloughed to employed.
## peri_dUnEm is the rate of transfer from unemployed to employed.
## post_dEmUn is the rate of transfer from employed to unemployed.
## post_dFuEm is the rate of transfer from furloughed to employed.
## post_dFuUN is the rate of transfer from furloughed to unemployed.
## peri_UnMin is the minimum proportion of the economically active that can be unemployed, based on historical data.
## post_dUnEm is the rate of transfer from unemployed to employed.

doSim <- function(inits,
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
                  ){
        ## Prepare the output data structure by copying the inits vector.
        results <-
                data.frame(matrix(
                        rep(0,length(inits)*simTime),
                        ncol=length(inits),
                        nrow=simTime
                ))
        colnames(results) <-
                c("Em","Fu","Un")
        ## E is the employed fraction.
        results[1,"Em"]<-inits["Em"]
        ## F is the furloughed fraction.
        results[1,"Fu"]<-inits["Fu"]
        ## U is the unemployed fraction.
        results[1,"Un"]<-inits["Un"]
        
        ## Now step through one day at a time to our time horizon 'simTime'.
        #print(paste("simTime:",simTime))
        
        
        for (i in 1:(simTime-1)){
                ## Fetch the existing compartment values.
                Em<-results[i, "Em"]
                Fu<-results[i, "Fu"]
                Un<-results[i, "Un"]
                
                ## Initialise the interim compartments used to calculate the daily changes.
                Emz <- Em
                Fuz <- Fu
                Unz <- Un
                
                ## Cycle through the (fidelity) slices of a day used to numerically estimate the integrals.
 
                for (j in 1:fidelity){
                        ##Here detect whether we are peri or post furlough.
                        if (i<endFur){
                                ##Proceed with peri-furlough changes.
                                ## First calculate compartment transfers.
                                Em_Fu <- (Emz*(peri_furMax-Fuz)*peri_dEmFu) / fidelity
                                Em_Un <- (Emz*peri_dEmUn) / fidelity
                                Fu_Em <- (Fuz*peri_dFuEm) / fidelity
                                Fu_Un <- 0
                                Un_Em <- ((Unz-peri_UnMin)*peri_dUnEm) / fidelity
                                Um_Fu <- 0
                                ## Then combine these to calculate compartment changes.
                                Emz <- Emz - Em_Un - Em_Fu + Fu_Em + Un_Em
                                Fuz <- Fuz - Fu_Em - Fu_Un + Em_Fu + Um_Fu
                                Unz <- Unz - Un_Em - Um_Fu + Em_Un + Fu_Un                               
                        }
                        else{
                                ##proceed with post-furlough changes.
                                ## First calculate compartment transfers.
                                Em_Fu <- 0
                                Em_Un <- (Emz*post_dEmUn) / fidelity
                                Fu_Em <- (Fuz*post_dFuEm) / fidelity
                                Fu_Un <- (Fuz*post_dFuUn) / fidelity
                                Un_Em <- ((Unz-post_UnMin)*post_dUnEm) / fidelity
                                Um_Fu <- 0
                                ## Then combine these to calculate compartment changes.
                                Emz <- Emz - Em_Un - Em_Fu + Fu_Em + Un_Em
                                Fuz <- Fuz - Fu_Em - Fu_Un + Em_Fu + Um_Fu
                                Unz <- Unz - Un_Em - Um_Fu + Em_Un + Fu_Un
                                
                        }
                        
                }
                
                ## Place the new compartment values into the results object
                results[i+1,"Em"] <-Emz
                results[i+1,"Fu"] <-Fuz
                results[i+1,"Un"] <-Unz
        }
        ## Add the dates against the entries for each day.
        rownames(results) <-as.Date(0:(simTime-1),origin = theStart)
                
        ##The results are calculated as a proportion. Now calculate the number of people.
        results <- results*pop
        
        ##return the results.
        return(results)
        
        

        
}
        
