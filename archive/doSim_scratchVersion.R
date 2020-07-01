###################################################################

## The simulation function
## inits is the list of starting values for the compartments,
        ## in the order Em,Fu,Un.
## simTime is the number of days the simulation should run for.
## pop is the number of people in the population.
## fidelity is the number of slices a day is divided into for numerical estimation of teh integrals
## theStart is the date the simulation begins. [It needs to match the date of the first day of observed data.]

doSim <- function(inits,
                  simTime,
                  pop,
                  fidelity,
                  theStart,
                  propEmFu,
                  propFuUn){
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
                        ## This is where the diff equations go.
                        Emz <- Emz - (Emz*propEmFu)/fidelity
                        Fuz <- Fuz + (Emz*propEmFu)/fidelity - (Fuz*propFuUn)/fidelity
                        Unz <- Unz + (Fuz*propFuUn)/fidelity
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
        
