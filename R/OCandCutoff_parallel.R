#' OCC.Table
#' 
#' This is the main function used to generate simulation data, and evaluate the operating characteristics and cutoff value based on the simulation data.
#' Note that mclapply is used for parallel programming, which works on Unix-style operating systems.
#' 
#' @param N.sim   Number of simulations.
#' @param N.max   Maximum number of patients to enroll.
#' @param mu.trt  True mean for treatment arm (log).
#' @param Sigma.trt  True variance for treatment arm (log).
#' @param mu.ctrl True mean for control arm (log).
#' @param Sigma.ctrl True variance for control arm (log).
#' @param cens_upper Upper limit for the censoring time.
#' @param design A numeric value indicating the type of design. 1 = proposed design, 2 = time to recurrence design, 3= time to death design, 4 = time to first event design.
#' @param cohort  Interim cohort. 
#' @param recruit.int = Recruitment interval.
#' @param m0  Prior mean for mu.trt/mu.ctrl.
#' @param L0  Prior variance for mu.trt/mu.ctrl. 
#' @param v0  For the proposed design, this is the prior degrees of freedom for Sigma.trt/Sigma.ctrl. For the traditional designs, v0/2 is the prior shape for Sigma.trt/Sigma.ctrl.
#' @param S0  For the proposed design, this is the prior scale matrix for Sigma.trt/Sigma.ctrl. For the traditional designs, v0*S0/2 is the prior scale for Sigma.trt/Sigma.ctrl.
#' @param time_max The upper limit for the recurrence and death time sampled from truncated normal. This will set the upper limit to to time_max rather than Inf.
#' @param eta A pre-specified lower bound of acceptable performance based on historical information.
#' @param lambda  Cutoff parameter.
#' @param thin_MCMC  Thinning degree.
#' @param Niter Number of iterations for gibbs sampler.
#' @returns A list with the following components:\tabular{ll}{
#'    \code{mysim1} \tab A numeric vector. mysim1 = c(PRN, PEN,  EN), where PRN is the average percentage of trials that are not stopped, PEN is the average percentage of early termination, and EN is the average number of patients.  \cr
#'    \tab \cr
#'    \code{mu} \tab A numeric vector. mu = c(mu.trt,mu.ctrl).  \cr
#'    \tab \cr
#'    \code{Sigma.trt} \tab True covariance matrix for treatment arm (log). \cr
#'    \tab \cr
#'   \code{Sigma.ctrl} \tab True covariance matrix for control arm (log). \cr
#'    \tab \cr
#'    \code{lambda} \tab Cutoff parameter. \cr
#' }
#' @examples
#' \dontrun{
#' OCC.Table(N.sim = 2, N.max = 100,mu.trt = c(log(2), log(5)), 
#' Sigma.trt = matrix(c(1, 0.5, 0.5, 1), ncol=2),
#' mu.ctrl= c(log(2), log(5)),
#' Sigma.ctrl = matrix(c(1, 0.5, 0.5, 1), ncol=2), 
#' cens_upper = 20,
#' design = 1, 
#' cohort = c(40,60,80,100),
#' recruit.int  = 0.25,
#' m0 = c(0,0),L0 = diag(10^6, 2), 
#' v0 = 4,S0 = diag(10^(-6), 2),
#' time_max = 20,eta = 1,lambda = 0.25,
#'  thin_MCMC = 5,Niter = 100000)
#' }




#' @export
OCC.Table<- function(N.sim,N.max,mu.trt,Sigma.trt,
                     mu.ctrl,Sigma.ctrl,cens_upper,design, cohort, recruit.int,
                     m0,L0, v0, S0,time_max,eta, lambda, thin_MCMC,
                     Niter){
  
myData<-array(dim=c(N.max,7,N.sim));   

for (j in 1:N.sim){
  
  # Simulate randomization groups
  rand.list = sample(c(0,1), N.max, replace = T)
  
  # Simulatetime to recurrence and death
  trt.events <- 	as.data.frame(MASS::mvrnorm(sum(rand.list),
                                       mu= mu.trt,
                                       Sigma= Sigma.trt))     
  trt.events$id =  which(rand.list==1); trt.events$grp = 1; 
  
  ctrl.events <- as.data.frame(MASS::mvrnorm(sum(rand.list==0),
                                       mu= mu.ctrl,
                                       Sigma= Sigma.ctrl))     
  
  ctrl.events$id = which(rand.list==0); ctrl.events$grp = 0
  
  # simulate time to censor
  time.censor = stats::runif(N.max, 0, cens_upper); 
  
  eventsall = rbind(trt.events, ctrl.events); eventsall = eventsall[order(eventsall$id),]
  
  myData[,1,j]<- exp(eventsall[,1])                 # time to recurrence 
  myData[,2,j]<- exp(eventsall[,2])                      # time to death
  myData[,3,j]<- time.censor
  myData[,4,j]<- eventsall[,4]                      # treatment group assnments
  myData[,5,j]<-  eventsall[,3]                     # id
}

# Create dimension names
dim_names <- list(
  NULL,
  Variables = c("recurrence_t", "death_t", "censor_t",   "group", "id", "delta1", "delta2"),
  Simulations = paste("Simulation", 1:N.sim)
)

# Assign the dimension names to the array
dimnames(myData) <- dim_names


  
  
  
  stop.all<-stop.early <-pts.all <-rep(0, N.sim);
  
 
  
  run_simulation <- function(i.sim) {      
    
    tryCatch({
  Time.entry <-Time.current <-n.current <-j.cohort<-0;
  trial.stop<- trialER.stop<-pts.stop<-  0;
  
  
  while (n.current < N.max &&  trial.stop !=1)
  {
    j.cohort<-j.cohort+1;
    n.current<- cohort[j.cohort];
    Time.entry <- recruit.int*(c(1:n.current)-1)        ; ## take into account recruit interval;
    Time.current <- recruit.int*(n.current-1)  ; 
    

    if (design == 1) {
      
      currentData =  myData[1:n.current-1,,i.sim] 
      result = winratio(currentdd = currentData,
                        n_current = n.current,
                        Time_current = Time.current,
                        Time_entry = Time.entry,
                        m0 = m0,
                        L0 = L0,
                        v0 = v0,
                        S0 = S0,
                        time_max = time_max,
                        eta = eta,
                        lambda = lambda,
                        N.max = N.max,
                        thin_MCMC = thin_MCMC,
                        Niter = Niter)
     
      if (!is.na(result$probs) && result$probs <= result$cutoff) { trial.stop<-1; trialER.stop <-1 }
    }else if (design %in% c(2,3,4)){
      if(design==2|design == 3){
        currentData = myData[1:n.current-1,c(design-1, 3:6),i.sim]
      }else{
        currentData = data.frame(evnt.time = pmin(myData[1:n.current-1,1,i.sim],myData[1:n.current-1,2,i.sim] ), 
                                 myData[1:n.current-1,c(3:6),i.sim])
      }
     
      
      # update actual censoring time at this interim look
      currentData[,2] = pmin(Time.current-Time.entry[-n.current], currentData[,2])
      currentData[,5] <- as.integer(currentData[,1] <= currentData[,2])
      
      # update event time to reflect the observed data at this interim look 
      currentData[,1] = ifelse(currentData[,5]==1, currentData[,1], NA)
      
      # number enrolled to each group
      n.current.trt = sum(currentData[,3])
      n.current.ctrl =  n.current - 1 - n.current.trt 
      
      # separate data to trt and control
      currentData.trt = currentData[currentData[,3]==1,]
      currentData.ctrl = currentData[currentData[,3]==0,]

      # update theta
      trt.post = update_theta_univariate(N_iter = Niter, dd = currentData.trt, 
                                                     n = n.current.trt, L0 = L0,m0 = m0, v0 = v0, S0 = S0,
                                                     time_max = time_max)
      ctrl.post = update_theta_univariate(N_iter = Niter, dd = currentData.ctrl, 
                                                      n = n.current.ctrl, L0 = L0, m0 = m0, v0 = v0, S0 = S0,
                                                      time_max = time_max)
      burn_MCMC = as.integer(0.3*Niter)
      idxs <- seq(burn_MCMC, Niter, by = thin_MCMC)
      
      probs = mean(trt.post$theta[ idxs,1] - ctrl.post$theta[ idxs,1] > eta)
      cutoff =  (n.current/N.max)**lambda
      if (!is.na(probs) && probs <= cutoff) { trial.stop<-1; trialER.stop <-1 }
   }
  }
  
  pts.stop<-n.current;
  if (pts.stop==N.max) {trialER.stop<-0}      
  ### early stop for the last cohort; trial.stop =1, while not early stoped.
  
  # Return the results
  list(trial.stop = trial.stop, trialER.stop = trialER.stop, pts.stop = pts.stop)
}, error = function(e) {
  # Handle the error (e.g., print a message or log it)
  cat("Error in iteration", i.sim, ":", conditionMessage(e), "\n")
})
  }


  sim_results <- parallel::mclapply(1:N.sim, run_simulation, mc.cores = parallel::detectCores())
  stop.all <- sapply(sim_results, function(res) res$trial.stop)
  stop.early <- sapply(sim_results, function(res) res$trialER.stop)
  pts.all <- sapply(sim_results, function(res) res$pts.stop)

  
  MyRaw <- data.frame(PRN=1-stop.all, PEN=stop.early, EN=pts.all);
  mysim1<-(apply(MyRaw,2,mean)); 
  mysim1<- round(mysim1,3);
  
  return(list(mysim1 = mysim1, mu = c(mu.trt,mu.ctrl), Sigma.trt = Sigma.trt, Sigma.ctrl = Sigma.ctrl,
              lambda = lambda))
  
}
