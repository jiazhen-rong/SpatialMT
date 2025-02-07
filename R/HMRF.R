# this script contains code for HMRF-ICM process

#' Calculate pairwise distance matrix for 
#' a given spatial coordinate of transcriptomic data.
#'
#' @param spatial_coords A matrix or data frame where each row represents a point in space with coordinates.
#' @return A distance matrix where entry (i, j) represents the distance between spatial coordinates i and j.
#' @import distances
#' @export
#' @examples
#' spatial_coords <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2, byrow = TRUE)
#' calc_pairwise_distance(spatial_coords)
calc_pairwise_distance <- function(spatial_coords){
  if (!is.matrix(spatial_coords) && !is.data.frame(spatial_coords)) {
    stop("Error: spatial_coords must be a matrix or data frame.")
  }
  D= distances(spatial_coords)
  return(D)
}

#' Expectation-Maximization (EM) Algorithm for HMRF-ICM
#'
#' This function performs the EM algorithm to estimate model parameters and 
#' assign latent variables (Z_j) using an Iterated Conditional Modes (ICM) approach.
#'
#' @param X_j Numeric vector of observed variant allele counts.
#' @param N_j Numeric vector of total read counts.
#' @param pi_j0 Initial probability of background signals.
#' @param gamma Initial gamma parameter, representing contamination rate.
#' @param pi_1j Initial probability of carrier.
#' @param iterations Maximum number of EM iterations.
#' @param stop_diff Convergence threshold; the algorithm stops if the parameter change is below this.
#' @param neighbor_s_list List of nearest neighbor indices for each spot.
#' @param neighbor_s_vaf Numeric vector of neighborhood variant allele frequencies.
#' @param Wc Numeric vector representing weighted cell-type proportions from tools like RCTD.
#' @param prior_energy Logical; whether to include prior energy function in calculations.
#' @param save_path Optional directory path to save model results.
#' @param verbose Logical; whether to print progress messages.
#'
#' @return A list containing estimated parameter trajectories and final Z_j assignments.
#' \describe{
#'   \item{pi_j0_list}{Vector of estimated pi_j0 values across iterations.}
#'   \item{gamma_list}{Vector of estimated gamma values across iterations.}
#'   \item{pi_1j_list}{Vector of estimated pi_1j values across iterations.}
#'   \item{LPP_list}{Vector of log posterior probabilities across iterations.}
#'   \item{Z_j}{Final estimated latent variable assignments.}
#' }
#' @export
#'
EM <- function(X_j,N_j,pi_j0=0.1,gamma=0.5,pi_1j=0.5,iterations=10,stop_diff=1e-05,
               neighbor_s_list,neighbor_s_vaf, Wc,neighbor_num=4,icm_iter=10,prior_energy=F,
               save_path=NULL,verbose=T){
  
  # initialize parameter
  init_pi_j0=pi_j0;init_gamma=gamma;init_pi_1j=pi_1j;
  Z_j=sample(c(0,1),size=length(N_j),replace=T);
  #Z_j = ifelse(X_j / (N_j + 1e-6) > 0.2, 1, 0)  # initialize based on data
  
  # Track parameter updates
  pi_j0_list = c()
  gamma_list = c()
  pi_1j_list = c()
  LPP_list = c()
  prev_param = list(0,0,0)
  params = list(pi_j0,gamma,pi_1j)
  
  start.time0 <- Sys.time()
  
  for(i in 1:iterations){
    print(paste0("EM iter:",i))
    diff = sum(abs(unlist(prev_param) - unlist(params))) #sum(unlist(prev_param) - unlist(params))
    if(diff < stop_diff){
      print("difference small enough")
      break
    }

    # Estep - Update latent variables (Z_j) using ICM
    prev_param = params
    #delta_ij_hat = Estep(X_sub,N_sub,alpha,beta,delta_gc)
    # get Z identiy from ICM updates
    Z_j = ICM_update(X_j,N_j,Z_j,pi_j0,gamma,pi_1j,neighbor_s_vaf,Wc,neighbor_s_list, 
                     neighbor_num=neighbor_num,icm_iter=icm_iter, 
                     prior_energy=prior_energy,verbose=verbose)
    # Maximization Step in the EM, estimate each parameter based on expectation of Zij
    params = Mstep(X_j,N_j,Z_j,pi_j0,gamma,pi_1j,
                   neighbor_s_vaf,Wc,
                   prior_energy=prior_energy,opt_method="optim",neighbor_num = neighbor_num)
    # save all parameters for diagnostic
    #list(pi_j0, gamma, pi_1j) <- params
    pi_j0=params[[1]]
    pi_j0_list = append(pi_j0_list,pi_j0)
    gamma=params[[2]]
    gamma_list = append(gamma_list,gamma)
    pi_1j=params[[3]]
    pi_1j_list = append(pi_1j_list,pi_1j)
    
    # log posterior
    new_LPP = log_posterior_prob(N_j,X_j,pi_j0,gamma,Z_j,pi_1j,
                                 neighbor_s_vaf,Wc,prior_energy=prior_energy,
                                 neighbor_num=neighbor_num, neighbor_s_list=neighbor_s_list)
    LPP_list = append( LPP_list,new_LPP)
    
    if((verbose==T) & (i%%round(iterations/10)==0)){
      #end.time <- Sys.time()
      #time.taken <- end.time - start.time
      end.time0 <- Sys.time()
      time.taken0 <- difftime(end.time0, start.time0, units='mins') #end.time0 - start.time0
      message(paste0("Iteration till ", i, "th takes ", format(time.taken0,digits=2)," s"))
    }
  }
  ## get final identity
  Z_j = ICM_update(X_j,N_j,Z_j,pi_j0,gamma,pi_1j,neighbor_s_vaf,Wc,neighbor_s_list, 
                   neighbor_num=neighbor_num,
                   icm_iter=icm_iter, prior_energy=prior_energy,verbose=verbose)
  
  end.time0 <- Sys.time()
  time.taken0 <- difftime(end.time0, start.time0, units='mins') #end.time0 - start.time0
  message(paste0("Total EM time for ", i, " iterations take ",format(time.taken0,digits=2)," mins"))
  
  # save the estimated parameter for cluster c (gc)
  if(is.null(save_path)){
    dir.create("hmrf_icm_fit/")
    save(pi_j0_list,gamma_list, pi_1j_list, LPP_list, Z_j,
         file = paste0("hmrf_icm_fit/start_pij0_",init_pi_j0,"_gamma_",init_gamma,"_pij1_",init_pi_1j,".RData"))
  }else{
    dir.create(save_path,recursive = T)
    save(pi_j0_list,gamma_list, pi_1j_list, LPP_list, Z_j,
         file = paste0(save_path,"/start_pij0_",init_pi_j0,"_gamma_",init_gamma,"_pij1_",init_pi_1j,".RData"))
  }
  return(list(pi_j0_list=pi_j0_list,gamma_list=gamma_list,pi_1j_list=pi_1j_list,
              LPP_list=LPP_list,Z_j=Z_j))
}


#' Iterated Conditional Modes (ICM) Update for HMRF
#'
#' This function updates the latent variable \eqn{Z_j} using the Iterated Conditional Modes (ICM) algorithm.
#' It iteratively optimizes \eqn{Z_j} to maximize the posterior probability, considering both 
#' the likelihood of the observed data and the influence of spatial neighbors.
#'
#' @param X_j Numeric vector of observed variant allele counts for each spatial location.
#' @param N_j Numeric vector of total read counts for each spatial location.
#' @param Z_j Numeric vector of current latent variable assignments (0 = background, 1 = carrier).
#' @param pi_j0 Background probability parameter, controlling the baseline variant allele frequency.
#' @param gamma Scaling parameter for the neighborhood variant allele frequency effect.
#' @param pi_1j Carrier probability parameter, affecting locations labeled as carriers.
#' @param neighbor_s_vaf Numeric vector of neighborhood variant allele frequencies.
#' @param Wc Numeric vector of weighted cell-type proportions at each location.
#' @param neighbor_s_list List where each element contains the indices of neighboring spots.
#' @param neighbor_num Integer specifying the number of neighbors to consider for each spot.
#' @param icm_iter Maximum number of ICM iterations to perform.
#' @param prior_energy Logical; if \code{TRUE}, incorporates a prior energy function in the update.
#' @param verbose Logical; if \code{TRUE}, prints debug information during iterations.
#'
#' @return Updated numeric vector \eqn{Z_j} containing the refined latent variable assignments.
#' 
#' @details
#' The ICM algorithm updates each \eqn{Z_j} sequentially by maximizing its posterior probability:
#' \deqn{P(Z_j | X_j, \text{neighbors}) \propto P(X_j | Z_j) P(Z_j | \text{neighbors})}
#' It iterates until convergence or until the specified maximum number of iterations is reached.
#' If \code{prior_energy = TRUE}, the function incorporates a spatial prior that penalizes 
#' differences between neighboring \eqn{Z_j} values.
#'
#' @export
#'
#' @examples
#' # Example usage:
#' X_j <- c(10, 5, 8)
#' N_j <- c(100, 50, 80)
#' Z_j <- c(0, 1, 0)
#' neighbor_s_list <- list(c(2,3), c(1,3), c(1,2))
#' neighbor_s_vaf <- c(0.1, 0.05, 0.08)
#' Wc <- c(0.2, 0.1, 0.3)
#' 
#' updated_Z_j <- ICM_update(X_j, N_j, Z_j, pi_j0 = 0.1, gamma = 0.5, pi_1j = 0.5, 
#'                            neighbor_s_vaf, Wc, neighbor_s_list, neighbor_num = 4, 
#'                            icm_iter = 10, prior_energy = FALSE, verbose = TRUE)
#' print(updated_Z_j)
ICM_update <- function(X_j,N_j,Z_j,pi_j0,gamma,pi_1j,#beta,
                       neighbor_s_vaf,Wc,neighbor_s_list,
                       neighbor_num=4,icm_iter=10, prior_energy=F,verbose=F){
  #print(paste0("pi_1j",pi_1j))
  # initialization
  Z_j_old =  rep(-1, length(Z_j)) # Ensures first iteration runs
  Z_j_new = Z_j #sample(c(0,1),replace=T,length(N_j))
  iter=1
  # ICM Update
  while(any(Z_j_new != Z_j_old) & (iter <= icm_iter)){ # keep updating if spot identity changing
    if(verbose){
      print(paste0("ICM iter:",iter))
      print(paste0("equal num:",sum(Z_j_new == Z_j_old)))
    }
    Z_j_old = Z_j_new
    for(s in 1:length(N_j)){
      #print(s)
      # Z_j[s] == 0 - background
      vaf_s0 = min(0.999,pi_j0+gamma*neighbor_s_vaf[s])
      prob_0 = dbinom(x=X_j[s],size=N_j[s],prob=vaf_s0,log=T)
      # Z_-sj =1 - carrier
      vaf_s1 = min(0.999,pi_j0+gamma*neighbor_s_vaf[s]+pi_1j*Wc[s])
      prob_1 = dbinom(x=X_j[s],size=N_j[s],prob=vaf_s1,log=T)
      if(prior_energy){
        prob_0 = prob_0 + log(prior_energy_U(neighbor_num,neighbor_s_list[[s]], Z_j_new,Z=0))
        prob_1 = prob_1 + log(prior_energy_U(neighbor_num,neighbor_s_list[[s]], Z_j_new,Z=1))
      }
      # hard assignment
      #else{
      #  if(prob_0 < prob_1){
      #    Z_j_new[s] = 1
      #  }else{
      #    Z_j_new[s] = 0
      #  }
      #}
      # Probabilistic update
      prob_exp_0 <- exp(prob_0)
      prob_exp_1 <- exp(prob_1)
      #print(paste0(prob_exp_0,",",prob_exp_1))
      Z_j_new[s] <- sample(c(0, 1), size=1, prob=c(prob_exp_0, prob_exp_1))
    }
    # Early stopping if 99% of values are unchanged
    if (sum(Z_j_new == Z_j_old) / length(Z_j_old) > 0.99) {
      message("ICM early stopping: 99% of values are unchanged.")
      break
    }
    
    iter = iter+1
  }
  return(Z_j_new)
}


#' Maximization Step (M-Step) for HMRF-ICM
#'
#' This function performs the Maximization step in the Expectation-Maximization (EM) algorithm for 
#' Hidden Markov Random Fields (HMRF) using either numerical optimization (`optim()`) 
#' or alternative estimation methods.
#'
#' @param X_j Numeric vector of observed variant allele counts for each spatial location.
#' @param N_j Numeric vector of total read counts for each spatial location.
#' @param Z_j Numeric vector of current latent variable assignments (0 = background, 1 = carrier).
#' @param pi_j0 Background probability parameter, controlling the baseline variant allele frequency.
#' @param gamma Scaling parameter for the neighborhood variant allele frequency effect.
#' @param pi_1j Carrier probability parameter, affecting locations labeled as carriers.
#' @param neighbor_s_vaf Numeric vector of neighborhood variant allele frequencies.
#' @param Wc Numeric vector of weighted cell-type proportions at each location.
#' @param prior_energy Logical; if \code{TRUE}, incorporates a prior energy function in parameter estimation.
#' @param opt_method Character vector specifying the optimization method. 
#'        Can be \code{"optim"} (default, numerical optimization) or \code{"MoM"} (method of moments).
#'
#' @return A list containing the updated parameter estimates:
#' \describe{
#'   \item{pi_j0_new}{Updated background probability parameter.}
#'   \item{gamma_new}{Updated neighborhood scaling parameter.}
#'   \item{pi_1j_new}{Updated carrier probability parameter.}
#' }
#'
#' @details
#' This function updates the model parameters by maximizing the posterior probability.
#' By default, it uses `optim()` with the `"L-BFGS-B"` method to find the optimal values 
#' of \eqn{\pi_{j0}}, \eqn{\gamma}, and \eqn{\pi_{1j}} by minimizing the negative log-posterior probability.
#' If \code{prior_energy = TRUE}, it incorporates a spatial prior energy term.
#' 
#' The method can be extended to include alternative estimation approaches such as 
#' Newton-Raphson (`"NR"`) or the Method of Moments (`"MoM"`) in future implementations.
#'
#' @export
#'
#' @examples
#' # Example usage:
#' X_j <- c(10, 5, 8)
#' N_j <- c(100, 50, 80)
#' Z_j <- c(0, 1, 0)
#' neighbor_s_list <- list(c(2,3), c(1,3), c(1,2))
#' neighbor_s_vaf <- c(0.1, 0.05, 0.08)
#' Wc <- c(0.2, 0.1, 0.3)
#' 
#' # Run M-step with numerical optimization
#' updated_params <- Mstep(X_j, N_j, Z_j, pi_j0 = 0.1, gamma = 0.5, pi_1j = 0.5, 
#'                         neighbor_s_vaf, Wc, prior_energy = FALSE, opt_method = "optim")
#' 
#' print(updated_params)
Mstep <- function(X_j, N_j, Z_j, pi_j0, gamma, pi_1j,
                  neighbor_s_vaf, Wc, prior_energy=FALSE, 
                  opt_method=c("optim", "MoM"),neighbor_num = 4) {
  
  if (opt_method=="optim") {  
    data_list <- list(X_j=X_j, N_j=N_j, neighbor_s_vaf=neighbor_s_vaf, 
                      Wc=Wc, neighbor_s_list=neighbor_s_list, Z_j=Z_j)
    
    if (prior_energy) {
      param_fit <- optim(par=c(pi_j0, gamma, pi_1j), 
                         fn=neg_posterior_prob_prior,
                         data=data_list,
                         neighbor_num = neighbor_num,
                         method="L-BFGS-B", 
                         lower=c(1e-5, 1e-5, 0.1), upper=c(0.99, 0.5, 0.99))
    } else {
      # param_fit <- optim(par=c(pi_j0, gamma, pi_1j), 
      #                    fn=neg_posterior_prob,
      #                    data=data_list,
      #                    #method="L-BFGS-B", 
      #                    method="BFGS")#,
      #                    #lower=c(1e-5, 1e-5, 0.1), upper=c(0.99, 0.5, 0.99))
      #                    #lower = c(1e-6, 1e-3, 0.05),  # More relaxed lower bounds
      #                    #upper = c(0.999, 0.8, 0.99))  # Higher gamma and pi_1j bounds
      param_fit <- optim(par=logit(c(pi_j0, gamma, pi_1j)),  # Convert params to logit space
                         fn=neg_posterior_prob_logit,
                         data=data_list, 
                         method="BFGS")
                         
    }
    print(param_fit$convergence)  # 0 means successful optimization
    print(param_fit$par)  # Print updated parameters
    updated_params <- sigmoid(param_fit$par)  # Convert back to [0,1]
    
    #pi_j0_new <- param_fit$par[1]
    #gamma_new <- param_fit$par[2]
    #pi_1j_new <- param_fit$par[3]

    pi_j0_new <- updated_params[1]
    gamma_new <- updated_params[2]
    pi_1j_new <- updated_params[3]
    
  } else if ("NR" %in% opt_method) {  
    # Implement Newton-Raphson if needed
  } 
  
  return(list(pi_j0_new, gamma_new, pi_1j_new))
}


logit <- function(x) log(x / (1 - x))
sigmoid <- function(x) 1 / (1 + exp(-x))


neg_posterior_prob_logit <- function(param_logit, data){
  # Convert logit-transformed parameters back to probability space
  param <- sigmoid(param_logit)  
  return(neg_posterior_prob(param, data))  # Compute negative log posterior
 
}

#' Compute Negative Log Posterior Probability (Without Prior Energy)
#'
#' This function calculates the negative log posterior probability for a given set of parameters
#' without incorporating a spatial prior energy function. It is used in the M-step of the EM algorithm 
#' for optimizing \eqn{\pi_{j0}}, \eqn{\gamma}, and \eqn{\pi_{1j}}.
#'
#' @param param Numeric vector containing the model parameters: 
#'   \describe{
#'     \item{\code{param[1]}}{Initial probability of background \eqn{\pi_{j0}}.}
#'     \item{\code{param[2]}}{Neighborhood effect scaling parameter \eqn{\gamma}.}
#'     \item{\code{param[3]}}{Probability of carrier identity \eqn{\pi_{1j}}.}
#'   }
#' @param data A list containing the following elements:
#'   \describe{
#'     \item{\code{X_j}}{Numeric vector of observed variant allele counts.}
#'     \item{\code{N_j}}{Numeric vector of total read counts.}
#'     \item{\code{Z_j}}{Numeric vector of latent variable assignments (0 = background, 1 = carrier).}
#'     \item{\code{neighbor_s_vaf}}{Numeric vector of neighborhood variant allele frequencies.}
#'     \item{\code{Wc}}{Numeric vector of weighted cell-type proportions at each location.}
#'   }
#'
#' @return The negative log posterior probability (a scalar value).
#'
#' @details
#' This function computes:
#' \deqn{NPP = - \sum \log P(X_j | Z_j, \pi_{j0}, \gamma, \pi_{1j})}
#' where \eqn{P(X_j | Z_j)} is modeled using the binomial likelihood.
#' The function is used for parameter estimation during the M-step in the EM algorithm.
#'
#' @export
#'
#' @examples
#' # Define parameter values
#' params <- c(0.1, 0.5, 0.5)
#'
#' # Define input data
#' data_list <- list(
#'   X_j = c(10, 5, 8),
#'   N_j = c(100, 50, 80),
#'   Z_j = c(0, 1, 0),
#'   neighbor_s_vaf = c(0.1, 0.05, 0.08),
#'   Wc = c(0.2, 0.1, 0.3)
#' )
#'
#' # Compute negative log posterior probability
#' neg_log_posterior <- neg_posterior_prob(params, data_list)
#' print(neg_log_posterior)
neg_posterior_prob  <-function(param,data){
  pi_j0<-param[1]; gamma <- param[2]; pi_1j <- param[3]
  
  # negative posterior probability
  prob_all = c() # log probability
  for(s in 1:length(data$N_j)){
    if(data$Z_j[s] == 0){
      vaf_s0 = min(0.999,pi_j0+gamma*data$neighbor_s_vaf[s])
      prob = dbinom(x=data$X_j[s],size=data$N_j[s],prob=vaf_s0,log=T)
    }else{ # Z_-sj =1 - carrier
      vaf_s1 = min(0.999,pi_j0+gamma*data$neighbor_s_vaf[s]+pi_1j*data$Wc[s])
      prob = dbinom(x=data$X_j[s],size=data$N_j[s],prob=vaf_s1,log=T)
    }
    prob_all <- append(prob_all,prob)
    # capping largest and smallest probability
  }
  
  # Debugging: Check for NA values
  if (any(is.na(prob_all))) warning("Total of ",sum(is.na(prob_all)), " NA values detected in prob_all")
  NPP = -sum(prob_all,na.rm=T)
  
  return(NPP)
}

#' Compute Negative Log Posterior Probability (With Prior Energy)
#'
#' This function calculates the negative log posterior probability while incorporating
#' a spatial prior energy function to encourage spatial smoothness in the HMRF model.
#'
#' @param param Numeric vector containing the model parameters: 
#'   \describe{
#'     \item{\code{param[1]}}{Initial probability of background \eqn{\pi_{j0}}.}
#'     \item{\code{param[2]}}{Neighborhood effect scaling parameter \eqn{\gamma}.}
#'     \item{\code{param[3]}}{Probability of carrier identity \eqn{\pi_{1j}}.}
#'   }
#' @param data A list containing the following elements:
#'   \describe{
#'     \item{\code{X_j}}{Numeric vector of observed variant allele counts.}
#'     \item{\code{N_j}}{Numeric vector of total read counts.}
#'     \item{\code{Z_j}}{Numeric vector of latent variable assignments (0 = background, 1 = carrier).}
#'     \item{\code{neighbor_s_vaf}}{Numeric vector of neighborhood variant allele frequencies.}
#'     \item{\code{Wc}}{Numeric vector of weighted cell-type proportions at each location.}
#'     \item{\code{neighbor_s_list}}{List where each element contains the indices of neighboring spots.}
#'   }
#' @param neighbor_num Integer specifying the number of neighbors to consider for each spot.
#'
#' @return The negative log posterior probability (a scalar value).
#'
#' @details
#' This function computes:
#' \deqn{NPP = - \sum \left[ \log P(X_j | Z_j, \pi_{j0}, \gamma, \pi_{1j}) + \log P(Z_j | \text{neighbors}) \right]}
#' where:
#' - \eqn{P(X_j | Z_j)} is modeled using the binomial likelihood.
#' - \eqn{P(Z_j | \text{neighbors})} is a spatial prior energy function that encourages smoothness.
#'
#' The function is used for parameter estimation during the M-step in the EM algorithm with prior energy.
#'
#' @export
#'
#' @examples
#' # Define parameter values
#' params <- c(0.1, 0.5, 0.5)
#'
#' # Define input data
#' data_list <- list(
#'   X_j = c(10, 5, 8),
#'   N_j = c(100, 50, 80),
#'   Z_j = c(0, 1, 0),
#'   neighbor_s_vaf = c(0.1, 0.05, 0.08),
#'   Wc = c(0.2, 0.1, 0.3),
#'   neighbor_s_list = list(c(2,3), c(1,3), c(1,2))
#' )
#'
#' # Compute negative log posterior probability with prior energy
#' neg_log_posterior <- neg_posterior_prob_prior(params, data_list, neighbor_num = 4)
#' print(neg_log_posterior)
neg_posterior_prob_prior <- function(param, data, neighbor_num=4) {
  pi_j0 <- param[1]; gamma <- param[2]; pi_1j <- param[3]
  neighbor_num <- length(data$neighbor_s_list[[1]])  # Corrected from `len()`
  
  # Preallocate vector for efficiency
  prob_all <- numeric(length(data$N_j))
  
  for (s in seq_along(data$N_j)) {
    if (data$Z_j[s] == 0) {
      vaf_s0 <- min(0.999, pi_j0 + gamma * data$neighbor_s_vaf[s])
      prob <- dbinom(x = data$X_j[s], size = data$N_j[s], prob = vaf_s0, log = TRUE) +
        log(prior_energy_U(neighbor_num, data$neighbor_s_list[[s]], data$Z_j, Z = 0))
    } else {  # Z_j[s] == 1 (carrier)
      vaf_s1 <- min(0.999, pi_j0 + gamma * data$neighbor_s_vaf[s] + pi_1j * data$Wc[s])
      prob <- dbinom(x = data$X_j[s], size = data$N_j[s], prob = vaf_s1, log = TRUE) +
        log(prior_energy_U(neighbor_num, data$neighbor_s_list[[s]], data$Z_j, Z = 1))
    }
    prob_all[s] <- prob
  }
  
  # Debugging: Check for NA values
  if (any(is.na(prob_all))) warning("NA values detected in prob_all")
  
  # Negative log posterior probability
  NPP <- -sum(prob_all, na.rm = TRUE)
  return(NPP)
}

#' Compute Log Posterior Probability of the Model
#'
#' This function calculates the log posterior probability for the HMRF model, 
#' incorporating binomial likelihood and optional spatial prior energy terms.
#'
#' @param N_j Numeric vector of total read counts.
#' @param X_j Numeric vector of observed variant allele counts.
#' @param pi_j0 Initial probability of background.
#' @param gamma Neighborhood effect scaling parameter.
#' @param Z_j Numeric vector of latent variable assignments (0 = background, 1 = carrier).
#' @param pi_1j Carrier probability parameter.
#' @param neighbor_s_vaf Numeric vector of neighborhood variant allele frequencies.
#' @param Wc Numeric vector of weighted cell-type proportions at each location.
#' @param prior_energy Logical; if \code{TRUE}, incorporates a spatial prior energy function.
#' @param neighbor_num Integer specifying the number of neighbors (needed if \code{prior_energy = TRUE}).
#' @param neighbor_s_list List where each element contains the indices of neighboring spots.
#'
#' @return The log posterior probability (a scalar value).
#'
#' @details
#' This function computes:
#' \deqn{LPP = \sum \log P(X_j | Z_j, \pi_{j0}, \gamma, \pi_{1j}) + \log P(Z_j | \text{neighbors})}
#' where:
#' - \eqn{P(X_j | Z_j)} is modeled using the binomial likelihood.
#' - \eqn{P(Z_j | \text{neighbors})} is a spatial prior energy function (if enabled).
#'
#' @export
#'
#' @examples
#' # Example usage
#' log_posterior_prob(N_j, X_j, pi_j0=0.1, gamma=0.5, Z_j, pi_1j=0.5,
#'                    neighbor_s_vaf, Wc, prior_energy=TRUE,
#'                    neighbor_num=4, neighbor_s_list=neighbor_s_list)
#'                    
log_posterior_prob <- function(N_j,X_j,pi_j0=0.1,gamma=0.5,Z_j,pi_1j=0.5,
                               neighbor_s_vaf,Wc,prior_energy=F,
                               neighbor_num=4, neighbor_s_list=NULL){
  # log posterior probability
  prob_all = c() # log probability
  for(s in 1:length(N_j)){
    if(Z_j[s] == 0){
      vaf_s0 = min(0.999,pi_j0+gamma*neighbor_s_vaf[s])
      prob = dbinom(x=X_j[s],size=N_j[s],prob=vaf_s0,log=T)
    }else{ # Z_-sj =1 - carrier
      vaf_s1 = min(0.999,pi_j0+gamma*neighbor_s_vaf[s]+pi_1j*Wc[s])
      prob = dbinom(x=X_j[s],size=N_j[s],prob=vaf_s1,log=T)
    }
    if(prior_energy){ # double check
      prob = prob + log(prior_energy_U(neighbor_num,neighbor_s_list[[s]], Z=0)) +
        log(prior_energy_U(neighbor_num,neighbor_s_list[[s]], Z=1))
    }
    prob_all <- append(prob_all,prob)
    # capping largest and smallest probability
  }
  # Debugging: Check for NA values
  if (any(is.na(prob_all))) warning("NA values detected in prob_all")
  
  LPP = sum(prob_all,na.rm=T)
  return(LPP)
}


#' Compute Prior Energy Function for HMRF Model
#'
#' This function calculates the spatial prior probability \eqn{P(Z_j | \text{neighbors})}
#' given the latent variable \eqn{Z_j} and its neighboring spots.
#'
#' @param neighbor_spots Numeric vector of indices representing neighboring spots.
#' @param Z_j Numeric vector of latent variable assignments for all spots (0 = background, 1 = carrier).
#' @param Z Scalar (0 or 1), the latent state being considered for the given spot s.
#'
#' @return The prior energy probability value.
#'
#' @details
#' The function computes:
#' \deqn{U_x = \frac{1}{(2\pi)^{n/2}} \exp(- V_c)}
#' where:
#' - \eqn{V_c} is the sum of pairwise neighbor mismatches.
#' - \eqn{n} is the number of neighbors.
#'
#' If no neighbors are present, the function returns 1 (neutral effect).
#'
#' @export
#'
#' @examples
#' Z_j <- c(0, 1, 1, 0, 1)
#' neighbor_spots <- c(2, 3, 5)
#' prior_energy_U(neighbor_spots, Z_j, Z = 1)
prior_energy_U <- function(neighbor_num=4,neighbor_spots, Z_j,Z){
  total_n = length(neighbor_spots) 
  #neighbor_spots = neighbor_s_list[[s]]
  # Compute mismatch penalty
  Vc_sum = sum(1/2 *(1-(Z == Z_j[neighbor_spots])))
  U_x = 1/((2*pi)^(total_n/2)) * exp(-Vc_sum)
  return(U_x)
}

# This function calculates probability of being 0 or being 1
spotwise_prob  <-function(X_j,N_j,Z_j,pi_j0,gamma,pi_1j,
                          neighbor_s_vaf,neighbor_s_list,Wc,
                          prior_energy=F){
  
  # negative posterior probability
  prob_0_all = c() # log probability
  prob_1_all = c() # log probability
  for(s in 1:length(N_j)){
    vaf_s0 = min(0.999,pi_j0+gamma*neighbor_s_vaf[s])
    prob_0 = dbinom(x=X_j[s],size=N_j[s],prob=vaf_s0,log=T)
    
    vaf_s1 = min(0.999,pi_j0+gamma*neighbor_s_vaf[s]+pi_1j*Wc[s])
    prob_1 = dbinom(x=X_j[s],size=N_j[s],prob=vaf_s1,log=T)
    if(prior_energy==T){
      prob_0 = prob + log(prior_energy_U(s,neighbor_num,neighbor_s_list, Z=0))
      prob_1 = prob + log(prior_energy_U(s,neighbor_num,neighbor_s_list, Z=1))
    }
    
    prob_0_all <- append(prob_0_all,prob_0)
    prob_1_all <- append(prob_1_all,prob_1)
    # capping largest and smallest probability?
  }
  return(list( prob_0_all,prob_1_all))
}


# This function plots estimated Z_j
# @ import patchwork
spatial_plot <- function(spots,Z_j,save_path=NULL,plot_ratio=T,spot_prob=NULL,val_cap=3,
                         vaf_cov=T,N_j=NULL,vaf_j=NULL,thresh=0,joint_title=""){
  # create plotting data frame
  if(plot_ratio){
    plot_df = as.data.frame(cbind(spots,as.factor(Z_j),exp(spot_prob[[1]]),exp(spot_prob[[2]])))
    colnames(plot_df) = c("x","y","Carrier","log_Prob_0","log_Prob_1")
    plot_df$log_Ratio = spot_prob[[2]] - spot_prob[[1]]
    # cap the values
    plot_df$log_Ratio = pmax(plot_df$log_Ratio,-val_cap)
    plot_df$log_Ratio = pmin(plot_df$log_Ratio,val_cap)
    
  }else{
    plot_df = as.data.frame(cbind(spots,as.factor(Z_j)))
    colnames(plot_df) = c("x","y","Carrier")
  }
  plot_df = data.frame(plot_df)
  if(vaf_cov){
    print("Plotting VAF and coverage")
    plot_df$vaf=vaf_j
    plot_df$cov=N_j
    plot_df = data.frame(plot_df)
    # plot VAF and cov
    g5 <-  ggplot(plot_df,aes(x=-y,y=-x)) + geom_point(size=0.5,color="grey") +
      geom_point(data=plot_df[plot_df$vaf > 0,], aes(x=-y,y=-x,color=vaf),size=1) +
      scale_color_gradient(low="grey",high="red") +
      theme_bw() +  ggtitle(paste0("VAF >0: ",sum(plot_df$vaf > 0)," spots"))+
      theme(plot.title=(element_text(size=20)))
    
    g6 <-  ggplot(plot_df,aes(x=-y,y=-x)) + geom_point(size=0.5,color="grey") +
      geom_point(data=plot_df[plot_df$cov > 0,], aes(x=-y,y=-x,color=cov),size=1) +
      scale_color_gradient(low="grey",high="red") +
      theme_bw() +  ggtitle(paste0("Coverage(N_j) >0: ",sum(plot_df$N_j > 0)," spots"))+
      theme(plot.title=(element_text(size=20)))
  }
  
  
  # plot Zj
  g1 <- ggplot(plot_df,aes(x=-y,y=-x)) + geom_point(color="grey",size=1) +
    geom_point(data=plot_df[plot_df$Carrier == 1,],
               aes(x=-y,y=-x),color="red",size=1) +
    theme_bw() +  ggtitle(paste0("Final Z_j (threshold=0):",sum(Z_j)," spots"))+
    theme(plot.title=(element_text(size=20)))
  if(plot_ratio){
    print("Plotting Ratios ...")
    g2 <- ggplot(plot_df,aes(x=-y,y=-x,color=log_Prob_0)) + geom_point(size=1) +
      scale_color_gradient(low="grey",high="red") +
      theme_bw() +  ggtitle("Prob Z=0")+
      theme(plot.title=(element_text(size=20)))
    g3 <- ggplot(plot_df,aes(x=-y,y=-x,color=log_Prob_1)) + geom_point(size=1) +
      scale_color_gradient(low="grey",high="red") +
      theme_bw() +  ggtitle("Prob Z=1")+
      theme(plot.title=(element_text(size=20)))
    
    myPalette <- colorRampPalette(c("blue","grey","red"))
    sc <- scale_colour_gradientn(colours = myPalette(100), limits=c(-3,3))
    
    g4 <- ggplot(plot_df,aes(x=-y,y=-x,color=log_Ratio)) + geom_point(size=1,) +
      sc+
      #scale_color_gradient2(low="blue",mid="grey",high="red") +
      geom_point(data=plot_df[plot_df$log_Ratio > 0,], aes(x=-y,y=-x,color=log_Ratio),size=1) +
      theme_bw() +  ggtitle("Log Ratio of Z=1 over Z=0")+
      theme(plot.title=(element_text(size=20)))
    
    if(vaf_cov){
      g=plot_grid(g5,g6,g1,g2,g3,g4)
    }else{
      g=plot_grid(g1,g2,g3,g4)
    }
    # add joint title
    title_theme <- ggdraw() +
      draw_label(joint_title, x = 0.05, hjust = 0,size=25)
    g_final<- plot_grid(title_theme, g, ncol = 1, rel_heights = c(0.05, 1))
    
    if(!is.null(save_path)){
      print("DEBUG-yes")
      png(paste0(save_path,"/final_carrier_identity.png"),
          height=500*floor(sqrt(length(g$layers))),width=500*ceiling(sqrt(length(g$layers))))
      print(g_final)
      dev.off()
    }
    return(g_final)
    
  }else{
    if(!is.null(save_path)){
      png(paste0(save_path,"/final_carrier_identity.png"),height=800,width=800)
      print(g1)
      dev.off()
    }
    return(g1)
  }
}


