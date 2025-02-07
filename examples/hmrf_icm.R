# this script contains code for HMRF-ICM process
library(distances)
calc_pairwise_distance <- function(spatial_coords){
  D= distances(spatial_coords)
}

# this function calculates posterior probability of the model
# D is the distance matrix
# s represents the number of closest neighbors to consider
log_posterior_prob <- function(N_j,X_j,pi_j0=0.1,gamma=0.5,Z_j,pi_1j=0.5,neighbor_s_vaf,Wc,prior_energy=F){
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
      prob + log(prior_energy_U(neighbor_num,neighbor_s_list[[s]], Z=1))
    }
    prob_all <- append(prob_all,prob)
    # capping largest and smallest probability
  }

  LPP = sum(prob_all,na.rm=T)
  return(LPP)
}


# This function icm updates Zsj (varrier identity)
ICM_update <- function(X_j,N_j,Z_j,pi_j0,gamma,pi_1j,#beta,
                       neighbor_s_vaf,Wc,neighbor_s_list,
                       neighbor_num=4,icm_iter=10, prior_energy=F,verbose=F){
  #print(paste0("pi_1j",pi_1j))
  Z_j_old = Z_j
  Z_j_new = sample(c(0,1),replace=T,length(N_j))
  iter=1
  while(any(Z_j_new != Z_j_old) | (iter > icm_iter)){ # keep updating if spot identity changing
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
      }else{
        if(prob_0 < prob_1){
          Z_j_new[s] = 1
        }else{
          Z_j_new[s] = 0
        }
      }
    }
    iter = iter+1
  }
  return(Z_j_new)
}

# calculates P(theta| x_neighbor)
prior_energy_U <- function(neighbor_num=4,neighbor_spots, Z_j,Z){
  total_n = neighbor_num
  #neighbor_spots = neighbor_s_list[[s]]
  Vc_sum = sum(1/2 *(1-(Z == Z_j[neighbor_spots])))
  U_x = 1/((2*pi)^(total_n/2)) * exp(-Vc_sum)
  return(U_x)
}


# This function calculates negative log posterior probability
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
# This function calculates negative log posterior probability
# with prior field density
neg_posterior_prob_prior  <-function(param,data){
  pi_j0<-param[1]; gamma <- param[2]; pi_1j <- param[3]
  neighbor_num = len(data$neighbor_s_list[[1]])
  # negative posterior probability
  prob_all = c() # log probability
  for(s in 1:length(N_j)){
    if(Z_j[s] == 0){
      vaf_s0 = min(0.999,pi_j0+gamma*data$neighbor_s_vaf[s])
      prob = dbinom(x=data$X_j[s],size=data$N_j[s],prob=vaf_s0,log=T)+
        log(prior_energy_U(neighbor_num,data$neighbor_s_list[[s]], Z_list,Z=0))
    }else{ # Z_-sj =1 - carrier
      vaf_s1 = min(0.999,pi_j0+gamma*data$neighbor_s_vaf[s]+pi_1j*data$Wc[s])
      prob = dbinom(x=data$X_j[s],size=data$N_j[s],prob=vaf_s1,log=T)+
        log(prior_energy_U(s,neighbor_num,data$neighbor_s_list[[s]], Z_list,Z=1))
    }
    prob_all <- append(prob_all,prob)
    # capping largest and smallest probability
  }
  NPP = -sum(prob_all,na.rm=T)
  return(NPP)
}

Mstep <- function(X_j,N_j,Z_j,pi_j0,gamma,pi_1j,
                  neighbor_s_vaf,Wc,
                  prior_energy=F,opt_method=c("optim","MoM")){

  # optim to find alpha and beta -- way too slow
  if(opt_method == "optim"){ # R's optimization function
    #start.time <- Sys.time()
    if(prior_energy){
      param_fit = optim(par=c(pi_j0, gamma, pi_1j), fn=neg_posterior_prob_prior,
                        data=data.frame(X_j=c(as.matrix(X_j)),N_j=c(as.matrix(N_j)),
                                        neighbor_s_vaf=c(as.matrix(neighbor_s_vaf)),
                                        Wc=c(as.matrix(Wc)),#c(as.matrix(neighbor_s_vaf)),
                                        neighbor_s_list=neighbor_s_list,
                                        Z_j=Z_j),
                        method = "L-BFGS-B",lower=c(1e-5,1e-5,0.1), upper=c(0.99,0.5,0.99))

    }else{
      param_fit = optim(par=c(pi_j0, gamma, pi_1j), fn=neg_posterior_prob ,
                        data=data.frame(X_j=c(as.matrix(X_j)),N_j=c(as.matrix(N_j)),
                                        neighbor_s_vaf=c(as.matrix(neighbor_s_vaf)),
                                        Wc=c(as.matrix(Wc)),#c(as.matrix(neighbor_s_vaf)),
                                        #neighbor_s_list=neighbor_s_list,
                                        Z_j=Z_j),
                        method = "L-BFGS-B",lower=c(1e-5,1e-5,0.1), upper=c(0.99,0.5,0.99))
    }
    #end.time <- Sys.time()
    #time.taken <- end.time - start.time
    pi_j0_new = param_fit$par[1]
    gamma_new = param_fit$par[2]
    p1_1j_new =  param_fit$par[3]
  } else if(opt_method == "NR"){ # Newton-Raptson for estimating alpha and beta

  } # else if("MoM"){ # methods of moments for alpha and beta
  #}
  return(list(pi_j0_new,gamma_new,p1_1j_new))
}

# opt
EM <- function(X_j,N_j,pi_j0=0.1,gamma=0.5,pi_1j=0.5,iterations=10,stop_diff=1e-05,
               neighbor_s_list,neighbor_s_vaf, Wc,prior_energy=F,
               save_path=NULL,verbose=T){
  # intialize parameter
  init_pi_j0=pi_j0;init_gamma=gamma;init_pi_1j=pi_1j;
  Z_j=sample(c(0,1),size=length(N_j),replace=T);
  pi_j0_list = c()
  gamma_list = c()
  pi_1j_list = c()
  LPP_list = c()
  prev_param = list(0,0,0)
  params = list(pi_j0,gamma,pi_1j)
  start.time0 <- Sys.time()
  for(i in 1:iterations){
    print(paste0("EM iter:",i))
    diff = sum(unlist(prev_param) - unlist(params))
    if(abs(diff) < stop_diff){
      print("difference small enough")
      break
    }
    #start.time <- Sys.time()
    # Estep to get
    prev_param = params
    #delta_ij_hat = Estep(X_sub,N_sub,alpha,beta,delta_gc)
    # get Z identiy from ICM updates
    Z_j = ICM_update(X_j,N_j,Z_j,pi_j0,gamma,pi_1j,neighbor_s_vaf,Wc,neighbor_s_list,
                     neighbor_num=4,icm_iter=10, prior_energy=prior_energy,verbose=verbose)
    # Maximization Step in the EM, estimate each parameter based on expectation of Zij
    params = Mstep(X_j,N_j,Z_j,pi_j0,gamma,pi_1j,
                    neighbor_s_vaf,Wc,
                   prior_energy=prior_energy,opt_method="optim")
    # save all parameters for diagnostic
    pi_j0=params[[1]]
    pi_j0_list = append(pi_j0_list,pi_j0)
    gamma=params[[2]]
    gamma_list = append(gamma_list,gamma)
    pi_1j=params[[3]]
    pi_1j_list = append(pi_1j_list,pi_1j)

    # log posterior
    new_LPP = log_posterior_prob(N_j,X_j,pi_j0,gamma,Z_j,pi_1j,neighbor_s_vaf,Wc,prior_energy=prior_energy)
    LPP_list = append( LPP_list,new_LPP)
    if((verbose==T) & (i%%round(iterations/10)==0)){
      #end.time <- Sys.time()
      #time.taken <- end.time - start.time
      end.time0 <- Sys.time()
      time.taken0 <- difftime(end.time0, start.time0, units='mins') #end.time0 - start.time0
      message(paste0("Iteration till ", i, "th takes ", format(time.taken0,digits=2)," s"))
    }
  }
  # get final identity
  Z_j = ICM_update(X_j,N_j,Z_j,pi_j0,gamma,pi_1j,neighbor_s_vaf,Wc,neighbor_s_list,
                   neighbor_num=4,icm_iter=10, prior_energy=prior_energy,verbose=verbose)

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

