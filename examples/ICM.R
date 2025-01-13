# this script contains code for HMRF-ICM process
library(distances)
calc_pairwise_distance <- function(spatial_coords){
  D= distances(spatial_coords)
}
# This function icm updates Zsj (varrier identity)
ICM_update <- function(X_j,N_j,Z_j,pi_j0,gamma,pi_1j,neighbor_s_vaf,Wc,neighbor_s_list, 
                       neighbor_num=4,icm_iter=10, prior_energy=T){
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
        prob_0 = prob_0 * prior_energy_U(s,neighbor_num,neighbor_s_list, Z=0)
        prob_1 = prob_1 * prior_energy_U(s,neighbor_num,neighbor_s_list, Z=1)
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

# calculates P(theta| x_neighbor)? Maybe include in the next version after discussion.
prior_energy_U <- function(s,neighbor_num=4,neighbor_s_list, Z_j,Z){
  total_n = neighbor_num
  neighbor_spots = neighbor_s_list[[s]]
  Vc_sum = sum(1/2 *(1-(Z == Z_j[neighbor_spots])))
  Ux = 1/((2*pi)^(total_n/2)) * exp(-Vc_sum)
  return(U_x)
}

