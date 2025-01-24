#' Regression Test for Identifying Significant Mutations and Carrier Celltypes.
#'
#' `celltype_test()` implements the statistical tests for identifying significant celltypes given spatial decomposition and varaint information.
#'
#' @param celltypes A string vector containing all unique celltypes.
#' @param voi A string vector containing names of variant of interest. Usually a subset of variants that went through basic filtering.
#' @param N A MT variant x cell total coverage matrix.
#' @param X A MT variant x cell alternative allele count matrix. Optional if vaf matrix is provided.
#' @param vaf A MT variant x cell variant allele frequnecy (VAF) matrix. Equivalent as af.dm matrix.
#' @param Ws A cell x celltype weight matrix. Representing spatial decomposition results from tools like RCTD.
#' @param spatial_coords A cell x spatial coordinates matrix. Representing the actual spot/cell location in space.
#' @param test_type Type of Celltype test. Options from c("linear","weighted").
#' @param permute_num A numeric number representing permuation number for each variant.
#'
#' @param ... Other options used to control matching behavior between duplicate
#'   strings. Passed on to [stringi::stri_opts_collator()].
#' @returns A list object containing fitted parameters (intercept, coefficient and p-values).
#' @import Matrix
#'
#' @examples
#' celltype_test(celltypes=c("BE","DYSP","Normal"), voi=voi_list,
#' N=N,vaf=af.dm, test_type="linear",permute_num=1000)
#'
#' @export
#'
celltype_test <- function(celltypes=NULL,voi=NULL,N=NULL,vaf=NULL,X=NULL,Ws=NULL,spatial_coords=NULL,
                          test_type=c("linear","weighted"),permute_num=1000,plot=F){
  # check for input and conditions
  if(test_type == "linear"){
    message("Performing Linear Rgression ...")
  }else if(test_type == "weighted"){
    message("Performing Linear Rgression weighted by coverage...")
  }else{
    error("Not a valid test. Please choose from documented tests.")
  }

  if(is.null(vaf)){
    if(!is.null(N)){
      vaf_j = X[j,]/N[j,]
    }else{
      error("Missing alternative allelle count matrix N or frequency matrix vaf.")
    }
  }
  # check if celltypes aligned with weight matrix
  if(is.null(celltypes)){
    message("Missing celltypes...Set as the column names of celltype weight matrix Ws.")
    celltypes = colnames(Ws)
  }else if(any(celltypes != colnames(Ws))){
    error("Celltypes and the column names of celltype weight matrix Ws not matched...Please check the variables celltypes or Ws.")
  }

  # create plotting data frame
  if(plot==T){
    intersect_bc=intersect(intersect(colnames(N),rownames(Ws)),rownames(spatial_coords))
    plot_df=cbind(spatial_coords[intersect_bc,],
                  Ws[intersect_bc,],
                  as.data.frame(t(vaf[,intersect_bc])))
    # add in coverage of variant of interest
    cov = N[,intersect_bc];rownames(cov) =  paste0(rownames(N),"_cov")
    plot_df = cbind(plot_df,t(cov))
  }

  # for each variant j, fit a Xs/Ns vs Ws Linear regression
  intercept_df = Matrix(NA,nrow=length(voi),ncol=length(celltypes))
  coef_df = Matrix(NA,nrow=length(voi),ncol=length(celltypes))
  pval_df=Matrix(NA,nrow=length(voi),ncol=length(celltypes))

  if(plot==T){
    pdf(paste0(test_type,"_regression_anova_p_VAF_vs_celltype_plot.pdf"),height=10,width=10)
  }
  for(j in 1:length(voi)){
    var = voi[j]
    print(paste0(j,": ", var))
    N_j = N[j,]
    vaf_j = vaf[j,]


    for(k in 1:length(celltypes)){ # for each celltype

      data = data.frame(vaf_j=vaf_j,N_j,W_sk=Ws[,k])
      cell_ratio_bins = seq(from=0,to=1,by=0.05)

      if(test_type == "linear"){
        fit=lm(vaf_j ~ W_sk,data=data)# weights equals to inverse of variance (1/sigma_i^2)
      }else if(test_type == "weighted"){
        fit=lm(vaf_j ~ W_sk,data=data, weights=sqrt(N_j))# weights equals to inverse of variance (1/sigma_i^2)
      }
      res=summary(fit)
      intercept_df[j,k]=res$coefficients[1,1]
      coef_df[j,k]= res$coefficients[2,1]

      if(test_type == "linear"){
        pval_df[j,k] = res$coefficients[2,4]
      }else if(test_type == "weighted"){
        # doing permutations and see where the data locates
        set.seed(42)
        coef_list <- c()
        #  decouple VAF with Ws
        for(l in 1:permute_num){
          idx=sample(1:length(vaf_j),size=length(vaf_j))
          vaf_shuffled = vaf_j[idx]
          #vaf_shuffled = pmin(X_shuffled/(N_j+1e-6),1) # cap max VAF to be by 1
          fit=lm(vaf_shuffled ~ Ws[,k],weights=sqrt(N_j[idx]))# weights equals to inverse of variance (1/sigma_i^2)
          res_shuf=summary(fit)
          coef_list <- append(coef_list,res_shuf$coefficients[2,1])
        }
        pval_df[j,k] = 1- sum(res$coefficients[2,1] >= coef_list)/length(coef_list)
      }
    }
    # for variant j, plot out the fitted values on the diagnostic plots
    if(plot==T){
      intercept=intercept_df[j,]
      coef=coef_df[j,]
      pval=pval_df[j,]
      plot_vaf_cellprop(j,af.dm,Ws,plot_df,intercept,coef,pval)
    }
  }
  if(plot==T){
    dev.off()
    message("Diagnostic plots generated.")
  }

  return(res=list(intercept=intercept_df,coef=coef_df,pval=pval_df))
}


#' Power Analysis for Identifying Negative Mutations and Carrier Celltypes.
#' Simulate alternative allele count values Xi through Binomial distribution Xi ~ Binom(Ni,beta_0 + beta_1 * Wik)
#'
#' @param beta A numeric value. Equivalent as beta_1 in the simulation.
#' @param Nj A total coverage counts for a given variant j
#' @param X A MT variant x cell alternative allele count matrix. Optional if vaf matrix is provided.
#' @param vaf A MT variant x cell variant allele frequnecy (VAF) matrix. Equivalent as af.dm matrix.
#' @param Ws A cell x celltype weight matrix. Representing spatial decomposition results from tools like RCTD.
#' @param spatial_coords A cell x spatial coordinates matrix. Representing the actual spot/cell location in space.
#' @param test_type Type of Celltype test. Options from c("linear","weighted").
#' @param permute_num A numeric number representing permuation number for each variant.
#'
calc_power <- function(beta,Nj,Wj,vaf_j,alpha,n_sim,verbose=T,report_frac=3){
  rejections = 0
  p_vals = c()
  for(i in 1:n_sim){
    if(verbose){
      if(i%%(round(n_sim/report_frac))){
       print(paste0("sim:",i,"/",n_sim))
      }
    }
    X_sim = rbinom(n=length(Nj),size=Nj,prob=beta*Wj)
    data = data.frame(vaf_j=X_sim/(Nj+0.001),Nj,Wj=Wj)
    fit=lm(vaf_j ~ Wj,data=data)
    res=summary(fit)
    p_vals = append(p_vals,res$coefficients[2,4])
  }
  rejections  = sum(p_vals <= alpha,na.rm=T)
  power = rejections/(n_sim)
  return(power)
}

#' Power Analysis for Identifying Negative Mutations and Carrier Celltypes.
#' for all variants of interest & all celltypes
#' set permute_num=20 for now, shall be 1000
#' @export
power_analysis_all <- function(voi=NULL,celltypes=NULL,Ws=NULL,N=NULL,vaf=NULL,X=NULL,
                               sample_num=100,alpha=0.05,n_sim=10,beta_threshold=0.5,plot=T){
  message("Running Power Analysis for Negative Celltypes and Mutations")
  intersect_bc = intersect(colnames(N),rownames(Ws))
  beta_list = c()

  if(plot){
    pdf("example_power_analysis.pdf",height=8,width=8)
    par(mfrow = c(ceiling(sqrt(length(celltypes))), ceiling(sqrt(length(celltypes))))) # Create a 2 x 2 plotting matrix
  }
  for(j in 1:length(voi)){
    var = voi[j];
    print(paste0("voi ",j,":",var))
    N_j = N[j,intersect_bc]
    vaf_j = vaf[j,intersect_bc]

    for(k in 1:length(celltypes)){
      sample_idx = sample(intersect_bc,sample_num)
      Nj = N_j[sample_idx] # sampled coverage
      Wj = Ws[sample_idx,k] # celltype ratio
      beta_grid <- seq(0.1,1,by=0.2)

      power_results = sapply(beta_grid,function(beta){
        #print(beta)
        calc_power(beta,Nj,Wj,vaf_j,alpha,n_sim,verbose=F)
      })
      names(power_results) = beta_grid
      beta_list[paste0(voi[j],"_",celltypes[k])] <- power_results

      if(plot){
        plot(beta_grid,power_results,type="l",col="blue",lwd=2,
             xlab="Beta",ylab="Power",main=paste0(var,", ",celltypes[k]),
             ylim=c(-0.1,1),xlim=c(0,1)
        )
        points(c(0.1,0.5,1),power_results[as.character(c(0.1,0.5,1))],pch=16)
        abline(h=0.8,col="orange",lty=2)
      }
    }
  }
  if(plot){
    dev.off()
  }
  return(beta_list)

}
