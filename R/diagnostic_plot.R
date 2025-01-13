#'
#' Helper functions of the spatial MT model
#'
#' `plot_vaf_cellprop()` plots the VAF vs coverage for a given celltype and given variant.
#'
#' @param i A numeric indicator of variant.
#' @param voi A string vector containing names of variant of interest. Usually a subset of variants that went through basic filtering.
#' @param af.dm A MT variant x cell variant allele frequnecy (VAF) matrix. Equivalent as vaf matrix.
#' @param norm_weights A cell x celltype weight matrix, representing spatial decomposition results from tools like RCTD.
#' @param plot_df A plotting dataframe.
#'
#' @import dplyr ggplot2 grid gridExtra
#' @export
# plotting VAF vs celltype proportions, colored by coverage
plot_vaf_cellprop <- function(i,af.dm,norm_weights,plot_df,intercept,coef,pval,permute=F){
  if(permute==T){
    p_str = "permute_p"
  }else{
    p_str = "ANOVA_p"
  }
  # specify a variant
  var = rownames(af.dm)[i]
  p_list = list()
  for(k in 1:length(colnames(norm_weights))){
    celltype=colnames(norm_weights)[k]
    fit_df=data.frame(x=seq(0,1,0.01),y=intercept[k]+coef[k]*seq(0,1,0.01))
    p_list[[celltype]]  = local({
      celltype=celltype
      var=var
      var_cov=paste0(var,"_cov")
      temp=plot_df[order(plot_df[,var],decreasing=F),]
      #ggplot(temp %>% arrange(get(var_cov)),
             #aes(x=get(celltype),y=get(var),color=get(var_cov))) + #,size=get(var_cov))) +
      ggplot(temp %>% arrange(get(var_cov))) +
        geom_point(size=1,alpha=0.5,aes(x=get(celltype),y=get(var),color=get(var_cov))) +
        geom_line(data=fit_df,aes(x=x,y=y)) +
        theme_bw()+
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.text.x =  element_text(angle = 90,vjust=0.5,hjust=1))+
        scale_color_continuous(low="grey",high="red") +
        ggtitle(paste0(celltype," ",p_str,": ",formatC(pval[k],format="e",digits=2))) +
        xlab("Celltype Proportions") +  ylab("VAF") +
        guides(color=guide_legend(title="Coverage"))#,size=guide_legend(title="Coverage"))
    })

  }
  # to speed up, can choose to save in pdf
  grid.arrange(grobs=p_list,top=var)
}
