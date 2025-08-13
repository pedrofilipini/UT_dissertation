
## BART: Bayesian Additive Regression Trees
## Copyright (C) 2017-2018 Robert McCulloch and Rodney Sparapani

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, a copy is available at
## https://www.R-project.org/Licenses/GPL-2

plbart = function(
    x.test,		#x matrix to predict at
    model,		#lbart model
    mc.cores=1L,         #thread count
    dodraws=TRUE,
    nice=19L             #mc.plbart only
)
{
  mu = model$ymean #mean of y to be added at the end
  
  if((model$linear == "control")||(model$linear == "both")){
    res_con = .Call("cplbart",
                    model$treedraws_con,	#trees list
                    t(x.test),      #the test x
                    t(cbind(rep(1,nrow(x.test)),x.test)),
                    mc.cores   	#thread count
    )
  }else{
    res_con = .Call("cpwbart",
                    model$treedraws_con,	#trees list
                    t(x.test),      #the test x
                    mc.cores   	#thread count
    )
  }
  
  if((model$linear == "moderate")||(model$linear == "both")){
    res_mod = .Call("cplbart",
                    model$treedraws_mod,	#trees list
                    t(x.test),      #the test x
                    t(cbind(rep(1,nrow(x.test)),x.test)),
                    mc.cores   	#thread count
    )
  }else{
    res_mod = .Call("cpwbart",
                    model$treedraws_mod,	#trees list
                    t(x.test),      #the test x
                    mc.cores   	#thread count
    )
  }
  
  if((dim(res_con$yhat)[1]) == length(model$eta_con)){
    for(i in 1:dim(res_con$yhat)[1]){
      res_con$yhat[i,] = res_con$yhat[i,]*(model$eta_con[i])
    }
  }
  
  if((dim(res_mod$yhat)[1]) == length(model$eta_mod)){
    for(i in 1:dim(res_mod$yhat)[1]){
      res_mod$yhat[i,] = res_mod$yhat[i,]*(model$eta_mod[i])
    }
  }
  
  if(dodraws) return(list(mu_post = res_con$yhat+mu, b_post = res_mod$yhat))
  else return(list(mu_post = apply(res_con$yhat, 2, mean)+mu, b_post = apply(res_mod$yhat, 2, mean)))
}
