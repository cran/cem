`att` <-
function (obj, formula, data, model="lm",family="binomial") 
{
 if((class(obj)[1] != "cem.match") && (class(obj)[1] != "multicem"))
  stop("Argument `obj' must be a `cem.match' or `multicem' object") 
 if(!(model %in% c("lm","glm")))
  stop("please, only use `lm' or `glm'")
   
 if(class(obj)[1] == "cem.match"){
  if(model=="lm")
   out <- do.call(model, list(formula=formula, data=data, weights=obj$w))
  else
   out <- do.call(model, list(formula=formula, data=data, weights=obj$w,family=family))
  aa <- summary(out)
  out <- list(model = t(aa$coefficients))
  class(out) <- "cem.att"
  return(out)
 }
 
 if(class(data) != "list")
  stop("Argument `data' must be a list of `data.frame's")

 n.cems <- length(obj) - 1

 if(length(data) != n.cems)
  stop("lengths of `multicem' object and `data' do not match")
  
 est <- vector(n.cems, mode="list")
 for(i in 1:n.cems){ 
  if(model=="lm")
   out <- do.call(model, list(formula=formula, data=data[[i]], weights=obj[[i]]$w))
  else
   out <- do.call(model, list(formula=formula, data=data[[i]], weights=obj[[i]]$w,family=family))
  aa <- summary(out)
  est[[i]] <- aa$coefficients 
 }
 qoi <- numeric(dim(est[[1]])[1])
 seq <- numeric(dim(est[[1]])[1])
 for(i in 1:n.cems){
  qoi <- qoi + est[[i]][,1]
  seq <- seq + est[[i]][,2]^2
 }
 
 qoi <- qoi/n.cems
 seq <- seq/n.cems

 s2 <- numeric(dim(est[[1]])[1])
 for(i in 1:n.cems){
  s2 <- s2 + (est[[i]][,1] - qoi)^2
 }
 s2 <- s2/(n.cems-1)

 S <- sqrt(seq + s2*(1+1/n.cems))

 mat <- rbind(qoi, S)
 rownames(mat) <- c("Estimate", "Std. Error")
 
 out <- list(mult=est, model = mat) 
 class(out) <- "cem.att"
 out
}


print.cem.att <- function(x,...){
 if(!is.null(x$model)){
  cat("\n")
  print(x$model)
  }
}

