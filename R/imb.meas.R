L1.meas <- function(group, data, drop=NULL, breaks=NULL, weights){
  if(!is.null(drop)){
   drop <- unique(drop)
   dropped <- match(drop, colnames(data))
   dropped <- dropped[!is.na(dropped)]

   if(length(dropped)>0) 
		data <- data[-dropped]
  }

 if(is.null(breaks)){
  vars <- colnames(data)
  nv <- length(vars)
  breaks <- vector(nv, mode="list")
  for(i in 1:nv){
   if(is.numeric(data[i]) | is.integer (data[i]))
    breaks[[i]] <- pretty(range(data[i],na.rm=TRUE), n=nclass.scott(data[[i]]), 1)
   names(breaks) <- vars
  }
 }  

 n <- dim(data)[1]

 if(missing(weights)){
  weights <- rep(1, n)
 }
 rem <- which(weights<=0)
 if(length(rem)>0){
  data <- data[-rem,]
  weights <- weights[-rem]
  group <- group[-rem]
 }


 cem.imbalance(group, data, collapsed=FALSE, reduced=FALSE, breaks=breaks, weights=weights)
}


cem.imbalance <- function(group, data, collapsed=TRUE, reduced=TRUE, breaks=NULL, weights){
 lv <- unique(na.omit(group))


 if(length(lv)>2)
  stop("more than 2 groups")
   
 if(is.null(breaks)){
  vars <- colnames(data)
  nv <- length(vars)
  breaks <- vector(nv, mode="list")
  for(i in 1:nv){
   if(is.numeric(data[i]) | is.integer (data[i]))
    breaks[[i]] <- pretty(range(data[i],na.rm=TRUE), n=nclass.scott(na.omit(data[[i]])), 1)
   names(breaks) <- vars
  }
 }  
 
 if(!reduced){
  tmp <- reduce.data(data, collapse=TRUE, breaks=breaks)
  data <- tmp$data
  new.breaks <- tmp$breaks
  collapsed <- TRUE
 }  
 if(!collapsed)
  data <- collapse.data(data)
  
 n <- length(data)
 if(missing(weights)){
  weights <- rep(1, n)
 }

 idx1 <- which(group==lv[1])
 idx2 <- which(group==lv[2])

 n1 <- sum(weights[idx1])
 n2 <- sum(weights[idx2])

 tab <- table(data)
 cell <- attr(tab, "dimnames")[[1]]
 L1 <- 0

 for(i in cell){
  idx <- as.numeric(which(data == i))   
  jdx1 <- intersect(idx, idx1)
  jdx2 <- intersect(idx, idx2)
  m1 <- sum(weights[jdx1])/n1
  m2 <- sum(weights[jdx2])/n2
  L1 <- L1 + abs(m1-m2) 
 }
 out <- list(L1=L1, breaks=new.breaks)
 class(out) <- "L1.meas"
 return(out)
}


imbalance <- function(group, data, drop=NULL, breaks=NULL, weights){
 if (!is.data.frame(data))
        stop("Data must be a dataframe", call. = FALSE)

  if(!is.null(drop)){
   drop <- unique(drop)
   dropped <- match(drop, colnames(data))
   dropped <- dropped[!is.na(dropped)]

   if(length(dropped)>0) 
		data <- data[-dropped]
  }
 if(is.null(breaks)){
  vars <- colnames(data)
  nv <- length(vars)
  breaks <- vector(nv, mode="list")
  for(i in 1:nv){
   if(is.numeric(data[i]) | is.integer (data[i]))
    breaks[[i]] <- pretty(range(data[i],na.rm=TRUE), n=nclass.scott(na.omit(data[[i]])), 1)
   names(breaks) <- vars
  }
 }  

 n <- dim(data)[1]
 if(missing(weights)){
  weights <- rep(1, n)
 }
 rem <- which(weights<=0)
 if(length(rem)>0){
  data <- data[-rem,]
  weights <- weights[-rem]
  group <- group[-rem]
 }



 lv <- unique(na.omit(group))
 if(length(lv)>2)
  stop("more than 2 groups")
 idx1 <- which(group==lv[1])
 idx2 <- which(group==lv[2])
 tmp <- reduce.data(data, breaks=breaks)$data
 
 stdmean <- function(x, wh) weighted.mean(x, w=wh, na.rm=TRUE)

 qt.w <-function(x, probs=c(0,0.25,0.5, 0.75, 1), wh){
  q <- NULL
  which(is.na(x)) -> id
  if(length(id)>0)
   x <- x[-id]
  ord <- order(x)
  x <- x[ord]
  wh <- wh[ord]
  F <- cumsum(wh)/sum(wh)
  
  for(i in probs){
    q <- c(q, x[which(F>=i)[1]])  
  }
  idx <- which(is.na(q))
  q[idx] <- max(x)
  q
}

 vnames <- colnames(data)
 nv <- length(vnames)
 tab <- as.data.frame(matrix(NA, nv, 8)) 
 rownames(tab) <- vnames
 colnames(tab) <- c("statistic", "type","L1", "min", "25%",  "50%",  "75%", "max")

 for (i in 1:nv){
  tab[i,3] <- L1.meas(group, tmp[i], breaks=breaks,weights=weights)$L1
  if((is.numeric(data[,i]) | is.integer(data[,i])) ){
    tab[i,1] <- stdmean(data[idx1,i],weights[idx1]) - stdmean(data[idx2,i],weights[idx2])  
    tab[i,4:8] <- qt.w(x=data[idx1,i], wh=weights[idx1]) - qt.w(x=data[idx2,i], wh=weights[idx2])  
    tab[i,2] <- "(diff)"
  } else {
   t1 <- table(data[idx1,i])
   t2 <- table(data[idx2,i])
   keep <- which(t1>0 & t2>0)
    tab[i,1] <- as.numeric(chisq.test(cbind(t1[keep],t2[keep]))$statistic)   
    tab[i,4:8] <- NA  
    tab[i,2] <- "(Chi2)"
  }
 }
 
 out <- list(tab=tab, L1=L1.meas(group=group, data=data, breaks=breaks, weights=weights))
 class(out) <- "imbalance"
 return( out )
 
}


print.imbalance <- function(x,...){
 print(x$L1)
 cat("Univariate Imbalance Measures:\n\n")
 print(x$tab,...)
}


print.L1.meas <- function(x,...){
 cat(sprintf("\nMultivariate Imbalance Measure: L1=%f\n\n", x$L1))
}
