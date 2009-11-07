plot.imbalance.space <- function(...) imbspace.plot(...)

imbspace.plot2 <- function(obj,group="1"){
	if(class(obj) != "imbalance.space")
	stop("obj must be of class `imbalance.space'")

	g <- sprintf("G%s",group)
	n <-  obj$space[[g]]
	ML1 <- obj$space$ML1
	Relaxed <- obj$space$Relaxed
	class <- rep("relax", length(n))
	class[ which(Relaxed=="<raw>") ] <- "raw"
	class[ which(Relaxed=="<start>") ] <- "start"
	
	tab <- data.frame(n=n, ML1=ML1, class=class)

	print(xyplot(ML1 ~ 1/sqrt(n), group=class, data=tab, xlab="1/sqrt(matched units)", ylab="median of L1 profile", 
	   pch=20, alpha=c(1,0.5,1), col=c("red","gray","green")))
}


imbspace.plot <- function(obj,group="1"){
	if(!interactive()){
	   imbspace.plot2(obj, group)
	   return(NULL)	
	}
	 

	if(class(obj) != "imbalance.space")
		stop("obj must be of class `imbalance.space'")
	g <- sprintf("G%s",group)
	n <-  obj$space[[g]]
	ML1 <- obj$space$ML1
	Relaxed <- obj$space$Relaxed
	class <- rep("relax", length(n))
	class[ which(Relaxed=="<raw>") ] <- "raw"
	class[ which(Relaxed=="<start>") ] <- "start"
	name.vars <- names(obj$coars[[1]])
	n.vars <- length(name.vars)
	tab <- data.frame(n=n, ML1=ML1, class=class)
	 
	
	update.imbplot <- function(mT){
	if(length(mT)==0)
 	 return(FALSE)

	idx <- which(n==n[mT])

	 class <- rep("relax", length(n))
	 class[ which(Relaxed=="<raw>") ] <- "raw"
	 class[ which(Relaxed=="<start>") ] <- "start"
	 class[idx] <- "selected"
	 tab$class <- class

	y <- lapply(obj$coars[idx], function(a) unlist(lapply(a, length))) 
	x <- matrix(unlist(y), length(y), length(y[[1]]), byrow=TRUE) 
	colnames(x) <- names(y[[1]])
	tmp <- as.data.frame(x)

	p1 <- xyplot(ML1 ~ 1/sqrt(n), group=class, data=tab, 
		xlab="1/sqrt(matched units)", ylab="median of L1 profile", 
	   	pch=20, alpha=c(1,0.5,0.5,1), col=c("red","gray","blue","green"),
		main=sprintf("matched units=%d", n[mT]))
	p2 <- parallel( ~tmp, tmp, col=gray(1-ML1[idx]),alpha=0.5, main="right click to exit")
	print(p2, split=c(2,1,2,1))
	print(p1, split=c(1,1,2,1), newpage=FALSE)
	nvars <- dim(tmp)[2]

	tmp2 <- data.frame(tmp, ML1=ML1[idx])	
	rownames(tmp2) <- idx

	old.idx <<- list(data=tmp2, coars=obj$coars[idx])
	return(TRUE)

       }

	old.idx <- NULL

	goOn <- update.imbplot( 1 )

	while(goOn){
	 trellis.focus("panel", 1, 1)
	 goOn <- update.imbplot(panel.identify(n=1)) 
	}
	
	class(old.idx) <- "selected.cem"
  	old.idx
	
}

print.selected.cem <- function(x, ...){
	print(x$data)

}


