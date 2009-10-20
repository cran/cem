# ref.cut : number of cutpoints of other variables when cutting on one variable
# max.cut : maximal number of cut points
# min.cut : minimal number of cut points

L1.profile <- function (group, data, drop = NULL, min.cut=2, max.cut=6, ref.cut=3, 
 weights, plot=TRUE, add=FALSE, col="red", lty=1) 
{

	if(min.cut<2)
	 min.cut <- 2
	min.cut <- as.integer(min.cut)
	max.cut <- as.integer(max.cut)
	ref.cut <- as.integer(ref.cut)
	
	if (!is.null(drop)) {
        drop <- unique(drop)
        dropped <- match(drop, colnames(data))
        dropped <- dropped[!is.na(dropped)]
        if (length(dropped) > 0) 
		data <- data[-dropped]
    }
	
    vnames <- colnames(data)
    nv <- length(vnames)
	idx <- sapply(1:nv, function(x) is.numeric(data[[x]]) | is.integer(data[[x]])) 
	numvar <- vnames[which(idx)]
	catvar <- vnames[which(!idx)]
	
	R <- sapply(numvar, function(x) range(data[x], na.rm=TRUE),simplify=FALSE)

	br <- vector(length(numvar), mode="list")
	names(br) <- numvar
	out <- NULL
	cat("\nprocessing all variables together")
	for(i in catvar)
	 br[i] <- NULL
	for(theta in min.cut:max.cut){
		cat(".")
		for(i in numvar){
			br[[i]] <- unique(seq(R[[i]][1], R[[i]][2], length=theta))
		}
		out <- c(out, L1.meas(group = group, data = data, breaks = br, weights = weights)$L1)
	}
	names(out) <- sprintf("all(%d)", min.cut:max.cut)

	
	for(j in numvar){
		cat(sprintf("\nprocessing variable(%3.0f%%): %s...",100*which(numvar==j)/length(numvar),j))
		nm.out <- names(out)
		for(i in numvar)
			br[[i]] <- unique(seq(R[[i]][1], R[[i]][2], length=ref.cut))
		for(theta in min.cut:max.cut){
			cat(".")
			br[[j]] <- unique(seq(R[[j]][1], R[[j]][2], length=theta))
			out <- c(out, L1.meas(group = group, data = data, breaks = br, weights = weights)$L1)
		}
		names(out) <- c(nm.out,sprintf("%s(%d)",j,min.cut:max.cut))
	}
	
	class(out) <- "L1profile"
	if(plot)
		plot(out, lty=lty, col=col, add=add)

	return(invisible(out))
}


plot.L1profile <- function(x, add=FALSE, cex.axis=0.3, ...)
{
	if(!add){
		plot(unclass(x),type="n",axes=F,xlab="cut points",ylab=expression(L[1]),ylim=c(0,1))
		axis(1, 1:length(x), names(x), padj=0.5, las=3, xpd=TRUE, cex.axis=cex.axis)
		axis(2)
	}
	lines(x,...)
}
