`cem.summary` <-
function (obj, verbose=0) 
{
    if(is.null(obj$matched))
	 return(NULL)
	 
    tab <- matrix(,3,obj$n.groups)
	tab[1,] <- table(obj$groups)
	tab[3:2,] <- table(obj$matched, obj$groups)
	tab[3,] <- tab[1,] - tab[2,]

    colnames(tab) <- paste("G", obj$g.names,sep="")
	rownames(tab) <- c("All","Matched","Unmatched")
	if(verbose > 1){
	 cat(sprintf("\nCEM Subclasses: %d\n", length(obj$mstrataID) ))
	 cat("\nSample sizes:\n")
	 print(tab)
	}
	return(tab)
}

