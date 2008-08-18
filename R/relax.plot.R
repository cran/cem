plot.relax.cem <- function(...) relax.plot(...)

relax.plot <- function(tab, group="1", max.terms=50, perc=.5, unique=FALSE, colors=TRUE){
     if(class(tab) != "relax.cem")
	  stop("tab must be of class `relax.cem'")
     par(xpd=TRUE)
	 par(mar = c(10, 5, 4, 4) + 0.1)
     g <- sprintf("G%s",group)
	 pg <- sprintf("PercG%s",group)
     n <- dim(tab[[g]])[1]
	 
	 idx.start <- which(tab[[g]]$Relaxed=="<start>")
	 n.start <- tab[[g]][[g]][idx.start]
	 p.start <- tab[[g]][[pg]][idx.start]
	 title <- sprintf("Pre-relax: %d matched (%.1f %%)", n.start, p.start) 
    
	 idx <- 1:n

	 if(unique){
      idx <- match(unique(tab[[g]][[pg]]), tab[[g]][[pg]])
	  idx1 <- match(p.start, tab[[g]][[pg]][idx])  
	  idx[ idx1 ] <- idx.start
	 }
	 G1 <-  tab[[g]][[g]][idx]
	 PercG1 <-  tab[[g]][[pg]][idx]
	 Relaxed <- tab[[g]]$Relaxed[idx]
	 L1 <- tab[[g]]$L1[idx]
	 var <-  tab[[g]]$var[idx]
	 var <- as.integer(var)
	 n.var <- max(unique(var))
     mycol <- rep(c("black","blue", "green", "yellow", "orange","red","violet", "salmon", "purple","brown","tan4"),10)
     var <- mycol[var] #rainbow(n.var, start=0.2, end=0.8)[var]
     if(!colors)
	  var <- rep("black", length(var))
	   
	 n <- length(idx)
	 
	 
	 p.terms <- which(PercG1>=perc*100)
	
	 max.terms <- min(length(p.terms), max.terms)
	 max <- min(max.terms, n)
     labx <- Relaxed
	
	
	 plot(1:max, PercG1[(n-max+1):n], type="b",axes=F,xlab="",ylab="",pch=16, main=title, col=var[(n-max+1):n])

     axis(1, 1:max, labels=F)    

     y1 <- sort(unique(PercG1[(n-max+1):n]))
	 y2 <- sort(unique(G1[(n-max+1):n]))

     axis(2, y1, labels=F)
     text(rep(par("usr")[1], length(y1)), y1, labels=sprintf("%.1f",y1), srt=0, adj=2,cex=0.75)

     axis(4, y1, labels=F)
     text(rep(par("usr")[2], length(y1)), y1, y2, srt=0, adj=-1,cex=0.75)
 
     text(1:max, par("usr")[3] - 0.35, srt = 45, adj = 1,
          labels = labx[(n-max+1):n], xpd = TRUE, cex=0.75, col=var[(n-max+1):n])
 
     text(1:max,  PercG1[(n-max+1):n], adj=-0.35, 
      labels= sprintf("%3.2f", L1[(n-max+1):n]), srt=-45, cex=0.75)
 
     mtext("number of matched", 4, line=3)
	 mtext("% matched",2,line=4)
     box()

}
