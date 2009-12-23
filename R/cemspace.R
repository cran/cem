`cemspace` <-
function (treatment=NULL, data = NULL, R=100, grouping = NULL, drop=NULL,
L1.breaks = NULL, plot = TRUE, fixed = NULL, minimal = 1, maximal = 5,
M=250, raw.profile=NULL) 
{
    if (!is.null(raw.profile) & class(raw.profile) != "L1profile") 
	 stop("raw.profile must be of class `L1profile'")

    gn <- NULL
    if (!is.null(grouping) & !is.null(names(grouping))) {
        gn <- names(grouping)
        n.gn <- length(gn)
        for (g in 1:n.gn) {
            if (!is.null(data)) 
			data[[gn[g]]] <- group.var(data[[gn[g]]], grouping[[g]])
        }
    }

	imb0 <- NULL
	medianL1 <- NULL
	medianCP <- NULL	

	drop <- unique(drop)
	dropped <- match(drop, colnames(data))
	dropped <- dropped[!is.na(dropped)]
	
	if(length(dropped)>0) 
	 data <- data[-dropped]
	vnames <- colnames(data)

	if(!is.null(treatment)){
		groups <- as.factor(data[[treatment]])
        idx <- match(treatment, colnames(data))
		if(length(idx)>0)
		 vnames <- vnames[-idx]
	}
	
	if(is.null(raw.profile)){
		cat("\nCalculating L1 profile for the raw data...\n")
		imb0 <- L1.profile(groups, data, drop=treatment, M=M, plot=FALSE)
		medianL1 <- median(imb0$L1)
		medianCP <- imb0$CP[[ which(imb0$L1>medianL1)[1] ]]
	} else {
		imb0 <- raw.profile
		medianL1 <- raw.profile$medianL1
		medianCP <- raw.profile$medianCP
	}

	mnames <- vnames
	if(!is.null(gn)){
	 idx <- match(gn, vnames)
	 if(length(idx)>0)
	  mnames <- mnames[-idx] 
	}

	if (!is.null(fixed)) {
        idx <- match(fixed, vnames)
        if(length(idx) > 0) 
		 mnames <- mnames[-idx]
    }
	
    nv <- length(mnames)
    v.num <- 1:nv

    b.seq <- vector(nv, mode = "list")
    names(b.seq) <- mnames
	tmp.min <- 2
	tmp.max <- 7
	if(!is.list(minimal)){
	 tmp.min <- minimal+1
	 minimal <- vector(nv, mode="list")
	 for (i in v.num) 
	  minimal[[i]] <- tmp.min
	}
	if(!is.list(maximal)){
		tmp.max <- maximal+1
		maximal <- vector(nv, mode="list")
		for (i in v.num) 
		maximal[[i]] <- tmp.max
	}
    for (i in v.num) {
		vna <- mnames[i]
		min.br <- tmp.min
		max.br <- tmp.max
		
		nuval <- length(unique(data[[vna]]))

        if( (nuval==2) | !(is.numeric(data[[vna]]) | is.integer(data[[vna]]) | is.logical(data[[vna]])) )
			max.br <- nuval+1
		
		if (!is.null(minimal[[vna]]))  
			min.br <- minimal[[vna]] + 1
		min.br <- max(tmp.min, min.br)
		if (!is.null(maximal[[vna]])) 
			max.br <- maximal[[vna]] + 1
		max.br <- min(tmp.max, max.br)
		b.seq[[i]] <- min.br:max.br
    }
    relax <- NULL

	
	g.names <- levels(groups)
	n.groups <- length(g.names)
	
    cat(sprintf("Executing %d different random CEM solutions\n", R))

    tab <- as.data.frame(matrix(NA, R + 1, 2 * n.groups + 2))
    colnames(tab) <- c(paste("G", g.names, sep = ""), 
					   paste("PercG", g.names, sep = ""), "ML1", "Relaxed")
	n.coars <- dim(tab)[1]
	coars <- vector(n.coars, mode="list")

	tab[1, 1:n.groups] <- as.numeric(table(groups))
    tab[1, (n.groups + 1):(2 * n.groups)] <- 100
    tab$Relaxed[1] <- "<raw>"
	coars[[1]] <- NULL
	tab[1, "ML1"] <- medianL1

    
	newcut <- vector(nv, mode="list")
    names(newcut) <- mnames

	pb <- txtProgressBar(min = 1, max = R, initial = 1, style = 3)

    for (r in 1:R) {
		setTxtProgressBar(pb, r)
		for(i in 1:nv)
			newcut[[i]] <- sample(b.seq[[i]], 1) 
		obj <- cem(treatment, data, cut=newcut, eval=FALSE)
		
		coars[[r+1]] <- obj$breaks
		tab[r+1, 1:(2 * n.groups)] <- as.numeric(c(obj$tab[2,], 
												   obj$tab[2, ]/obj$tab[1, ] * 100))
		tab$Relaxed[r+1] <- "random"
		tab[r+1, "ML1"] <- L1.meas(groups, data, drop=treatment, breaks=medianCP, 
								    weights=obj$w)$L1

    }
	close(pb)
    idx <- order(tab[, "ML1"])
    tab <- tab[idx, ]
    rownames(tab) <- 1:(dim(tab)[1])
	out <- list(space = tab)
	out$L1breaks <- L1.breaks
	out$raw.profile <- imb0
	out$tab <- obj$tab
    out$medianCP <- medianCP
	out$medianL1 <- medianL1
	out$coars <- coars[idx]
	out$n.coars <- n.coars
	out$match <- obj
    class(out) <- "imbalance.space"
	
    if (plot) 
	  plot(out,data=data,explore=interactive())
	
    return(invisible(out))
}


