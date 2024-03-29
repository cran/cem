pair <- function(obj, data, method=NULL, mpower=2, verbose=1){
		nm <- NULL
		strataID <- unique(obj$strata)
		reservoir <- NULL	
		
		cl <- unlist(lapply(data, class))
		anyF <- any(cl %in% c("character","factor"))
		if(anyF){
		  cat("\nTransforming factor/charater variable to numeric to calculate distance in pair matching...\n")
		  idf <- which(cl %in% c("character","factor"))
		  for(i in 1:length(idf)){
		    data[,idf[i]] <- as.numeric(as.character(data[,i]))
		  }
		}
		for(i in strataID){
		 idx <- which(obj$strata==i)
		 n <- length(idx)
		 m <- n/2
		 if(m != as.integer(m)) { # we have an odd number of units
	         reservoir <- c(reservoir, idx[1])
		   idx <- idx[-1]  
		   n <- n-1
		   m <- n/2
		  }
		  
		  if(length(idx)>0){
		   if(is.null(method)){
		    mat <- matrix(runif(n*n), n, n)
			} else {
		     mat <- as.matrix(dist(data[idx, obj$vars], method=method, p=mpower))
			}
			colnames(mat) <- rownames(data)[idx]
			rownames(mat) <- rownames(data)[idx]
		   mat1 <- matrix(mat[1:m,(m+1):n], m,m )
		   rownames(mat1) <- (rownames(data)[idx])[1:m]
		   colnames(mat1) <- (rownames(data)[idx])[(m+1):n]
		   if(dim(mat1)[1] == 1){
            nm <- rbind(nm, c(colnames(mat1), rownames(mat1)))
		   } else {
		    for(k in 1:m){
     #        mins <- apply(mat1, 2, function(x) min(x, na.rm=TRUE))
		#	 min.c <- min(mins, na.rm=TRUE)
		#      col <- which.min(mins)
		#      row <- which(mat1[,col]==min.c)[1]
		   if(all(is.na(mat1))) break ;
			 idxMin <- which(mat1 == min(mat1,na.rm=TRUE), arr.ind = TRUE)
			 col <- idxMin[1,"col"]
			 row <- idxMin[1,"row"]
			 mat1[row, 1:m] <- NA
			 mat1[1:m ,col] <- NA
			 nm <- rbind(nm, c(colnames(mat1)[col], rownames(mat1)[row])) 
	        }
		   }
		  } 
		  
		}
#  print(reservoir)
         idx <- reservoir
		 reservoir2 <- NULL
         n <- length(idx)
		 m <- n/2
		 if(m != as.integer(m)) { # we have an odd number of units in the reservoir
	
		   reservoir2 <- idx[1]
		   idx <- idx[-1]  
		   n <- n-1
		   m <- n/2
	}
		  
		 nm2 <- nm
		 if(length(idx)>0){
		   if(is.null(method)){
		    mat <- matrix(runif(n*n), n, n)
			} else {
		     mat <- as.matrix(dist(data[idx, obj$vars], method=method, p=mpower))
			}
			colnames(mat) <- rownames(data)[idx]
			rownames(mat) <- rownames(data)[idx]

		   mat1 <- matrix(mat[1:m,(m+1):n], m,m )
		   rownames(mat1) <- (rownames(data)[idx])[1:m]
		   colnames(mat1) <- (rownames(data)[idx])[(m+1):n]
		   if(dim(mat1)[1] == 1){
            nm2 <- rbind(nm2, c(colnames(mat1), rownames(mat1)))
		   } else {
		    for(k in 1:m){
#              mins <- apply(mat1, 2, function(x) min(x, na.rm=TRUE))
# 			 min.c <- min(mins, na.rm=TRUE)
# 			 col <- which.min(mins)
# 			 row <- which(mat1[,col]==min.c)[1]
			 if(all(is.na(mat1))) break ;
			 idxMin <- which(mat1 == min(mat1,na.rm=TRUE), arr.ind = TRUE)
			 col <- idxMin[1,"col"]
			 row <- idxMin[1,"row"]
			 mat1[row, 1:m] <- NA
			 mat1[1:m ,col] <- NA
			 nm2 <- rbind(nm2, c(colnames(mat1)[col], rownames(mat1)[row])) 
	        }
		   }
		  } 
  		
  paired <- rep(NA, NROW(data))
  names(paired) <- rownames(data)
  full.paired <- rep(NA, NROW(data))
  names(full.paired) <- rownames(data)
  k <- 0
  for(i in 1:NROW(nm)){
   k <- k+1
   paired[nm[i,]] <- k
  }		

  k <- 0
  for(i in 1:NROW(nm2)){
   k <- k+1
   full.paired[nm2[i,]] <- k
  }		
 
  if(verbose>=1)
    cat(sprintf("\nTotal number of units paired in CEM strata: %d\nTotal number of units matched: %d\n", 2*NROW(nm), 2*NROW(nm2)))
  if(length(reservoir2)>0){
   if(verbose>=1)  
     cat(sprintf("Unit corresponding to row `%s', not paired\n", rownames(data)[reservoir2]))
  }
  if(verbose>=1)
    cat("\n")
  return(invisible(list(reservoir=reservoir,reservoir2=reservoir2, paired=paired, full.paired=full.paired)))
}


