prepData_DCM <-
function(M1, M2, QN = FALSE, resid.full = FALSE, echo = FALSE, strict = "low"){

	# Make sure data is numeric
	
	M1 = t(apply(M1, 1, as.numeric))
	M2 = t(apply(M2, 1, as.numeric))
	
	
	# Simplify
			
	M1 = (M1 - min(M1, na.rm = TRUE))/max(M1, na.rm = TRUE)
	M2 = (M2 - min(M2, na.rm = TRUE))/max(M2, na.rm = TRUE)
		
	M1 = round(M1, 5)
	M2 = round(M2, 5)
	
	# Make sure same number of rows
	p1 = nrow(M1)
	p2 = nrow(M2)
	
	if(p1 != p2){
		print("Datasets must have the same number of rows.")
		break
	}else{
		p = p1
	}
	
	# Make sure 3 or more samples
	n1 = ncol(M1)
	n2 = ncol(M2)
	
	if(n1 < 3 | n2 < 3){
		print("Must have at least 3 samples in each group.")
		break
	}
			
	# Remove and report bad rows
	
	killrows = sanitize_DCM(M1, M2, strict)

	if(length(killrows) > 0){
		if(echo){print(sprintf("Removing %i rows - too many zero expressions.", length(killrows)))}
		
		# Record dimension
		idcs = (1:p)[-killrows]

		# Remove ignored rows
		M1 = M1[-killrows,]
		M2 = M2[-killrows,]	
	}else{
		idcs = 1:p
	}
		
	p = length(idcs)
	
	# If not quantile normalizing, check normality and warn if non-Gaussian
	if(!QN){
		
		rand = sample((1:p), min(1000, p))
		pvs1 = sapply(rand, function(x) shapiro.test(M1[x,])$p.value)
		pvs2 = sapply(rand, function(x) shapiro.test(M2[x,])$p.value)
		
		bad1 = length(bh(pvs1)) > 0
		bad2 = length(bh(pvs2)) > 0
		
		if((bad1 | bad2) & echo){
			print("Warning: Data is non-normal; consider quantile normalizing.  (Set QN = TRUE)")
		}

	}
		
	print("Done checking data.")
	
	################# Process data ###################
	
	# Quantile normalize
	if(QN){
		M1 = quantNorm(M1)
		M2 = quantNorm(M2)
	}
	
	# Standardize (for easier calculation of correlations)
	M1 = stdize(M1)
	M2 = stdize(M2)
	
	m1 = apply(M1, 2, mean)
	m2 = apply(M2, 2, mean)
	
	c1 = (sum(m1^2)*(n1^2) - n1)/(n1^2-n1)
	c2 = (sum(m2^2)*(n2^2) - n2)/(n2^2-n2)
	
	if(echo){
		print(sprintf("Overall cor, group 1: %2f", c1))
		print(sprintf("Overall cor, group 2: %2f", c2))
	}
	
	# Residualize full matrix, if worried about systematic correlation within groups
	if(resid.full){
		M1 = resid_DCM(M1, QN = QN)
		M2 = resid_DCM(M2, QN = QN)
	}else if(abs(c1 - c2) > .1){
		if(echo){print("Warning: Large overall correlation difference.  Consider residualizing full matrix first.")}
	}
	
	if(echo){print("Done preparing data.")}

	return(list(M1 = M1, M2 = M2, idcs = idcs))
	
}