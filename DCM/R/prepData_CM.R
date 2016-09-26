prepData_CM <-
function(M, QN = FALSE, resid.full = FALSE, echo = FALSE, strict = "low"){

	# Make sure data is numeric
	
	M = t(apply(M, 1, as.numeric))
	
	
	# Simplify
			
	M = (M - min(M, na.rm = TRUE))/max(M, na.rm = TRUE)
		
	M = round(M, 5)
	

	# Make sure 3 or more samples
	p = nrow(M)
	n = ncol(M)
	
	if(n < 3){
		print("Must have at least 3 samples.  (25+ recommended.)")
		break
	}
			
	# Remove and report bad rows
	
	killrows = sanitize_CM(M, strict)

	if(length(killrows) > 0){
		if(echo){print(sprintf("Removing %i rows - too many zero expressions.", length(killrows)))}
		
		# Record dimension
		idcs = (1:p)[-killrows]

		# Remove ignored rows
		M = M[-killrows,]
	}else{
		idcs = 1:p
	}
		
	p = length(idcs)
	
	# If not quantile normalizing, check normality and warn if non-Gaussian
	if(!QN){
		
		rand = sample((1:p), min(1000, p))
		pvs = sapply(rand, function(x) shapiro.test(M[x,])$p.value)
		
		bad = length(bh(pvs)) > 0
		
		if((bad) & echo){
			print("Warning: Data is non-normal; consider quantile normalizing.  (Set QN = TRUE)")
		}

	}
		
	print("Done checking data.")
	
	################# Process data ###################
	
	# Quantile normalize
	if(QN){
		M = quantNorm(M)
	}
	
	# Standardize (for easier calculation of correlations)
	M = stdize(M)
	
	mean_M = apply(M, 2, mean)
	
	mean_cor = (sum(mean_M^2)*(n^2) - n)/(n^2-n)
	
	if(echo){
		print(sprintf("Overall cor, group 1: %2f", mean_cor))
	}
	
	# Residualize full matrix, if worried about systematic correlation within groups
	if(resid.full){
		M = resid_CM(M, QN = QN)
	}else if(abs(mean_cor) > .1){
		if(echo){print("Warning: Large overall correlation.  Consider residualizing full matrix first.")}
	}
	
	if(echo){print("Done preparing data.")}

	return(list(M = M, idcs = idcs))
	
}