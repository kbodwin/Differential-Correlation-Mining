run_DCM <-
function(M1, M2, seed, del = c(), echo = FALSE, alpha = 0.05, max.iter = 50){
	
	# Calculate runtime
	starttime = Sys.time()

	# Find dimensions	
	n1 = ncol(M1) 
	n2 = ncol(M2)
	p = nrow(M1)
	idcs = 1:p
	
	# If we are ignoring certain rows, remove them
	if(length(del) > 0){
		idcs = idcs[-del]
		
		M1 = M1[-del,]
		M2 = M2[-del,]
	}

	# Find new dimensions
	p = nrow(M1)
	P = 1:p
	
	##### Prep Function

	# Initialize
	difference = 1
	A = which(idcs %in% seed)
	prevA = p+1
	it = 0
	osc = FALSE
	diags = FALSE
	startdels = c()
	
	# Continue searching until convergence on non-degenerate set
	while(difference > 0 & length(A) > 5){
		
		xA1 = M1[A,]
		xA2 = M2[A,]
		
		y1 = M1[-A,]
		y2 = M2[-A,]
		
		# Find mean vectors
		mean1 = colMeans(xA1)
		mean2 = colMeans(xA2)
			
		# Find norms
		n_m1 = sqrt(sum(mean1^2))
		n_m2 = sqrt(sum(mean2^2))
		
		# Length of A
		k = length(A)
		
		# Find test quantities for all variables
		corsm1 = t(cor(mean1, t(y1)))
		corsm2 = t(cor(mean2, t(y2)))
		
		# Test stat and variance
		obs = corsm1*n_m1 - corsm2*n_m2
		
		#sd = sqrt(makePhi(n_m1*corsm1, n_m1^2 - 1/k, k)/n1 + makePhi(n_m2*corsm2, n_m2^2 - 1/k, k)/n2)
		#sd = sqrt(makeYvars(y1, mean1, y2, mean2))
		sd = sqrt(makeVars(y1, xA1) + makeVars(y2, xA2))


		# p-values for rows not in A
		test_out = pt(-obs/sd, min(c(n1-1, n2-1)), 0)
		#test_out = pnorm(-obs, 0, sd)
		
			
		## Calculate p-values for rows in A
		
		# Adjust means to not include row
		mean1s = -t(t(xA1) - mean1*k)/(k-1)
		mean2s = -t(t(xA2) - mean2*k)/(k-1)

		# Find new mean norms
		n_m1s = sqrt(rowSums(mean1s*mean1s))
		n_m2s = sqrt(rowSums(mean2s*mean2s))
		
		# Find cors of rows with means
		corsm1 = rowMeans(stdize(mean1s)*xA1)*n1
		corsm2 = rowMeans(stdize(mean2s)*xA2)*n2

		# Make test stat and variance
					
		obss = corsm1*n_m1s - corsm2*n_m2s
		
		#sds = sqrt(makePhi(n_m1s*corsm1, n_m1s^2 - 1/(k-1), k-1)/n1 + makePhi(n_m2s*corsm2, n_m2s^2 - 1/(k-1), k-1)/n2)
		#sds = sapply(1:k, function(x) sqrt(makeYvar(xA1[x,], mean1s[x,], xA2[x,], mean2s[x,])))
		sds = sqrt(sapply(1:k, function(x) makeVar(xA1[x,], xA1[-x,]) + makeVar(xA2[x,], xA2[-x,])))

		
		# Find pvals
		test_in = pt(-obss/sds, min(c(n1-1, n2-1)), 0)
		#test_in = pnorm(-obss, 0, sds)
		
		# Combine all p-values
		test = P
		test[-A] = test_out
		test[A] = test_in
		
		# Update A to include significant rows
		newA = bhy(test, alpha = alpha)

		# Check for convergence
		if(it >= max.iter){
			
			# Check for max iterations reached
			if(echo){print("Reached iteration limit.")}
			difference = 0
			newA = newA[newA %in% A]
			
		}else{
			
			# Calculate difference between current and updated sets
			difference = sum(!(newA %in% A)) + sum(!(A %in% newA))
			
		}
		
		
		# Check to see if algorithm is oscillating between two similar sets or jumping sideways between sets

		if(sum(!(newA %in% prevA)) + sum(!(prevA %in% newA)) == 0){

			if(echo){print("Oscillating...")}
			
			# If already oscillated, don't keep going.
			if(osc){
				
				difference = 0
				overlap = newA[newA %in% A]
				
				if(min(length(setdiff(overlap, newA))/length(newA), length(setdiff(overlap, A))/length(A)) < .05){
					newA = overlap
				}else{
					newA = c()
				}
				startdels = c(startdels, newA, prevA[!(prevA %in% newA)])
				osc = FALSE
				
			}else{
				
				osc = TRUE
				newA = A[A %in% newA]
				
			}

		}
		
		# Save old A, update current A
		prevA = A
		A = newA
		
		#print(prevA)
		#print(newA)
		
		# Print progress if desired
		if(echo){print(sprintf("Size = %i", length(A)))}	
		
		# Count iterations
		it = it + 1
	} #while(difference > 0 & length(A) > 10)
	
	# New length of A
	k = length(A)	

	# If algorithm didn't degenerate, find mean correlations
	if(k > 1){
		mA1 = apply(M1[A,], 2, mean)
		mA2 = apply(M2[A,], 2, mean)
		
		meanA1 = (sum(mA1^2)*(k^2) - k)/(k^2-k)
		meanA2 = (sum(mA2^2)*(k^2) - k)/(k^2-k)
	}else{
		meanA1 = NA
		meanA2 = NA
	}
	
	# Runtime
	time = as.double(difftime(Sys.time(), starttime, units = "secs"))
	
	# Deal with diag sitch
	if(diags){
		found = c(0, cut, idcs[A])
	}else{
		found = idcs[A]
	}
		
	# Save converged set and properties
	return(list(found = found, mc1 = meanA1, mc2 = meanA2, its = it, time = time, pvals = test, startdels = startdels))
}
