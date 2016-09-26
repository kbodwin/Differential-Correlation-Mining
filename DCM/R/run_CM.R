run_CM <-
function(M, seed, del = c(), echo = FALSE, alpha = 0.05, max.iter = 50){
	
	print("Version with var est")
	
	# Calculate runtime
	starttime = Sys.time()

	# Find dimensions	
	n = ncol(M) 
	p = nrow(M)
	idcs = 1:p
	
	# If we are ignoring certain rows, remove them
	if(length(del) > 0){
		idcs = idcs[-del]
		
		M = M[-del,]
	}

	# Find new dimensions
	p = nrow(M)
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
		
		xA = M[A,]		
		y = M[-A,]
		
		# Find mean vector
		mean_vec = colMeans(xA)
			
		# Find norm
		n_m = sqrt(sum(mean_vec^2))
		
		# Length of A
		k = length(A)
		
		# Find test quantities for all variables
		corsm = t(cor(mean_vec, t(y)))
		
		# Test stat and variance
		obs = corsm*n_m
		sd = sqrt(makeVars(y, xA))


		# p-values for rows not in A
		test_out = pnorm(-obs, 0, sd)
		
			
		## Calculate p-values for rows in A
		
		# Adjust means to not include row
		means = -t(t(xA) - mean_vec*k)/(k-1)

		# Find new mean norms
		n_ms = sqrt(rowSums(means*means))
		
		# Find cors of rows with means
		corsm = rowMeans(stdize(means)*xA)*n

		# Make test stat and variance
					
		obss = corsm*n_ms

		sds = sqrt(sapply(1:k, function(x) makeVar(xA[x,], xA[-x,])))

		
		# Find pvals
		test_in = pnorm(-obss, 0, sds)
		
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
		mA = apply(M[A,], 2, mean)
		
		meanA = (sum(mA^2)*(k^2) - k)/(k^2-k)
	}else{
		meanA = NA
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
	return(list(found = found, mc = meanA, its = it, time = time, pvals = test, startdels = startdels))
}
