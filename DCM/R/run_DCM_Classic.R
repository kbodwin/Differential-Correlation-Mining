run_DCM_Classic <-
function(M1, M2, seed, del = c(), echo = FALSE, alpha = 0.05, max.iter = 50){
	
	print("Version with classic test")
	
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
		
		cors1 = tcrossprod(y1, xA1)
		cors2 = tcrossprod(y2, xA2)
		
		# Length of A
		k = length(A)
		
		# Test stat and variance

		test_out = sapply(1:(p-k), function(x) classicTest(cors1[x,], cors2[x,], n1 ,n2))
		
			
		## Calculate p-values for rows in A
		
		corsA1 = tcrossprod(xA1)
		corsA2 = tcrossprod(xA2)
		
		# Find pvals
		test_in = sapply(1:k, function(x) classicTest(corsA1[x,-x], corsA2[x,-x], n1 ,n2))
		
		# Combine all p-values
		test = P
		test[-A] = test_out
		test[A] = test_in
		
		# Update A to include significant rows
		newA = bhy(test, alpha = alpha)

		# Check for convergence
		if(it > max.iter){
			
			# Check for max iterations reached
			if(echo){print("Reached iteration limit.")}
			difference = 0
			newA = c()
			
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
