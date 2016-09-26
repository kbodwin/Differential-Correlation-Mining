init_CM <-
function(M, k, seed =c(), del = c()){
	# Takes pre-prepared Matrix (standardized and optionally Quantile Normalized). Looks for high correlation in M and low/negative correlation in M2, finds group of size k.
    
   if(length(del) != 0){
		# Record dimension
		real_p = dim(M)[1]
		idcs = (1:real_p)[-del]

		# Remove ignored rows
		M = M[-del,]
		
	}else{
		# Record dimension
		p = dim(M)[1]
		idcs = 1:p
	}
	
    
    # Lengths
    n = dim(M)[2]
    
    # Find size
    p = dim(M)[1]

	# Start with prespecified row set or with a random k rows.
	if(length(seed) == k){
		A = seed
	}else{
		A = sample(length(idcs), k)
	}
	
	# Save starting point
	orig_A = A

	# Make list of all indices, and of those not in seed set
	P = 1:p
	notA = P[!(P %in% A)]

	# Initialize for loop
	done = FALSE
	it = 0

	# Find correlations of all genes with only genes in A
	cross = round(M %*% t(M[A,]), digits = 10)

	# Note: resulting matrices have rows in order of data indices, columns correspond to values of A in the order listed by A

	# Rowsums represent total of pairwise correlations between [row] and A	
	rows = rowSums(cross)

	d = (p-k)*k # Number of in/out swap combos
	
	# Count iterations
	it = 0
	
	# Time it
	start = Sys.time()
	
	# Iterate until convergence
	while(!done){

			
		# Gain due to o is getting back the contribution of o from h, loss due to o is corrs for s
		# Similarly for in
		maxA <- function(i){			
			temp = rows[notA] - c(cross[notA, i])
			idx = which.max(temp)
			effect = temp[idx]
			return(c(idx, effect))
		}

		bestAs = t(sapply(1:k, function(x) maxA(x)))
		bestAs[,2] = bestAs[,2] - rows[A]
		best_out = which.max(bestAs[,2])
		best_in = bestAs[best_out,1]

		if(bestAs[best_out,2] <= 0){
			done = TRUE
		}else{
			# Find in and out data labels
			
			out = A[best_out]
			inn = notA[best_in]
		
			# Switch "out" and "in" indices from their lists
			A[best_out] = inn
			notA[best_in] = out
		
			# Find correlation vector for "in"
			new = round(M %*% M[inn,], digits = 10)
			new[inn] = 0 # So Fisher isn't inf
			new = fisher(new)*sqrt(n - 3)
			# Edit rowsums to reflect inclusion of new index, exclusion of old
			rows = rows - cross[,best_out] + new
			# Replace relevant column of corr matrix with new correlations
			cross[,best_out] = new
			
			# Increase iteration count
			it = it + 1	
		}
	}
	time = difftime(Sys.time(), start, units="secs")
	
	# Translate back to real indices
	A = idcs[A]
	
	return(list(seed = orig_A, found = A, iterations = it, time = time))	
}
