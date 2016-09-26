init_DCM <-
function(M1, M2, k, start =c(), del = c()){
	# Takes pre-prepared Matrix (standardized and optionally Quantile Normalized). Looks for high correlation in M1 and low/negative correlation in M2, finds group of size k.
    
   if(length(del) != 0){
		# Record dimension
		real_p = dim(M1)[1]
		idcs = (1:real_p)[-del]

		# Remove ignored rows
		M1 = M1[-del,]
		M2 = M2[-del,]
	}else{
		# Record dimension
		p = dim(M1)[1]
		idcs = 1:p
	}
	
    
    # Lengths
    n1 = dim(M1)[2]
    n2 = dim(M2)[2]
    
    # Find size
    p = dim(M1)[1]

	# Start with prespecified row set or with a random k rows.
	if(length(start) == k){
		A = start
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
	cross_1 = round(M1 %*% t(M1[A,]), digits = 10)
	cross_1[cross_1 == 1] = 0 # So that transformed value will be 0 and variance won't contribute to sum
	cross_1 = fisher(cross_1)*sqrt(n1 - 3) # Fisher transform
	cross_2 = round(M2 %*% t(M2[A,]), digits = 10)
	cross_2[cross_2 == 1] = 0
	cross_2 = fisher(cross_2)*sqrt(n2 - 3)

	# Note: resulting matrices have rows in order of data indices, columns correspond to values of A in the order listed by A

	# Rowsums represent total of pairwise correlations between [row] and A	
	rows_1 = rowSums(cross_1)
	rows_2 = rowSums(cross_2)

	d = (p-k)*k # Number of in/out swap combos
	
	# Count iterations
	it = 0
	
	# Time it
	start = Sys.time()
	
	# Iterate until convergence
	while(!done){
    
      	# Rowsum Differences
		diffs12 = rows_1 - rows_2
			
		# Gain due to o is getting back the contribution of o from h, loss due to o is corrs for s
		# Similarly for in
		maxA <- function(i){			
			temp = diffs12[notA] + c(cross_2[notA, i]) - c(cross_1[notA, i])
			idx = which.max(temp)
			effect = temp[idx]
			return(c(idx, effect))
		}

		bestAs = t(sapply(1:k, function(x) maxA(x)))
		bestAs[,2] = bestAs[,2] - diffs12[A]
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
			new_1 = round(M1 %*% M1[inn,], digits = 10)
			new_1[inn] = 0 # So Fisher isn't inf
			new_1 = fisher(new_1)*sqrt(n1 - 3)
			# Edit rowsums to reflect inclusion of new index, exclusion of old
			rows_1 = rows_1 - cross_1[,best_out] + new_1
			# Replace relevant column of corr matrix with new correlations
			cross_1[,best_out] = new_1
		
			# Same thing for Group 2
			new_2 = round(M2 %*% M2[inn,], digits = 10)
			new_2[inn] = 0
			new_2 = fisher(new_2)*sqrt(n2 - 3)
			rows_2 = rows_2 - cross_2[,best_out] + new_2
			cross_2[,best_out] = new_2
			
			# Increase iteration count
			it = it + 1	
		}
	}
	time = difftime(Sys.time(), start, units="secs")
	
	# Translate back to real indices
	A = idcs[A]
	
	return(list(seed = orig_A, found = A, iterations = it, time = time))	
}
