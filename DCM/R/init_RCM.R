init_RCM <-
function(M, resp, k, start =c(), del = c()){
	# Takes pre-prepared Matrix (standardized and optionally Quantile Normalized). Looks for increasing correlation with response.
	# Response vector should be centered
    
    print("Initializer: I NEED TESTING.")
    
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

	#d = (p-k)*k # Number of in/out swap combos
	
	# Count iterations
	it = 0
	
	# Time it
	start = Sys.time()
	
	# Iterate until convergence
	while(!done){
		
		# Find mean vec of A
		mean_vec = colMeans(M[A,])
		
		# Effect on mean vector from removing elt
		#(k by n)
		out_mat = (2*t(apply(M, 1, function(x) x*mean_vec)) - M^2/k)/k

		# How to check effect of swapping row x and elts of A
		max_eff <- function(x){
			
			# Changes for elts of A leaving, elt x going in
			changes = t(apply(out_mat[A,], 1, function(y) out_mat[x,] - y)) - t(t(M[A,])*M[x,])
			
			# Actual effect is crossprod with response
			effs = tcrossprod(changes, t(resp))
			
			return(c(which.max(effs), max(effs)))
			
		}

		# Find best swap for all of them
		swaps = t(sapply(notA, function(x) max_eff(x)))
		best = which.max(swaps[,2])

		if(swaps[best, 2] <= 0){
			done = TRUE
		}else{

			# Find true variable labels
			best_out = swaps[best, 1]
			out = A[best_out]
			inn = notA[best]
		
			# Switch "out" and "in" indices from their lists
			A[best_out] = inn
			notA[best] = out
					
			# Increase iteration count
			it = it + 1	
			
			print(A)
		}
	}
	time = difftime(Sys.time(), start, units="secs")
	
	# Translate back to real indices
	A = idcs[A]
	
	return(list(seed = orig_A, found = A, iterations = it, time = time))	
}
