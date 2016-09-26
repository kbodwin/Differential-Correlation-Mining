CM <-
function(M, max.groups = 5, max.iter = 50, max.time = 1, est.size = nrow(M)/10, alpha = 0.05, start = c(), QN = FALSE, resid.full = FALSE, echo = FALSE, strict = "low"){
	
	################### Prepare Data #####################
	
	prep = prepData_CM(M, QN, resid.full, echo, strict)
	M = prep$M
	idcs = prep$idcs
	
	p = nrow(M)
	
	################### Run Algorithm ####################

	# Prepare to store results
	DC_sets = list()
	iterations = list()
	time = list()
	meanCor = list()


	# Initialize
	starttime = Sys.time()
	tottime = 0
	ans = 1
	kval = est.size
	startdels = c()
	
	# Keep searching until all rows or exhausted, or max groups or time is reached
	while(length(startdels) < p-kval & ans <= max.groups & tottime < max.time){
		
		# Find initial starting set, ignoring already tried starting points
		if(length(start) < 5){
			tmp = init_CM(M, k=kval, del = startdels)$found
		}else{
			tmp = start
			start = c()
		}
	
		if(echo){
			print("Initialized")
			print(tmp)
		}
		
		diag = FALSE
		
		# Run search procedure with specified values
		CM = run_CM(M, seed = tmp, echo = echo, max.iter = max.iter)
		
		startdels = c(startdels, CM$startdels)

		res = CM$found
		
		k = length(res)
		
		# Ignore groups that are too small to be significant
		if(k > 10){
			
			# Check for diag blocks
			if(res[1] == 0){
				diag = TRUE
				cut = res[2]
				res = res[-c(1,2)]
			}
			
			if(diag){
				res1 = res[1:cut]
				res2 = res[-(1:cut)]
				
				# Find average correlations of result within each group
				#crossCor = mean(cor(t(M[res1,]), t(M[res2,])))
				
				if(echo){
					print(sprintf("Found an off-diagonal block group of size %i", k))
					#print(sprintf("cross cor 1 = %f", crossCor1))
					#print(sprintf("cross cor 2 = %f", crossCor2))
				}
				
				
			}else{
			
				# Announce result
				if(echo){
					print(sprintf("Found a group of size %i", k))
					print(sprintf("cor = %f", CM$mc))
				}
			
			}
			
			# Residualize remaining data for continued search
			M[res,] = resid_CM(M[res,], QN = QN)
			
			# Save results
			if(diag){
				DC_sets[[ans]] = list("Block Diag", idcs[res1], idcs[res2])
				#meanCor1[[ans]] = crossCor1
				#meanCor2[[ans]] = crossCor2				
				meanCor1[[ans]] = 0
				meanCor2[[ans]] = 0
			}else{
				DC_sets[[ans]] = idcs[res]
				meanCor[[ans]] = CM$mc
			}
			iterations[[ans]] = CM$its
			time[[ans]] = CM$time

			
			# Count number of sets found
			ans = ans+1
			
		}else{
			
			# If group is degenerate, do not try it again
			startdels = c(startdels, tmp)
			
		}# if(k > 10)
		
		# Check total runtime
		tottime = difftime(Sys.time(), starttime, units = "hours")
		
	} # while(length(dels) < p-kval & ans < max.groups)
	
	# Print reason for algorithm halting
	if(echo){
		if(tottime >= max.time){
			print("Timed out.")
		}else if(length(startdels) > p-kval){
			print("Exhausted all searches.")
		}else if(ans > max.groups){
			print(sprintf("%i groups found.", max.groups))
		}
	}
	
	# Return results and properties each result
	return(list(DC_sets = DC_sets, iterations = iterations, time = time, meanCor = meanCor, indices = idcs))

}
