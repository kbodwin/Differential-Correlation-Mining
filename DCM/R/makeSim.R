require(MASS)

makeSim <- function(p, k, rho, n1, n2, bg = "n", rho2 = 0, varinf = FALSE, large_k = 0){
	
	Sigma = diag(p)
	
	if(bg == "pos"){
		Sigma[] = rho/3
	}else if(bg == "noisy"){
		temp = array(rnorm(10*p), c(p, 10))
		Sigma = cor(t(temp))
	}
	
	Sigma2 = Sigma

			
	if(large_k != 0){
		Sigma[1:large_k, 1:large_k] = rho/2
		Sigma2[1:large_k, 1:large_k] = rho/2
	}else{
		Sigma2[1:k, 1:k] = rho2
	}
	
	Sigma[1:k, 1:k] = rho

	
	if(varinf){
		
		vars = runif(p, 0, 5)	
		Sigma[row(Sigma) == col(Sigma)] = vars
		Sigma2[row(Sigma2) == col(Sigma2)] = vars	
		
	}else{
			
		Sigma[row(Sigma) == col(Sigma)] = 1
		Sigma2[row(Sigma2) == col(Sigma2)] = 1
		
	}

	
	dat1 = mvrnorm(n = n1, mu = rep(0, p), Sigma = Sigma)
	dat2 = mvrnorm(n = n2, mu = rep(0, p), Sigma = Sigma2)
	
	
	# if(k > 0){
		# Sigma = diag(k)
		# Sigma[1:k, 1:k] = rho
		# Sigma[row(Sigma) == col(Sigma)] = 1
	
		# dat = t(mvrnorm(n, rep(0,k), Sigma))
		# dat = rbind(dat, array(rnorm(n1*(p-k)), c(p-k, n1)))
	# }else{
		# dat = array(rnorm(n1*p), c(p, n1))
	# }
	
	# dat1 = array(rnorm(n1*p), c(p, n1))
	# dat2 = array(rnorm(n2*p), c(p, n2))
	
	# if(bg == "pos"){
		
		# rhon = rho/3
		# vec1 = rnorm(n1)
		# vec2 = rnorm(n2)
		
		# dat1 = t(t(sqrt(1-rhon^2)*dat1) + rhon*vec1)
		# dat2 = t(t(sqrt(1-rhon^2)*dat2) + rhon*vec2)
		
	# }else if(bg == "noisy"){
		
		# rhon = runif(p, -1, 1)
		# vec1 = rnorm(n1)
		# vec2 = rnorm(n2)
		
		# dat1 = sqrt(1-rhon^2)*dat1 + rhon %*% t(vec1)
		# dat2 = sqrt(1-rhon^2)*dat2 + rhon %*% t(vec2)		
		
	# }
	
	# vec1 = rnorm(n1)
	# vec2 = rnorm(n2)
	
	# dat1[1:k,] = t(t(sqrt(1-rho^2)*dat1[1:k,]) + rho*vec1)
	# dat2[1:k,] = t(t(sqrt(1-rho2^2)*dat2[1:k,]) + rho2*vec2)
	
	
	return(list(dat1 = t(dat1), dat2 = t(dat2)))
	#return(bind(dat1, dat2))
}
