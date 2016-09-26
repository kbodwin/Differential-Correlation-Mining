sanitize_CM <- function(M, strict = "high"){
	p = nrow(M)
	n = ncol(M)
		
	na = apply(M, 1, function(x) sum(is.na(x)) > 0)	
	
	lb = apply(M, 1, function(x) sum(x == 0) > n/2)	
	
	ub = apply(M, 1, function(x) sum(x == 1) > n/2)	
	
	vars = apply(M, 1, function(x) var(x) < 1/(100*n))

		
	if(strict == "high"){
		
		killrows = (1:p)[na | lb | ub | vars]
	
	}else{
		
		killrows = (1:p)[na | lb | ub ]
	
	}

	return(killrows)
	
}