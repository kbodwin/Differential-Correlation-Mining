sanitize <- function(M1, M2, strict = "high"){
	p1 = nrow(M1)
	n1 = ncol(M1)
	n2 = ncol(M2)
		
	na1 = apply(M1, 1, function(x) sum(is.na(x)) > 0)	
	na2 = apply(M2, 1, function(x) sum(is.na(x)) > 0)
	
	lb1 = apply(M1, 1, function(x) sum(x == 0) > n1/10)	
	lb2 = apply(M2, 1, function(x) sum(x == 0) > n2/10)	
	
	ub1 = apply(M1, 1, function(x) sum(x == 1) > n1/10)	
	ub2 = apply(M2, 1, function(x) sum(x == 1) > n2/10)	
	
	var1 = apply(M1, 1, function(x) var(x) < 1/(100*n1))
	var2 = apply(M2, 1, function(x) var(x) < 1/(100*n2))

		
	if(strict == "high"){
		
		killrows = (1:p1)[na1 | na2 | lb1 | lb2 | ub1 | ub2 | var1 | var2]
	
	}else{
		
		killrows = (1:p1)[na1 | na2 | lb1 | lb2 | ub1 | ub2]
	
	}

	return(killrows)
	
}