fisher <-
function(r){
	# Set r=0 if outside range
	r = r*(abs(r) < 1)
	
	# Return fisher transform
	return(.5*log((1+r)/(1-r)))
}


# Version that truncates to 1 or -1, instead of zeroing values outside range
fisher_trunc <- function(r){
	
	# Return large number (approx infinity) for r outside range
	if(r <= -1){
		return(-10000000)
	}else if(r >= 1){
		return(10000000)
	}else{
		# Otherwise, fisher transform as usual
		return(.5*log((1+r)/(1-r)))
	}
}

# Vectorize function
fisher_trunc = Vectorize(fisher_trunc)