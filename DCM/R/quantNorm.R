# Quantile normalize a vector of data

quantNorm <-
function(Mat){
	
	# Find row length
	n = dim(Mat)[2]

	# Rank each row
	ranks = t(apply(Mat, 1, rank))
	#Mat = qnorm(ranks/(n+1))
	Mat = matrix(qnorm(ranks/(n+1)), ncol = n)
	
	return(Mat)
}
