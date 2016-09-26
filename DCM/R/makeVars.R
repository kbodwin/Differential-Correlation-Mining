makeVars <- function(Uis, M_A){
	
	# stdized to ss = 1
	
	n1 = ncol(Uis)
	
	r1s = M_A %*% t(Uis)
	k = length(r1s)
	
	Uis = Uis*sqrt(n1-1)
	W = colMeans(M_A)*sqrt(n1-1)
	Y = (t(r1s) %*% M_A^2)*(n1-1)/k
	rA = mean(r1s)
	
	mat = 1/4*rA^2*Uis^4 + rA*Y/2*Uis^2 + W^2*Uis^2 + Y^2/4 - W*Y*Uis - rA*W*Uis^3
	
	return(rowMeans(mat)/n1)
	
}

makeVar <- function(Ui, M_A){
	
	# stdized to ss = 1
	
	n1 = length(Ui)
	k = dim(M_A)[1]
	
	r1s = M_A %*% Ui
	
	Ui = Ui*sqrt(n1-1)
	W = colMeans(M_A)*sqrt(n1-1)
	Y = (t(r1s) %*% M_A^2)*(n1-1)/k
	rA = mean(r1s)
	
	mat = 1/4*rA^2*Ui^4 + rA*Y/2*Ui^2 + W^2*Ui^2 + Y^2/4 - W*Y*Ui - rA*W*Ui^3

	
	return(mean(mat)/n1)
	
}