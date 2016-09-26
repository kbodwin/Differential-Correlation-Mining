classicTest <- function(row_x, row_y, n1, n2){
	
	row_x = fisher(row_x)
	row_y = fisher(row_y)
	
	p = length(row_x)
	
	d = (n1 - 3)*(n2 - 3)*(n1 + n2 -6)^(-1)
	
	t3 = sqrt(d)*max(row_x - row_y)
	
	pval = (pnorm(t3))^p
	
	return(1- pval)
	
}