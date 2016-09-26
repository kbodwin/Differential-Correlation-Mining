# Removes single dimension of groupwise correlation in a dataset

resid_DCM <-
function(Mat, QN = FALSE){
    
    require(FAMT)
    
    # Make sure Mat is data frame save the ordering because it be reordered by famt
    Mat = data.frame(Mat)
    n = ncol(Mat)
    
    # Convert to FAMT data
    mat = as.FAMTdata(Mat)
    ord = sapply(names(Mat),function(x) which(x == names(mat$expression)))
    
    # Use EM to estimate B, Z
    fa = emfa_DCM(mat, nbf = 1)
    
    # Calculate group correlation effect
    add = tcrossprod(fa$B, fa$Factors[ord]/sqrt(n-1))
    
    # Subtract group correlation from data
    Mat = as.matrix(Mat - add)
    
    # If dataset was quantile normalize, re-normalize it
    if(QN){
        Mat = quantNorm(Mat)
    }
    
    # Re-standardize data
    Mat = stdize(Mat)
    
    # Return residualized data
    return(Mat)
}
