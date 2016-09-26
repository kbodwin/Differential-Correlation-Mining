bhy <-
function(pvals, alpha = 0.05){
    
    # Sorted p-vals
    sp = sort(pvals)
    
    # Save original order of p-vals
    ord = order(pvals)
    
    # Find bhy cutoff
    nums = 1:length(pvals)
    cms = cumsum(1/nums)
    
    # Find which p-vals are less than bh cutoff
    under = sp < (nums/(length(pvals)*cms)*alpha)
    
    # Return indices of significant p-vals
    if(sum(under) == 0){
        return(c())
    }else{
        cutoff = max(which(under))
        return(ord[1:cutoff])
    }
}


bh <-
function(pvals, alpha = 0.05){
    
    # Sorted p-vals
    sp = sort(pvals)
    
    # Save original order of p-vals
    ord = order(pvals)
    
    # Find bh cutoff
    nums = 1:length(pvals)
    
    # Find which p-vals are less than bh cutoff
    under = sp < (nums/(length(pvals))*alpha)
    
    # Return indices of significant p-vals
    if(sum(under) == 0){
        return(c())
    }else{
        cutoff = max(which(under))
        return(ord[1:cutoff])
    }
}