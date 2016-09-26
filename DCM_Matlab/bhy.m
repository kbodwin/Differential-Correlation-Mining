function[sigs] = bhy(pvals, alpha)

    if ~exist('alpha', 'var')
        alpha = 0.05;
    end
    
    
    % Sorted p-vals and order
    [sp, ord] = sort(pvals);
    
    % Find bhy cutoff
    nums = 1:length(pvals);
    cms = cumsum(nums.^(-1));
    
    % Find which p-vals are less than bh cutoff
    under = sp < (nums./(length(pvals)*cms)*alpha);
    
    % Return indices of significant p-vals
    if sum(under) == 0
        sigs = [];
    else
        cutoff = max(under.*nums);
        sigs = ord(1:cutoff);
    end
end