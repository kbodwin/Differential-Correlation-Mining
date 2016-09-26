% Standardize rows of matrix of data

function M = stdize(M)

    % Number of cols of matrix
    n = size(M,2);
	
	% Subtract mean of each row
	M = M - repmat(mean(M,2), 1, n);
	
	% Find sd of each row
	sds = std(M, 0, 2)*sqrt((n-1));
	
	% Return stdized
	M = M./repmat(sds, 1, n);
    
end