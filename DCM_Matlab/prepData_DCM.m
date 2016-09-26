function[M1, M2, idcs] = prepData_DCM(M1, M2, resid_full, echos, strict)

    % Optional args
    if ~exist('resid_full', 'var')
       resid_full = true;
    end
    
    if ~exist('echos', 'var')
        echos = false;
    end
    
    if ~exist('strict', 'var')
        strict = 'low';
    end

    % Simplify rows to a 0 to 1 scale
     
    [p, n1] = size(M1);
    [p2, n2] = size(M2);
    
    M1 = (M1 - repmat(min(M1, [], 2), 1, n1));
    M1 = M1./repmat(max(M1, [], 2), 1, n1);
	M2 = (M2 - repmat(min(M2, [], 2), 1, n2));
    M2 = M2./repmat(max(M2, [], 2), 1, n2);
		
	M1 = round(M1, 5);
	M2 = round(M2, 5);
    
    % Check that rows match
    if p ~= p2
        error('Datasets must have the same number of rows');
    end
    
    % Make sure 3 or more samples
    
    if n1 < 3 | n2 < 3
        error('Must have at least 3 samples in each group');
    end
    
    % Remove bad rows (NaN, too many at upper or lower bound, low variance)
    
    [M1, M2, idcs] = sanitize(M1, M2, strict, echos);
    
    p = size(idcs, 2);
    
    if echos
        disp('Done checking data.');
    end
    
    
    % Process data
    
    % Standardize
    M1 = stdize(M1);
    M2 = stdize(M2);
    
    
    % find cors
    m1 = mean(M1, 1);
    m2 = mean(M2, 1);
    
    c1 = (sum(m1.^2)*(n1^2) - n1)/(n1^2-n1);
	c2 = (sum(m2.^2)*(n2^2) - n2)/(n2^2-n2);
    	
	if echos
		disp(sprintf('Overall cor, group 1: %2f', c1));
		disp(sprintf('Overall cor, group 2: %2f', c2));
    end
   
    % Residualize full matrix, or check for large overall difference.
	
    if resid_full
        
        if echos
             disp('Residualizing...');
        end
        
        M1 = resid_DCM(M1);
        M2 = resid_DCM(M2);
        
        m1 = mean(M1, 1);
        m2 = mean(M2, 1);
    
        c1 = (sum(m1.^2)*(n1^2) - n1)/(n1^2-n1);
        c2 = (sum(m2.^2)*(n2^2) - n2)/(n2^2-n2);
        
        if echos
            disp(sprintf('Overall cor, group 1: %2f', c1));
            disp(sprintf('Overall cor, group 2: %2f', c2));
        end
     
    elseif abs(c1 - c2) > 0.1
            
        if(echos)
            disp('Warning: Large overall correlation difference.  Consider residualizing full matrix first.')   
        end
        
    end
    
            
   if echos
       disp('Done preparing data.');
   end
    
end
    
    
    
    
        