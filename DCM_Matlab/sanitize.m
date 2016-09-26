function [M1, M2, idcs] = sanitize(M1, M2, strict, echos)

    if ~exist('strict', 'var')
        strict = 'low';
    end
    
    if ~exist('echos', 'var')
       echos = true;
    end

    [p, n1] = size(M1);
    n2 = size(M2, 2);
    
    na1 = [];
    na2 = [];
    lb1 = [];
    lb2 = [];
    ub1 = [];
    ub2 = [];
    var1 = [];
    var2 = [];
    
    for i = 1:p
        
        row1 = M1(i,:);
        row2 = M2(i,:);
        
        na1 = [na1, sum(isnan(row1)) > 0];
        na2 = [na2, sum(isnan(row2)) > 0];
        
        lb1 = [lb1, sum(row1 == 0) > n1/10];
        lb2 = [lb2, sum(row2 == 0) > n2/10];
        
        ub1 = [ub1, sum(row1 == 1) > n1/10];
        ub2 = [ub2, sum(row2 == 1) > n2/10];
	
        var1 = [var1, var(row1) < 1/(100*n1)];
        var2 = [var2, var(row2) < 1/(100*n2)];
	
    end
    
	if(strcmp(strict,'high'))
		
		killrows = logical(na1 + na2 + lb1 + lb2 + ub1 + ub2 + var1 + var2);
	
    else
		
		killrows = logical(na1 + na2 + lb1 + lb2 + ub1 + ub2);
	
    end

    if echos
       disp(sprintf('Removing %i rows - too many zero expressions.', sum(killrows)));
    end

	M1 = M1(~killrows,:);
    M2 = M2(~killrows,:);
    P = [1:p];
    idcs = P(~killrows);
    
end