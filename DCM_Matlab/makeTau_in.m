function[vars] = makeTau_in(M)

   	%stdized to ss = 1
    
    	
    [k n1] = size(M);
    
    vars = 1:k;
    
    for(i = 1:k)
        
        Ui = M(i, :);
        M_A = M(1:end ~= i, :);

        r1s = M_A * Ui.';

        Ui = Ui.*sqrt(n1-1);
        W = mean(M_A,1).*sqrt(n1-1);
        Y = (r1s.' * M_A.^2).*(n1-1)./k;
        rA = mean(r1s);

        mat = 1/4.*rA^2.*Ui.^4 + rA.*Y./2.*Ui.^2 + W.^2.*Ui.^2 + Y.^2./4 - W.*Y.*Ui - rA.*W.*Ui.^3;

        var = mean(mat)./n1;
        
        vars(i) = var;
        
    end
    
end