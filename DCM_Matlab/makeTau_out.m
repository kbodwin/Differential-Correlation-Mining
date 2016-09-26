function[vars] = makeTau_out(M, M_A)

   	%stdized to ss = 1
    	
    [p n1] = size(M);
    k = size(M_A,1);
    
    vars = 1:p;
    
    for(i = 1:p)
        
        Ui = M(i, :);

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