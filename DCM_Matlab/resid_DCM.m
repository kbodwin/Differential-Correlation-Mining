function[M] = resid_DCM(M)

    lam = FFA(M, 1);
    %imagesc(lam*lam.' - w2*w2.', [-1 1]);

    slam = lam./(lam.'*lam);
    corrs = M.'*slam;
    
    add = corrs*lam.';

    M = M - add.';
    
    % Restandardize
    M = stdize(M);
    
end
    
    