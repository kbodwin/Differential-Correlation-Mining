function [found, meanA1, meanA2, it, time, test, startdels] = run_DCM(M1, M2, seed, del, echos, max_iter, alpha)
	
    disp('Searching...');

    % Optional arguments
     
    if ~exist('del', 'var')
        del = [];
    end
    
    if ~exist('max_iter', 'var')
        max_iter = 50;
    end
    
    if ~exist('echos', 'var')
        echos = true;
    end

    if ~exist('alpha', 'var')
        alpha = 0.05;
    end
    
	
	% Calculate runtime
	tic;

	% Find dimensions	
	[p, n1] = size(M1);
	n2 = size(M2,2);
	idcs = 1:p;
	
	% If we are ignoring certain rows, remove them

    if ~isempty(del)
		idcs(del) = [];
		
		M1(del,:) = []; 
		M2(del,:) = [];
    end

	% Find new dimensions
	p = size(M1,1);
	P = 1:p;
	
	%%%%%%% Prep Function

	% Initialize
	difference = 10;
	A = P(ismember(idcs, seed));
	prevA = p+1;
	it = 0;
	osc = false;
	diags = false;
	startdels = [];
	

       sigs = 0; %get rid of this
	% Continue searching until convergence on non-degenerate set
    while difference > floor(length(A)*.05) && length(A) > 5
		
		xA1 = M1(A,:);
		xA2 = M2(A,:);
		
		y1 = M1;
        y1(A,:) = [];
		y2 = M2;
        y2(A,:) = [];
		
		% Find mean vectors
		mean1 = mean(xA1, 1);
		mean2 = mean(xA2, 1);
			
		% Find norms
		n_m1 = sum(mean1.^2);
        n_m2 = sum(mean2.^2);
		
		% Length of A
		k = length(A);
		
		% Find test quantities for all variables

        corsm1 = stdize(mean1)*y1.';
        corsm2 = stdize(mean2)*y2.';
		
		% Observations and sd

        obs = corsm1*sqrt(n_m1) - max(corsm2.*sqrt(n_m2), 0);

       %sd = sqrt(makePhi(sqrt(n_m1)*corsm1, n_m1-1/k, k)./n1 + makePhi(sqrt(n_m2)*corsm2, n_m2-1/k, k)./n2);
		sd = sqrt(makeTau_out(y1, xA1) + makeTau_out(y2, xA2));
       
       
		% p-values for rows not in A
		test_out = normcdf(-obs, 0, sd);
		%test_out = tcdf(-obs./sd, min(n1-1, n2-1));
			
		%% Calculate p-values for rows in A
		
		% Adjust means to not include row
		mean1s = -((xA1 - ones(k,1)*mean1*k)/(k-1));
		mean2s = -((xA2 - ones(k,1)*mean2*k)/(k-1));

		% Find new mean norms
		n_m1s = sum(mean1s.^2, 2);
		n_m2s = sum(mean2s.^2, 2);
		
		% Find cors of rows with means
		corsm1 = mean(stdize(mean1s)*xA1.', 1).';
		corsm2 = mean(stdize(mean2s)*xA2.', 1).';

        obss = corsm1.*sqrt(n_m1s) - corsm2.*sqrt(n_m2s);
        %sds = sqrt(makePhi(corsm1.*sqrt(n_m1s), n_m1s - 1/(k-1), k-1)./n1 + makePhi(corsm2.*sqrt(n_m2s), n_m2s - 1/(k-1), k-1)./n2);
		sds = sqrt(makeTau_in(xA1) + makeTau_in(xA2));
		
		% Find pvals
		test_in = normcdf(-obss, 0, sds.');
        %test_in = tcdf(-obss./sds, min(n1-1, n2-1));
		
		% Combine all p-values
		test = P;
		test(~ismember(P, A)) = test_out;
		test(A) = test_in.';
		
		% Update A to include significant rows
		newA = bhy(test, alpha);

        length(newA)

        if echos
            %newA
        end

		% Check for convergence
        if it > max_iter
			
			% Check for max iterations reached
            if echos
                disp('Reached iteration limit.');
            end
            
			difference = 0;
			newA = [];
			
        else
			
			% Calculate difference between current and updated sets
			difference = sum(~ismember(newA, A)) + sum(~ismember(A, newA));
			
        end
		
		overl = sum(ismember(A, newA))
		% Check to see if algorithm is oscillating between two similar sets or jumping sideways between sets

        if(sum(~ismember(newA, prevA)) + sum(~ismember(prevA, newA)) < floor(length(newA)*.05))

            if echos
                disp('Oscillating...')
            end
			
			% If already oscillated, don't keep going.
            if osc
				
				difference = 0;
				overlap = newA(ismember(newA, A));
				
                if max(length(setdiff(overlap, newA))/length(newA), length(setdiff(overlap, A))/length(A)) < .1
					newA = overlap;
                else
					newA = [];
                end
                
				startdels = [startdels, newA, prevA(~ismember(prevA, newA))];
				osc = false;
				
            else
				
				osc = true;
				newA = A(ismember(A, newA));
                %newA = unique([A newA]);
				
            end

			
            
            
        end

		% Save old A, update current A
		prevA = A;
		A = newA;
		
		%disp(prevA)
		%disp(newA)
		
		% Print progress if desired
        if echos
            disp(sprintf('Size = %i', length(A)));
        end
		
		% Count iterations
		it = it + 1;
    end %while(difference > 0 & length(A) > 10)
	
	% New length of A
	k = length(A);
	if k > 40000
        pause
    end

	% If algorithm didn't degenerate, find mean correlations
    if k > 1
		mA1 = mean(M1(A,:), 1);
        mA2 = mean(M2(A,:), 1);
		
		meanA1 = (sum(mA1.^2)*(k^2) - k)/(k^2-k);
		meanA2 = (sum(mA2.^2)*(k^2) - k)/(k^2-k);
        
    else
        
		meanA1 = NaN;
		meanA2 = NaN;
        
    end
	
	
	% Runtime
	time = toc;
	
	% Deal with diag sitch
    if diags
		found = [0, cut, idcs(A)];
    else
		found = idcs(A);
    end

    %sigs
    
end
