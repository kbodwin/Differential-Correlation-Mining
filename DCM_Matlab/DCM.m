function[DC_sets, iterations, time, meanCor1, meanCor2, indices] = DCM(M1, M2, max_groups, echos, resid_full, max_iter, max_time, est_size, start, strict)

    % Options: max_groups, max_iter, max_time, est_size, start, resid_full, echos, strict
    % defaults = [5, 50, 1, 100, [], false, true, 'low']
    
    % optional args
    if ~exist('max_groups', 'var')
        max_groups = 5;
    end
    
    if ~exist('max_iter', 'var')
        max_iter = 100;
    end
    
    if ~exist('max_time', 'var')
        max_time = 160;
    end
    
    if ~exist('est_size', 'var')
        est_size = 100;
    end
    
    if ~exist('start', 'var')
        start = [];
    end
    
    if ~exist('resid_full', 'var')
       resid_full = false;
    end
    
    if ~exist('echos', 'var')
        echos = true;
    end
    
    if ~exist('strict', 'var')
        strict = 'low';
    end
    
    % In case we find nothing, initialize
    DC_sets = {};
    iterations = [];
    time = [];
    meanCor1 = [];
    meanCor2 = [];
    indices = [];

    alpha = 0.05;

    [M1, M2, idcs] = prepData_DCM(M1, M2, resid_full, echos, strict);
    
    tic;
    
    ans = 1;
    kval = est_size;
    startdels = [];
    p = size(M1, 1);
    
    while length(startdels) < p-kval & ans <= max_groups & toc < max_time
        
        % Initialize
        if length(start) < 5
            tmp = init_DCM(M1, M2, kval, startdels);
        else
            tmp = start;
            start = [];
        end
        
        if echos
            disp('Initialized');
            %disp(tmp);
        end
        
        diag = false;
        
        % Run search procedure
        [m_found, m_meanA1, m_meanA2, m_it, m_time, m_test, m_startdels] = run_DCM(M1, M2, tmp, startdels, echos, max_iter, alpha);

        startdels = [startdels, m_startdels];
       
		res = m_found;
		
		k = length(res);
		
		% Ignore groups that are too small to be significant
		if k > 10
			
			% Check for diag blocks
			if res(1) == 0
				diag = true;
				cut = res(2);
				res([1,2]) = [];
            end
			
			if diag
				res1 = res(1:cut);
				res2 = res((cut+1):end);
				
				% Find average correlations of result within each group
				crossCor1 = mean(mean(M1(res1, :) * M1(res2, :).'));
                crossCor2 = mean(mean(M2(res1, :) * M2(res2, :).'));
				
				if echos
					disp(sprintf('Found an off-diagonal block group of size %i', k));
					disp(sprintf('cross cor 1 = %f', crossCor1));
					disp(sprintf('cross cor 2 = %f', crossCor2));
                end

                 %startdels = [startdels, res1, res2];
				
				
            else
			
				% Announce result
				if echos
					disp(sprintf('Found a group of size %i', k));
					disp(sprintf('cor 1 = %f', m_meanA1));
					disp(sprintf('cor 2 = %f', m_meanA2));
                end

                 %startdels = [startdels, res];

            end
			
			% Residualize remaining data for continued search
			M1(res, :) = resid_DCM(M1(res, :));
			M2(res, :) = resid_DCM(M2(res, :));
			
			% Save results
			if diag
				DC_sets{ans} = {'Block Diag', idcs(res1), idcs(res2)};
				meanCor1{ans} = crossCor1;
				meanCor2{ans} = crossCor2;
            else
				DC_sets{ans} = idcs(res);
				meanCor1{ans} = m_meanA1;
				meanCor2{ans} = m_meanA2;
            end
			iterations{ans} = m_it;
			time{ans} = m_time;

			
			% Count number of sets found
			ans = ans+1;
			
        else
			
			% If group is degenerate, do not try it again

			startdels = [startdels, tmp.'];
			
        end % if(k > 10)
		
		% Check total runtime
		tottime = toc;
		
    end % while(length(dels) < p1-kval & ans < max_groups)
	
	% disp reason for algorithm halting
	if echos
		if toc >= max_time
			disp('Timed out.');
        elseif(length(startdels) > p-kval)
			disp('Exhausted all searches.');
        elseif(ans > max_groups)
			disp(sprintf('%i groups found.', max_groups));
        end
    end
	
	% Return results and properties each result
end
       