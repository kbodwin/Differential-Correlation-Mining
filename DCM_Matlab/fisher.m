function[fr] = fisher(r)
	% Set r=0 if outside range
	r = r.*(abs(r) < 1) + .999999*sign(r).*(abs(r) >= 1);
	
	% Return fisher transform
	fr = .5*log((1+r)./(1-r));
end