function [lower, upper] = lambdamin(A, reltol)
% Obtain rigorous bounds on the left-most eigenvalue of A
% through Cholesky factorization and bisection.
%
% This is meant to be accurate, though not necessarily efficient.
%
% Nicolas Boumal, 2016
	
	n = size(A, 1);
	assert(size(A, 2) == n, 'A must be square');
	assert(all(all(A == A')), 'A must be real, symmetric');
	
	if ~exist('reltol', 'var') || isempty(reltol)
		reltol = 1e-8;
	end
	
	% Approximately minimizing x'Ax over the unit norm vectors x
	% necessarily returns an upper bound on the left-most eigenvalue of A.
	upper = mineig_manopt(A);
	
	% Move to the left of the upper bound until we have a lower bound.
    % Here and below, this is tested using Cholesky factorization.
	lower = upper;
	coeff = 1;
	while ~is_psd(A - lower*eye(n))
		lower = lower - coeff*trace(A)/n;
		coeff = 2*coeff;
	end
	
	% We may have overshot quite a bit.
	% Now we zoom in using bisection.
	while abs((upper-lower)/((upper+lower)/2)) > reltol
		
		mid = (lower+upper)/2;
		
		if is_psd(A - mid*eye(n))
			lower = mid;
		else
			upper = mid;
		end
		
	end
	

end


function flag = is_psd(A)

	[~, p] = chol(A);
	flag = (p == 0);

end
