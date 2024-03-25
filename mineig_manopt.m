function [out1, out2] = mineig_manopt(A, k)
%     d  = mineig_manopt(A, k)
% [V, D] = mineig_manopt(A, k)
%
% Computes leftmost eigenspace of dimension k of A.
%
% Single output: returns the k leftmost eigenvalues of A
% Two outputs: returns V (orthonormal columns) and D (diagonal) such that V
% spans a leftmost eigenspace of A of dimension k and AV = VD.
%
% Nicolas Boumal, April 15, 2016.

    n = size(A, 1);
    assert(size(A, 2) == n);
    assert(nnz(A-A') == 0);
    
    if ~exist('k', 'var') || isempty(k)
        k = 1;
    end

    % TODO: get rid of redundant AV computation in cost and egrad.
    problem.M = grassmannfactory(n, k);
    problem.cost = @(V) .5*sum(sum(V .* (A*V)));
    problem.egrad = @(V) A*V;
    problem.ehess = @(V, Vdot) A*Vdot;
    
    V0 = [];
    opts.verbosity = 0;
    opts.tolgradnorm = 1e-8;
    VV = trustregions(problem, V0, opts);
    
    [Q, D] = eig(VV'*A*VV);
    
    V = VV*Q;
    
    if nargout <= 1
        out1 = diag(D);
    elseif nargout == 2
        out1 = V;
        out2 = D;
    else
        error('Either 1 or 2 outputs.');
    end
    
end
