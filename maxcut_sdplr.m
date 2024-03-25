function Y = maxcut_sdplr(C)

    n = size(C, 1);
    assert(size(C, 2) == n);
    assert(nnz(C-C') == 0);

    % Setting up the problem in SDPLR format (same as SeDuMi format).
    % The equality constraints are Ax = b, with x = X(:).
    % The cost is c'*x, with c = C(:).

    c = C(:);

    d = 1;
    m = n;
    
    % Build a mask for the constraints, to identify indices
    Z = sparse([], [], [], n, n, n*d);
    for i = 1 : m
        ii = (i-1)*d + (1:d);
        Z(ii, ii) = ones(d); %#ok<SPRIX>
    end
    indices = find(Z);
    n_constraints = length(indices);
    A = sparse(1:n_constraints, indices, ones(n_constraints, 1), ...
               n_constraints, n^2, n_constraints);

    Id = eye(d);
    b = repmat(Id(:), [m, 1]);

    % We further ask that X be symmetric, positive semidefinite.
    K.s = n;
    
    % Define options (see help sdplr)
    pars.feastol = 1e-6;
    pars.centol = 1e-1;
    pars.dir = 1; % 1 or 2
    pars.penfac = 2.0;
    pars.reduce = 0;
    pars.limit = 3600;
    pars.printlevel = 0;
    %pars.forcerank = p;
    pars.soln_factored = 1;
    
    % Solve with SDPLR. The actual computed solution is X = Y*Y'.
    [r, ~, ~] = sdplr(A, b, c, K, pars);
    Y = r{1};
    
end
