function [Y, problem] = maxcut_manopt(A, p, Y0)

    n = size(A, 1);
    assert(size(A, 2) == n && nnz(A-A') == 0, 'A must be symmetric.');
    
    if ~exist('p', 'var') || isempty(p)
        p = ceil(sqrt(8*n+1)/2);
    end

    manifold = obliquefactory(p, n, true);
    
    problem.M = manifold;
    

    % Products with A dominate the cost, hence we store the result.
    function store = prepare(Y, store)
        if ~isfield(store, 'AY')
            AY = A*Y;
            store.AY = AY;
            store.diagAYYt = sum(AY .* Y, 2);
        end
    end
    
    % Define the cost function to be /minimized/.
    problem.cost = @cost;
    function [f, store] = cost(Y, store)
        store = prepare(Y, store);
        f = .5*sum(store.diagAYYt);
    end

    % Define the Riemannian gradient.
    problem.grad = @grad;
    function [G, store] = grad(Y, store)
        store = prepare(Y, store);
        % G = store.AY - bsxfun(@times, Y, store.diagAYYt);
        G = store.AY - Y.*store.diagAYYt;
    end

    % If you want to, you can specify the Hessian as well.
    problem.hess = @hess;
    function [H, store] = hess(Y, Ydot, store)
        store = prepare(Y, store);
        % SYdot = A*Ydot - bsxfun(@times, Ydot, store.diagAYYt);
        SYdot = A*Ydot - Ydot.*store.diagAYYt;
        H = manifold.proj(Y, SYdot);
    end

    if ~exist('Y0', 'var') || isempty(Y0)
        Y0 = [];
    end

    % Call your favorite solver.
    opts = struct();
    opts.verbosity = 2;
    opts.maxinner = 500;
    Y = trustregions(problem, Y0, opts);

end
