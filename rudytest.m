function rudytest()
    
    % Store the results in "record": a 3-way array
    % where the idendxing goes:
    %   graph id ; performance metric ; algorithm
    record = zeros(81, 4, 5);

    for graphid = 1:1 %[1:67, 70, 72, 77, 81]

        use_mat_files = true;
        if ~use_mat_files

            data = importdata(sprintf('Gset/G%d', graphid), ' ', 1);
            header = sscanf(data.textdata{1}, '%d %d');
            n = header(1);
            m = header(2);
            I = data.data(:, 1);
            J = data.data(:, 2);
            W = data.data(:, 3);

            A = sparse([I;J], [J;I], [W;W], n, n, 2*m);

            save(sprintf('Gset/g%d.mat', graphid), 'A', 'n', 'm');

        else

            load(sprintf('Gset/g%d.mat', graphid), 'A', 'n', 'm'); %#ok<NASGU> 

        end
        
        

        fprintf('Graph id: %d,\tn: %5d,\tm: %10d\n', graphid, n, nnz(A)/2);

        solvers = {@local_maxcut_manopt, ...
                   @local_maxcut_manopt_incremental, ...
                   @local_maxcut_sdplr, ...
                   @local_maxcut_helmberg, ...
                   @local_maxcut_cvx};

        solver_names = {'Manopt         ', ...
                        'Manopt incr.   ' ...
                        'SDPLR          ', ...
                        'Helmberg et al.', ...
                        'CVX            '};

        % Laplacian of the graph
        L = spdiags(sum(A, 2), 0, n, n) - A;
        
        for k = 1 : numel(solvers)

            % Skip CVX for large graphs
            % (it's super slow, and can crash for lack of memory)
            if (k == 5 && n >= 10000)
                record(graphid, :, k) = NaN;
                fprintf('\tSolver: %s\tlambdamin(S) in [%8g, %8g]\ttime: %8g\n', ...
                        solver_names{k}, NaN, NaN, NaN);
                continue;
            end
            
            solver = solvers{k};

            % Call the solver on the graph
            [Y, time] = solver(A);

            % Normalize the rows to ensure constraint satisfaction.
            % Y = bsxfun(@times, Y, 1./sqrt(sum(Y.^2, 2)));
            % Newer versions of Matlab accept this code instead:
            Y = Y ./ sqrt(sum(Y.^2, 2));

            % Compute lambda min of S
            S = A - spdiags(sum((A*Y).*Y, 2), 0, n, n);
            %[~, D] = mineig_manopt(S, 1);
			[low, up] = lambdamin(S, 1e-2);

            fprintf('\tSolver: %s\tlambdamin(S) in [%8g, %8g]\ttime: %8g\n', ...
                    solver_names{k}, low, up, time);

            record(graphid, 1, k) = sum(sum((L*Y).*Y))/4;
            record(graphid, 2, k) = low;
            record(graphid, 3, k) = up;
            record(graphid, 4, k) = time;
            
        end
        
    
        save rudytest.mat record;
    
%     
% 
%     % Extract n cuts from Y and compute the values via the Laplacian.
%     best_cut = -inf;
%     L = spdiags(sum(A, 2), 0, n, n) - A;
% %     for repeat = 1 : n
% %         x = sign(Y*randn(size(Y, 2), 1));
% %         cut_value = x'*L*x/4;
% %         if cut_value > best_cut
% %             best_cut = cut_value;
% %         end
% %     end
% 
%     S = A - spdiags(sum((A*Y).*Y, 2), 0, n, n);
%     % eigs(S, 6, 'SA')
%     % eigs(S, 6, 'LA')
%     [V, D] = mineig_manopt(S, 1);
%     
%     fprintf('Best cut: %d,\tbound: %.12g,\tlambdamin(S) = %g\ttime: %g\n', best_cut, trace(Y'*L*Y)/4, D, time);% \tconstraints error: %g, norm(sum(Y.*Y, 2)-1));
% 
%     % [Ydot, lambda] = hessianextreme(problem, Y, 'min');
%     % lambda
        
    end

end


% Solve with Manopt
function [Y, time] = local_maxcut_manopt(A)    
    t = tic;
    Y = maxcut_manopt(A);
    time = toc(t);
end

% Solve with Manopt incremental
function [Y, time] = local_maxcut_manopt_incremental(A)
    t = tic;
    Y = maxcut_manopt_incremental(A);
    time = toc(t);
end

% Solve with SDPLR
function [Y, time] = local_maxcut_sdplr(A)    
    t = tic;
    Y = maxcut_sdplr(A);
    time = toc(t);
end

% Solve with CVX
function [Y, time] = local_maxcut_cvx(A)    
    t = tic;
    X = maxcut_cvx(A);
    time = toc(t);
    try
        Y = chol(X)';
    catch
        % emin = eigs(X, 1, 'SA');
        % Y = chol(X - emin*eye(size(X)))';
		[V, D] = eig(X);
		d = max(0, real(diag(D)));
		Y = V*diag(sqrt(d));
    end
end

% Solve with Helmberg et al.'s IPM
function [Y, time] = local_maxcut_helmberg(A)    
    n = size(A, 1);
    L = spdiags(sum(A, 2), 0, n, n) - A;
    t = tic;
    [~, X, ~] = helmberg(L, false);
    time = toc(t);
    try
        Y = chol(X)';
    catch
        % emin = eigs(X, 1, 'SA');
        % Y = chol(X - emin*eye(size(X)))';
		[V, D] = eig(X);
		d = max(0, real(diag(D)));
		Y = V*diag(sqrt(d));
    end
end
