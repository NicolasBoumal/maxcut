function [phi, X, y] = helmberg(L, output)
%
% This code is copy-pasted from the following paper:
%
% C. Helmberg, F. Rendl, R.J. Vanderbei, and H. Wolkowicz.
% An interior-point method for semidefinite programming.
% SIAM Journal on Optimization, 6(2):342–361, 1996.
% doi:10.1137/0806020.
%
% solves: max trace(LX) s.t. X psd, diag(X) = b; b = ones(n,1)
% min b'y s.t. Diag(y) - L psd, y unconstrained,
% input: L ... symmetric matrix
% output: phi ... optimal value of primal, phi = trace(LX)
% X ... optimal primal matrix
% y ... optimal dual vector
% call: [phi, X, y] = helmberg(L);
%
if ~exist('output', 'var') || isempty(output)
    output = true;
end
digits = 6; % # significant digits of phi
[n, ~] = size( L); % problem size
b = ones( n,1 ); % any b>0 works just as well
X = diag( b); % initial primal matrix is pos. def.
y = sum( abs( L))' * 1.1; % initial y is chosen so that
Z = diag( y) - L; % initial dual slack Z is pos. def.
phi = b'*y; % initial dual
psi = L(:)' * X( :); % and primal costs
mu = Z( :)' * X( :)/( 2*n); % initial complementarity
iter=0; % iteration count
if output
    disp(' iter alphap alphad gap lower upper');
end
while phi-psi > max([1,abs(phi)]) * 10^(-digits)
    iter = iter + 1; % start a new iteration
    Zi = inv( Z); % inv(Z) is needed explicitly
    Zi = (Zi + Zi')/2;
    dy = (Zi.*X) \ (mu * diag(Zi) - b); % solve for dy
    dX = - Zi * diag( dy) * X + mu * Zi - X; % back substitute for dX
    dX = ( dX + dX')/2; % symmetrize
    % line search on primal
    alphap = 1; % initial steplength
    [~,posdef] = chol( X + alphap * dX ); % test if pos.def
    while posdef > 0
        alphap = alphap * .8;
        [~,posdef] = chol( X + alphap * dX );
    end
    if alphap < 1, alphap = alphap * .95; end % stay away from boundary
    % line search on dual; dZ is handled implicitly: dZ = diag( dy);
    alphad = 1;
    [~,posdef] = chol( Z + alphad * diag(dy) );
    while posdef > 0
        alphad = alphad * .8;
        [~,posdef] = chol( Z + alphad * diag(dy) );
    end
    if alphad < 1, alphad = alphad * .95; end
    % update
    X = X + alphap * dX;
    y = y + alphad * dy;
    Z = Z + alphad * diag(dy);
    mu = X( :)' * Z( :) / (2*n);
    if alphap + alphad > 1.8, mu = mu/2; end % speed up for long steps
    phi = b' * y; psi = L( :)' * X( :);
    % display current iteration
    if output
        disp([ iter alphap alphad (phi-psi) psi phi ]);
    end
end % end of main loop
