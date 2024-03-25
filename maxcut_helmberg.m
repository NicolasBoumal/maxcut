function X = maxcut_helmberg(A)
    n = size(A, 1);
    L = spdiags(sum(A, 2), 0, n, n) - A;
    [~, X, ~] = helmberg(L, false);
end
