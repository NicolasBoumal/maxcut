function X = maxcut_cvx(A)

    n = size(A, 1);
    assert(size(A, 2) == n);
    assert(nnz(A-A') == 0);
    
    cvx_begin quiet
    
        variable X(n, n) symmetric;
        
        minimize(A(:)'*X(:))
        
        X(1:(n+1):end) == 1;  % diag(X) == 1
        
        X == semidefinite(n);
    
    cvx_end

end
