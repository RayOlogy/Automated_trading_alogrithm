function  [mu, Q, adj_Rsquare] = LASSO(returns, factRet, lambda, K)
    
    % Use this function for the LASSO model. Note that you will not use K 
    % in this model (K is for BSS).
    
    % for each asset, we will build LASSO to find the coefficients using
    % an appropriate lambda from our lambda candidates
    p = size(factRet,2); % number of factors
    T = size(factRet,1); % number of observations
    n = size(returns,2); % number of assets
    B = ones(p+1, n); % placeholder for coefficients
    for i =1:n
        for j=1:size(lambda,2)
            % lambda search
            lambda_temp=lambda(j);
        
            r = returns(:,i);
            
            X = [ones(T,1) factRet];
            X_au = [X zeros(T,p+1)]; % alpha and auxiliary
            
            % our optimal decision is B = [B1,...,Bp+1,y1,...,yp+1]'
            
            % our objective is ||r_i - X_au B_i||_2^2 + lambda 1^T y
            % after matrix calculation, it is equivalent to:
            %               min 2 B 'X_au' X_au B - 2r_i'X_au B + auxiliary terms * B
            Q = 2 * (X_au') * X_au; % Q = 2 X_au'X_au here
            auxiliary_term = [zeros(1,p+1) lambda_temp * ones(1,p+1)];
            c = (-2 * r' * X_au + auxiliary_term)';
            
            % constraint is: y >= B_i (B_i - y <= 0)
            %                y >= -B_i (-B_i - y <= 0)
            A = [eye(p+1)  -eye(p+1);
                -eye(p+1)  -eye(p+1)];
            b = zeros(size(A,1),1);
       
            % increase the tolerance of 'quadprog'
            options = optimoptions('quadprog','TolFun',1e-9,'Display','off');

            % use quadprog to solve optimal coefficients
            result = quadprog(Q,c,A,b,[],[],[],[],[], options);
            B_i = round(result(1:9),5); % pull coefficients and set too small coefficients to 0
            
            % Once find an appropriate alpha that 2 to 5 coefficients are
            % non-zero, stop the lambda search
            if (nnz(B_i)>=2 && nnz(B_i)<=5)
                B(:, i) = B_i;
                break
            end
        end
    end
    
    V = B(2:p+1, :); % coefficients
    
    e = returns - X * B; % residuals
    D = ones(1, size(e, 2)); % empty matrix to store variance of residuals
    for i = 1:size(e, 2)
        D(1, i) = norm(e(:,i))^2/(T-nnz(V(:, i))-1);
    end
    
    D = diag(D);
    
    mu = B' * mean(X,1)';         % n x 1 vector of asset exp. returns
    Q = V' * cov(factRet) * V + D;         % n x n asset covariance matrix
    
    % Sometimes quadprog shows a warning if the covariance matrix is not
    % perfectly symmetric.
    Q = (Q + Q')/2;
    
    % Adjusted R-square
    SSE = sum(e.^2, 1); % Sum of squared errors
    SST = sum((returns - mean(returns, 1)).^2, 1); % Total sum of squares
    adj_Rsquare = mean(1 - (SSE/(T-nnz(B(:, i))))./(SST/(T - 1)));
    %----------------------------------------------------------------------
    
end