function  [mu, Q, adj_Rsquare] = BSS(returns, factRet, lambda, K)
    
    % Use this function for the BSS model. Note that you will not use 
    % lambda in this model (lambda is for LASSO).

    p = size(factRet,2); % number of factors
    T = size(factRet,1); % number of observations
    n = size(returns,2); % number of assets
    B = ones(p+1, n); % placeholder for coefficients

    for i =1:n
        r = returns(:,i);
        
        X = [ones(T,1) factRet];
        X_au = [X zeros(T,p+1)]; % alpha and auxiliary
        
        % our optimal decision is B = [B1,...,Bp+1,y1,...,yp+1]'
            
        % our objective is ||r_i - XB_i||_2^2
        % after matrix calculation, it is equivalent to:
        %               min B'X_au' X_au B - 2r_i'X_au B
        Q = (X_au') * X_au; % Q = X_au'X_au here
        c = -2 * r' * X_au;
        
        % inequality constraint is: Ly <= B_i <= Uy, 
        % here we set L=-100, U=100
        % we have B_i - 100 y <= 0
        %         -B_i - 100 y <= 0
        A = [eye(p+1)  -100 * eye(p+1);
            -eye(p+1)  -100 * eye(p+1)];
        b = zeros(size(A,1),1);
        
        % equality constraint: 1^T y = K
        Aeq = [zeros(1,p+1) ones(1,p+1)];
        beq = K;
         
        % upper and lower bound: -100 <= B_i <= 100, 0 <= y_i <= 1
        lb = [-100 * ones(p+1,1); zeros(p+1,1)];
        ub = [100 * ones(p+1,1); ones(p+1,1)];
        
        % Define the variable types:'C' defines a continuous variable, 'B' defines
        % a binary variable
        varTypes=[repmat('C',p+1,1);repmat('B',p+1,1)];
        
        clear model;
        
        % Gurobi accepts an objective function of the following form:
        % f(x) = (1/2) x' Q x + obj' x 
        model.Q = sparse(Q);
        
        % define the c vector in the objective (which is a vector of zeros since
        % there is no linear term in our objective)
        model.obj = c;
        
        % Gurobi only accepts a single A matrix, with both inequality and equality
        % constraints
        model.A = [sparse(A); sparse(Aeq)];
        
        % Define the right-hand side vector b
        model.rhs = [b; beq];
        
        % Indicate whether the constraints are ">=", "<=", or "="
        model.sense = [ repmat('<', (2*p+2), 1); '=' ];
        
        % Define the variable type (continuous and binary)
        model.vtype = varTypes;
        
        % Define the variable upper and lower bounds
        model.lb = lb;
        model.ub = ub;
        
        % Set some Gurobi parameters to limit the runtime and to avoid printing the
        % output to the console. 
        clear params;
        params.TimeLimit = 100;
        params.OutputFlag = 0;

        results = gurobi(model,params);
        B(:, i) = results.x(1:9); % pull coefficients
    end
    
    V = B(2:p+1, :); % coefficients
    
    e = returns - X * B; % residuals
    D = ones(1, size(e, 2)); % empty matrix to store variance of residuals
    for i = 1:size(e, 2)
        D(1, i) = norm(e(:,i))^2/(T-K);
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
    adj_Rsquare = mean(1 - (SSE/(T-K))./(SST/(T - 1)));
    %----------------------------------------------------------------------
    
end