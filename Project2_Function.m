function x = Project2_Function(periodReturns, periodFactRet, x0)

    % Use this function to implement your algorithmic asset management
    % strategy. You can modify this function, but you must keep the inputs
    % and outputs consistent.
    %
    % INPUTS: periodReturns, periodFactRet, x0 (current portfolio weights)
    % OUTPUTS: x (optimal portfolio)
    %
    % An example of an MVO implementation with OLS regression is given
    % below. Please be sure to include comments in your code.
    %
    % *************** WRITE YOUR CODE HERE ***************
    %----------------------------------------------------------------------

    % Example: subset the data to consistently use the most recent 3 years
    % for parameter estimation
    returns = periodReturns(end-35:end,:);
    factRet = periodFactRet(end-35:end,:);
    
    FMList = {'OLS' 'FF' 'LASSO' 'BSS'};
    FMList = cellfun(@str2func, FMList, 'UniformOutput', false);
    NoModels = length(FMList);
    
    % search for the estimates of Q and mu with highest adjusted R^2
    lambda = 0.005:0.005:1; 
    K = 6;
    [mu, Q, adj_Rsquare] = OLS(returns, factRet, lambda, K);
    for i = 2 : NoModels
        
        [mu_new, Q_new, adj_Rsquare_new] = FMList{i}(returns, factRet, lambda, K);
        
        if adj_Rsquare_new > adj_Rsquare
            mu = mu_new;
            Q = Q_new;
        end
        
    end
    
    kappa = 5;

    x = RP(Q, kappa);
    
    %----------------------------------------------------------------------
end
