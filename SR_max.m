function  x = SR_max(mu, Q)
    
    % This function presents an implementation Sharpe Ratio Maximization
    %
    % min   y^T Q y
    % s.t.  mu^ T y = 1
    %       y >= 0

    Aeq = mu';
    beq = 1;
    lb = zeros(size(Aeq,2),1);
    options = optimoptions('quadprog','TolFun',1e-9,'Display','off');
    y = quadprog(Q,[],[],[],Aeq,beq,lb,[],[],options);
    
    % by definition, the weight x_i = y_i / sum(y)
    x = y/sum(y);
end