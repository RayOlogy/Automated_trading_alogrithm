% This script is mainly dealing with the portfolio performance for each of
% our optimization alogrithm.

clc
clear all
format short

% Program Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Read input files 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input file names
assetData  = 'MMF1921_AssetPrices_3.csv';
factorData = 'MMF1921_FactorReturns_3.csv';

% Initial budget to invest ($100,000)
initialVal = 100000;

% Length of investment period (in months)
investPeriod = 6;

% Load the stock weekly prices
adjClose = readtable(assetData);
adjClose.Properties.RowNames = cellstr(datetime(adjClose.Date));
adjClose.Properties.RowNames = cellstr(datetime(adjClose.Properties.RowNames));
adjClose.Date = [];

% Load the factors weekly returns
factorRet = readtable(factorData);
factorRet.Properties.RowNames = cellstr(datetime(factorRet.Date));
factorRet.Properties.RowNames = cellstr(datetime(factorRet.Properties.RowNames));
factorRet.Date = [];

riskFree = factorRet(:,9);
factorRet = factorRet(:,1:8);

% Identify the tickers and the dates 
tickers = adjClose.Properties.VariableNames';
dates   = datetime(factorRet.Properties.RowNames);

% Calculate the stocks' weekly EXCESS returns
prices  = table2array(adjClose);
returns = ( prices(2:end,:) - prices(1:end-1,:) ) ./ prices(1:end-1,:);
returns = returns - ( diag( table2array(riskFree) ) * ones( size(returns) ) );
returns = array2table(returns);
returns.Properties.VariableNames = tickers;
returns.Properties.RowNames = cellstr(datetime(factorRet.Properties.RowNames));

% Align the price table to the asset and factor returns tables by
% discarding the first observation.
adjClose = adjClose(2:end,:);

% Start of out-of-sample test period 
testStart = datetime(returns.Properties.RowNames{1}) + calyears(5);

% End of the first investment period
testEnd = testStart + calmonths(investPeriod) - days(1);

% End of calibration period (note that the start date is the first
% observation in the dataset)
calEnd = testStart - days(1);

% Total number of investment periods
NoPeriods = ceil( days(datetime(returns.Properties.RowNames{end}) - testStart) / (30.44*investPeriod) );

% Number of assets      
n = size(adjClose,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Run your program
% 
% This section will run your Project1_Function in a loop. The data will be
% loaded progressively as a growing window of historical observations.
% Rebalancing will take place after every loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Preallocate space for the portfolio per period value and turnover
currentVal = zeros(NoPeriods, 5);

% Initiate counter for the number of observations per investment period
toDay = 0;

for t = 1 : NoPeriods
  
    % Subset the returns and factor returns corresponding to the current
    % calibration period.
    periodReturns = table2array( returns( dates <= calEnd, :) );
    periodFactRet = table2array( factorRet( dates <= calEnd, :) );
    currentPrices = table2array( adjClose( ( calEnd - calmonths(1) - days(5) ) <= dates & dates <= calEnd, :) )';
    
    % Subset the prices corresponding to the current out-of-sample test 
    % period.
    periodPrices = table2array( adjClose( testStart <= dates & dates <= testEnd,:) );
    
    % Set the initial value of the portfolio or update the portfolio value
    if t == 1
        
        currentVal(t,:) = initialVal;
        
    else
        for i = 1 : 5
            
            currentVal(t,i) = currentPrices' * NoShares{i};
            
        end
    end
    
    %----------------------------------------------------------------------
    % Portfolio optimization
    % You must write code your own algorithmic trading function 
    %----------------------------------------------------------------------
    returns_temp = periodReturns(end-35:end,:);
    factRet_temp = periodFactRet(end-35:end,:);
    
    FMList = {'OLS' 'FF' 'LASSO' 'BSS'};
    FMList = cellfun(@str2func, FMList, 'UniformOutput', false);
    NoModels = length(FMList);
    
    % search for the estimates of Q and mu with highest adjusted R^2
    lambda = 0.005:0.005:1; 
    K = 6;
    [mu, Q, adj_Rsquare] = OLS(returns_temp, factRet_temp, lambda, K);
    R{1}(t) = adj_Rsquare; % store adjusted R squares 
    for i = 2 : NoModels
        
        [mu_new, Q_new, adj_Rsquare_new] = FMList{i}(returns_temp, factRet_temp, lambda, K);
        R{i}(t) = adj_Rsquare_new;
        if adj_Rsquare_new > adj_Rsquare
            mu = mu_new;
            Q = Q_new;
        end
        
    end
    
    % Example: Use MVO to optimize our portfolio
    alpha = 0.95;
    lambda = 0.02;
    kappa = 5;

    x{1}(:,t) = MVO(mu, Q);
    x{2}(:,t) = CVaR(mu, returns_temp, alpha);
    x{3}(:,t) = robustMVO(mu, Q, returns_temp, alpha, lambda);
    x{4}(:,t) = RP(Q, kappa);
    x{5}(:,t) = SR_max(mu, Q);
        
    % Update counter for the number of observations per investment period
    fromDay = toDay + 1;
    toDay   = toDay + size(periodPrices,1);
    
    for i = 1 : 5
        % Number of shares your portfolio holds per stock
        NoShares{i} = x{i}(:,t) .* currentVal(t, i) ./ currentPrices;

        % Weekly portfolio value during the out-of-sample window
        portfValue(fromDay:toDay, i) = periodPrices * NoShares{i};
    
    end
    % Update your calibration and out-of-sample test periods
    testStart = testStart + calmonths(investPeriod);
    testEnd   = testStart + calmonths(investPeriod) - days(1);
    calEnd    = testStart - days(1);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% 3.1 Portfolio wealth evolution plot
%--------------------------------------------------------------------------

% Calculate the dates of the out-of-sample period
plotDates = dates(dates >= datetime(returns.Properties.RowNames{1}) + calyears(5) );

fig1 = figure(1);

for i = 1 : 5
    
    plot( plotDates, portfValue(:,i) )
    hold on
    
end

tags = {'MVO' 'CVaR' 'Robust MVO' 'Risk Parity' 'Sharpe Ratio Maximize'};

legend(tags, 'Location', 'eastoutside','FontSize',12);
datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('Portfolio value', 'FontSize', 14)
ylabel('Value','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig1,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig1,'Position');
set(fig1,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig1,'fileName','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig1,'Portfolio Value','-dpng','-r0');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program End