
clc
close all
clear all

%% Matrices

load("Data/PricesData.mat");
load("Data/numStocks.mat");
load("Data/LargComps.mat");
load("Data/MCData.mat");

%% Returns

stock_returns = zeros(size(PricesData,1)-1,size(PricesData,2));

for i = 1:size(PricesData,2)
    for j = 1:size(PricesData,1)-1
    stock_returns(j,i) = PricesData(j+1,i)/PricesData(j,i) - 1;
    end
end
% Set NaN and Inf to 0
stock_returns(isnan(stock_returns))=0;
stock_returns(isinf(stock_returns))=0;


%% Weights

weights = zeros(size(MCData)-1);
for year = 1:size(MCData,1)-1
    % Calculate cov matrix for monthly returns of last 5 years, starting
    % for weights of year 1981 (calculated end of 1980)
    sigma = cov(stock_returns(year*12-11:year*12+48,LargComps(year,:)));
    
    
    minvar = @(w) w' * sigma * w;
    x0=ones(numStocks,1)/numStocks;      
    A2=[];              
    b2=[];              
    Aeq=ones(1,numStocks);     
    beq=[1];             
    lb=zeros(1,numStocks);      
    ub=ones(1,numStocks);
    [min_weights,min_var]=fmincon(minvar,x0,A2,b2,Aeq,beq,lb,ub);
    weights(year,LargComps(year,:)) = min_weights;
    
end
%weights(1)  weights calculated at end of 1980, for 1981
%weights(39) = weight at end of 2018

%% Weighted Returns
% First Weights calculated end of 1979 for year 1980
stock_returns_1980 = stock_returns(49:end,:);
stock_returns_1981 = stock_returns_1980(13:end,:);

weightedreturns = zeros(size(stock_returns_1981));
for stock  = 1:size(PricesData,2)
    for year = 1:size(weights,1)-1         %no Data for 2019 yet, but weights calculated for it
    weightedreturns(12*year-11:12*year,stock) = weights(year,stock).*stock_returns_1981(12*year-11:12*year,stock);
    end
end
weightedreturns(isnan(weightedreturns))=0;


%% MC Index Investment Course (absolute)
%Most Data starts at 1980, Index starts at 01-01-1981

sum_returns = sum(weightedreturns,2);    %weighted returns starting at 1981

MinVarInvest = zeros(size(sum_returns,1),1);
MinVarInvest(1) = 1;
for i = 1:size(sum_returns,1)
    MinVarInvest(i+1) = MinVarInvest(i) * (1+sum_returns(i));
end

%% Performance and Risk measures
load("Data/riskfree.mat");
numyears = size(sum_returns,1)/12;
total_stocks = size(weights,2);


% Annual index returns
ann_index_returns = zeros(numyears,1);
for year = 1:numyears
      ann_index_returns(year) = prod(1+sum_returns(year*12-11:year*12))-1;
end

% Average geometric return
avr_ann_return = ((MinVarInvest(end))^(1/numyears) -1);

% Volatility
ann_vol = std(sum_returns)*sqrt(12);

% Sharpe Ratio
rfexcess_returns = ann_index_returns - yearly_riskfree;
avr_rfexcess_return = mean(rfexcess_returns);
sharpe = avr_rfexcess_return / ann_vol;

% Maximum Drawdown
maxdraw = maxdrawdown(MinVarInvest);

% Number of positive months
count = length(find(sum_returns>0));
pos_months = count/size(sum_returns,1)*100;

% Annual returns (starting 1981 to 2018
ann_returns = zeros(numyears,total_stocks);
for year = 1:numyears
    for stock = 1:total_stocks
        ann_returns(year,stock) = prod(1+stock_returns_1981(year*12-11:year*12,stock))-1;
    end
end

% Annual Turnover

ann_turnover = zeros(size(weights,1)-1,1);
changed_weights = zeros(1,total_stocks);
for year = 1:numyears %weights(2) = weights at end of 1980
    for stock=1:total_stocks
        changed_weights(1,stock) = weights(year,stock)*(1+ann_returns(year,stock));
    end
    norm_weights = changed_weights./sum(changed_weights);
    change = abs(weights(year+1,:)-norm_weights);
    ann_turnover(year,1) = sum(change)/2;
end
turnover=mean(ann_turnover);

% Weighted average size of index components in 2018
weights2018 = weights(end,:);
weighted_sizes = weights2018 .* MCData(end,:);
MinVaraveragesize2018 = mean(weighted_sizes);

%Excess return to Benchmark (MC_Index), Information Ratio
load("Data/MCInvest.mat");
ret_diff = avr_ann_return- MC_annreturn;
benchex_returns = sum_returns-MC_returns;
tracking_error = std(benchex_returns)*sqrt(12);
information_ratio = ret_diff/tracking_error;


% Results Vector
MinVarresults = [avr_ann_return*100 ann_vol*100 sharpe maxdraw*100 pos_months turnover*100];
MinVarexresults = [ret_diff*100 tracking_error*100 information_ratio];

%% Factor Decomposition (Fama French)
load("Data/FFData.mat");
m_ex_returns = sum_returns(1:426)*100 - m_rf;
tbl = table(Mkt_Rf,SMB,HML,WML,m_ex_returns,'VariableNames',{'Mkt_Rf','SMB','HML','WML','Excess_Returns'});
lm = fitlm(tbl,'Excess_Returns~Mkt_Rf+SMB+HML+WML')

%% Last 10 years only (01-01-2009 to 31-12-2018)

recent_returns = sum_returns(337:end);
avr_recent_return = prod(1+recent_returns)^(1/10)-1;
recent_vol = std(recent_returns)*sqrt(12);
rec_rfex = ann_index_returns(29:end) - yearly_riskfree(29:end);
avr_rec_rfex = mean(rec_rfex);
rec_sharpe = avr_rec_rfex / recent_vol;
MinVar_recent_results = [avr_recent_return*100 recent_vol*100 rec_sharpe];

%% Save Data for Total.m

save("Data/MinVarInvest.mat","MinVarInvest","MinVarresults","MinVarexresults","MinVar_recent_results","MinVaraveragesize2018");
