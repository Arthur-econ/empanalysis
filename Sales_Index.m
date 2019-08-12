clc
close all
clear all

%% Import

% importSalesactive = xlsread("Data/Raw/Sales_active.xlsx");
% importSalesdeadRCI = xlsread("Data/Raw/Sales_deadRCIY.xlsx");
% importSalesdeadRCI = [importSalesdeadRCI,zeros(size(importSalesdeadRCI,1),7)];
% importSalesdeadnoRCI= xlsread("Data/Raw/Sales_deadRCIN.xlsx");
% save("Data/importedSalesData.mat","importSalesactive","importSalesdeadRCI","importSalesdeadnoRCI"); 

load("Data/importedSalesData.mat");

%% Matrices

load("Data/PricesData.mat");
SalesData = [importSalesactive(:,2:end),importSalesdeadRCI(:,2:end),importSalesdeadnoRCI(:,2:end-1)];
save("Data/SalesData.mat","SalesData");
load("Data/LargComps.mat");
load("Data/numStocks.mat");
load("Data/MCData.mat");


%% Returns

stock_returns = zeros(size(PricesData,1)-1,size(PricesData,2));

for i = 1:size(PricesData,2)
    for j = 1:size(PricesData,1)-1
    stock_returns(j,i) = PricesData(j+1,i)/PricesData(j,i) - 1;
    end
end
% Set NaN to 0
stock_returns(isnan(stock_returns))=0;
stock_returns(isinf(stock_returns))=0;

%% Cleaning
SalesData(isnan(SalesData))=0;

%% Weights
weights = zeros(size(SalesData));

for i = 1:size(SalesData,1)
total = sum(SalesData(i,LargComps(i,:)));
weights(i,LargComps(i,:)) = SalesData(i,LargComps(i,:))./total;
end

%% Weighted Returns
% First Weights calculated end of 1979 for year 1980
stock_returns_1980 = stock_returns(49:end,:);

weightedreturns = zeros(size(stock_returns_1980));
for stock  = 1:size(PricesData,2)
    for year = 1:size(weights,1)-1         %no Data for 2019 yet, but weights calculated for it
    weightedreturns(12*year-11:12*year,stock) = weights(year,stock).*stock_returns_1980(12*year-11:12*year,stock);
    end
end
weightedreturns(isnan(weightedreturns))=0;


%% MC Index Investment Course (absolute)
%Most Data starts at 1980, Index starts at 01-01-1981

sum_returns = sum(weightedreturns(13:end,:),2);    %weighted returns starting at 1981

SalesInvest = zeros(size(sum_returns,1),1);
SalesInvest(1) = 1;
for i = 1:size(sum_returns,1)
    SalesInvest(i+1) = SalesInvest(i) * (1+sum_returns(i));
end

%% Performance and Risk measures
load("Data/riskfree.mat");
numyears = size(sum_returns,1)/12;
total_stocks = size(weights,2);
stock_returns_1981 = stock_returns_1980(13:end,:);

% Annual index returns
ann_index_returns = zeros(numyears,1);
for year = 1:numyears
      ann_index_returns(year) = prod(1+sum_returns(year*12-11:year*12))-1;
end

% Average geometric return
avr_ann_return = ((SalesInvest(end))^(1/numyears) -1);

% Volatility
ann_vol = std(sum_returns)*sqrt(12);

% Sharpe Ratio
rfexcess_returns = ann_index_returns - yearly_riskfree;
avr_rfexcess_return = mean(rfexcess_returns);
sharpe = avr_rfexcess_return / ann_vol;

% Maximum Drawdown
maxdraw = maxdrawdown(SalesInvest);

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

ann_turnover = zeros(size(weights,1)-2,1);
changed_weights = zeros(1,total_stocks);
for year = 1:numyears %weights(2) = weights at end of 1980
    for stock=1:total_stocks
        changed_weights(1,stock) = weights(year+1,stock)*(1+ann_returns(year,stock));
    end
    norm_weights = changed_weights./sum(changed_weights);
    change = abs(weights(year+2,:)-norm_weights);
    ann_turnover(year,1) = sum(change)/2;
end
turnover=mean(ann_turnover);

% Weighted average size of index components in 2018
weights2018 = weights(40,:);
weighted_sizes = weights2018 .* MCData(40,:);
Salesaveragesize2018 = mean(weighted_sizes);


%Excess return to Benchmark (MC_Index), Information Ratio
load("Data/MCInvest.mat");
ret_diff = avr_ann_return- MC_annreturn;
benchex_returns = sum_returns-MC_returns;
tracking_error = std(benchex_returns)*sqrt(12);
information_ratio = ret_diff/tracking_error;


% Results Vector
Salesresults = [avr_ann_return*100 ann_vol*100 sharpe maxdraw*100 pos_months turnover*100];
Salesexresults = [ret_diff*100 tracking_error*100 information_ratio];

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

Sales_recent_results = [avr_recent_return*100 recent_vol*100 rec_sharpe];

%% Save Data for Total.m

save("Data/SalesInvest.mat","SalesInvest","Salesresults","Salesexresults","Sales_recent_results","Salesaveragesize2018");