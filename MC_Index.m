clc
close all
clear all

%% Import
% MCs are annually starting from 1979 (year end)
% importMCactive = xlsread("Data\Raw\MC_active.xlsx");
% importMCdeadRCI = xlsread("Data\Raw\MC_deadRCIY.xlsx");
% importMCdeadnoRCI = xlsread("Data\Raw\MC_deadRCIN.xlsx");

% Prices are monthly starting from 01-01-1976
% importPricesactive = xlsread("Data\Raw\Prices_Active.xlsx");
% importPricesdeadRCI = xlsread("Data\Raw\Prices_deadRCIY.xlsx");
% importPricesdeadnoRCI = xlsread("Data\Raw\Prices_deadRCIN.xlsx");

% save("Data\importedMCData.mat","importMCactive","importMCdeadRCI","importMCdeadnoRCI");
% save("Data\importedPricesData.mat","importPricesactive","importPricesdeadRCI","importPricesdeadnoRCI");

load('Data/importedMCData.mat');
load("Data/importedPricesData.mat");


%% Matrices

MCData = [importMCactive(:,2:end),importMCdeadRCI(:,2:end),importMCdeadnoRCI(:,2:end)];

PricesData = [importPricesactive,importPricesdeadRCI,importPricesdeadnoRCI];
save("Data\PricesData.mat","PricesData");
numStocks = 100;

%% Returns (beginning of next month / beginning of previous month)

stock_returns = zeros(size(PricesData,1)-1,size(PricesData,2));

for i = 1:size(PricesData,2)
    for j = 1:size(PricesData,1)-1
    stock_returns(j,i) = PricesData(j+1,i)/PricesData(j,i) - 1;
    end
end

% Set NaN and Inf to 0
stock_returns(isnan(stock_returns))=0;
stock_returns(isinf(stock_returns))=0;
PricesData(isnan(PricesData))=0;

%% Cleaning
MCData(isnan(MCData))=0;

% set dead companies to 0 MC
for i = size(importMCactive,2)+1:size(MCData,2)
    for j = 1:size(MCData,1)-1
        if MCData(j+1,i) == MCData(j,i) && MCData(j,i) ~= 0
            MCData(j:end,i) = 0;
            break
        end
    end
end
save("Data\MCData.mat","MCData");

%% numStocks largest mature companies by MC

% Mature companies have at least 1 Year of Price Data, otherwise MC = 0
MCData_mat = MCData;
for year = 1:size(MCData_mat,1)
    for stock = 1:size(MCData_mat,2)
        if sum(PricesData(year*12+25:year*12+36,stock) == 0) >= 1
            MCData_mat(year,stock) = 0;
        end
    end
end

% LargComps(1) = Largest Companies at end of 1979
% LargComps(40) = Largest Companies at end of 2018
LargComps = zeros(size(MCData,1),numStocks);
for i = 1:size(MCData,1)
[bigcapvalue,bigindex] = maxk(MCData_mat(i,:),numStocks);
LargComps(i,:) = bigindex;
end

save("Data\LargComps.mat","LargComps");

%% Weights
% weights(1) = weights calculated at end of 1979, used for index year 1980
weights = zeros(size(MCData_mat));

for i = 1:size(MCData,1)
totalcap = sum(MCData(i,LargComps(i,:)));
weights(i,LargComps(i,:)) = MCData(i,LargComps(i,:))./totalcap;
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
% Most Data starts at end of 1980, Index starts at 01-01-1981

sum_returns = sum(weightedreturns(13:end,:),2);    %weighted returns starting at 1981

MCInvest = zeros(size(sum_returns,1),1);
MCInvest(1) = 1;
for i = 1:size(sum_returns,1)
    MCInvest(i+1) = MCInvest(i) * (1+sum_returns(i));
end

%% Performance and Risk measures
% importriskfree = xlsread("Data/Raw/RiskFreeGer.xlsx");
% yearly_riskfree = importriskfree(:,2)/100;
% save("Data/riskfree.mat","yearly_riskfree");
load("Data/riskfree.mat");
numyears = size(sum_returns,1)/12;
total_stocks = size(weights,2);
stock_returns_1981 = stock_returns_1980(13:end,:);

% Annual index returns
ann_index_returns = zeros(numyears,1);
for year = 1:numyears
      ann_index_returns(year) = prod(1+sum_returns(year*12-11:year*12))-1;
end

% Average Geometric Return
avr_ann_return = ((MCInvest(end))^(1/numyears) -1);

% Volatility
ann_vol = std(sum_returns)*sqrt(12);

% Sharpe Ratio
rfexcess_returns = ann_index_returns - yearly_riskfree;
avr_rfexcess_return = mean(rfexcess_returns);
sharpe = avr_rfexcess_return / ann_vol;

% Maximum Drawdown
maxdraw = maxdrawdown(MCInvest);

% Number of positive months
count = length(find(sum_returns>0));
pos_months = count/size(sum_returns,1)*100;

% Annual returns for all stocks (starting 1981 to 2018)
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
MCaveragesize2018 = mean(weighted_sizes);

% Results Vector
MCresults = [avr_ann_return*100 ann_vol*100 sharpe maxdraw*100 pos_months turnover*100];
MCexresults = [0 0 0];

%% FF Yearly Decomposition with MCIndex as market returns

% FFData_table = readtable("Data\Raw\Europe_3_Factors.CSV");
% FFData = table2array(FFData_table);
% Mkt_Rf = rfexcess_returns*100;
% SMB = FFData(:,3);
% HML = FFData(:,4);
% save("Data\FFData.mat","Mkt_Rf","SMB","HML");
% load("Data\FFData.mat");
% 
% tbl = table(Mkt_Rf,SMB,HML,rfexcess_returns*100,'VariableNames',{'Mkt_Rf','SMB','HML','Excess_Returns'});
% lm = fitlm(tbl,'Excess_Returns~Mkt_Rf+SMB+HML')

%% Ger FF Factors
% FFData_table = readtable("Data\Raw\GerFactors.xlsx");
% FFData = table2array(FFData_table);
% % Mkt_Rf = FFData(:,2);
% SMB = FFData(:,3);
% HML = FFData(:,4);
% WML = FFData(:,5);
% m_rf = FFData(:,1);
% Mkt_Rf = sum_returns(1:426)*100-m_rf;
% save("Data\FFData.mat","Mkt_Rf","SMB","HML","WML","m_rf");
load("Data/FFData.mat");
m_ex_returns = sum_returns(1:426)*100 - m_rf;
tbl = table(Mkt_Rf,SMB,HML,WML,m_ex_returns,'VariableNames',{'Mkt_Rf','SMB','HML','WML','Excess_Returns'});
lm = fitlm(tbl,'Excess_Returns~Mkt_Rf+SMB+HML+WML')

% %% European FF Factors
% FFData_table = readtable("Data\Raw\European_4Factors.xlsx");
% FFData = table2array(FFData_table);
% Mkt_Rf = FFData(:,2);
% SMB = FFData(:,3);
% HML = FFData(:,4);
% WML = FFData(:,5);
% m_rf = FFData(:,6);
% save("Data\FFData.mat","Mkt_Rf","SMB","HML","WML","m_rf");
% load("Data\FFData.mat");
% m_ex_returns = sum_returns(119:end)*100 - m_rf;
% tbl = table(Mkt_Rf,SMB,HML,WML,m_ex_returns,'VariableNames',{'Mkt_Rf','SMB','HML','WML','Excess_Returns'});
% lm = fitlm(tbl,'Excess_Returns~Mkt_Rf+SMB+HML+WML')

%% Last 10 years only (01-01-2009 to 31-12-2018)

recent_returns = sum_returns(337:end);
avr_recent_return = prod(1+recent_returns)^(1/10)-1;
recent_vol = std(recent_returns)*sqrt(12);
rec_rfex = ann_index_returns(29:end) - yearly_riskfree(29:end);
avr_rec_rfex = mean(rec_rfex);
rec_sharpe = avr_rec_rfex / recent_vol;

MC_recent_results = [avr_recent_return*100 recent_vol*100 rec_sharpe];

%%  Save Data for Total.m

MC_returns = sum_returns;
MC_annreturn = avr_ann_return;
save("Data/MCInvest.mat","MCInvest","MC_returns","MCresults","MC_annreturn","MCexresults","MC_recent_results","MCaveragesize2018");




