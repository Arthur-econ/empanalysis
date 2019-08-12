clc
clear all
close all

%% Run Files
% run("MC_Index.m");
% run("BV_Index.m");
% run("Employee_Index.m");
% run("Equal_Index.m");
% run("MinVar_Index.m");
% run("Random_Index.m");
% run("Sales_Index.m");
% run("Div_Index.m");

clc


%% Load Data

load("Data/MCInvest.mat");
load("Data/DivInvest.mat");
load("Data/BVInvest.mat");
load("Data/EmpInvest.mat");
load("Data/SalesInvest.mat");
load("Data/EqInvest.mat");
load("Data/MinVarInvest.mat");
load("Data/RandInvest.mat");
load("Data/numStocks.mat");

%% Matrices

Investments = [MCInvest, DivInvest, BVInvest, EmpInvest, SalesInvest, EqInvest, RandInvest , MinVarInvest];

%% Comparison Table
% Total results
results = [MCresults; BVresults; Divresults; Empresults; Salesresults; Eqresults; Randresults; MinVarresults;];
Results_table = table(results(:,1),results(:,2),results(:,3),results(:,4),results(:,5),...
                results(:,6),'VariableNames',{'Average_Annual_Return','Annual_Volatility','Sharpe_Ratio',...
                'Maximum_Drawdown','Positive_Months','Turnover'});
Results_table.Properties.RowNames = {'Market Cap','Dividend','Book Value','Employees','Sales','Equal',...
                                     'Random', 'Min Var'}
                                 
% Benchmarked with MC Index                          
Benchmark_results = [MCexresults; BVexresults; Divexresults; Empexresults; Salesexresults;...
                    Eqexresults;  Randexresults; MinVarexresults;];
Benchmark_table = table(Benchmark_results(:,1),Benchmark_results(:,2),Benchmark_results(:,3),'VariableNames',...
                    {'Excess_Return_over_Benchmark','Tracking_Error','Information_Ratio'});
Benchmark_table.Properties.RowNames = {'Market Cap','Dividend','Book Value','Employees','Sales','Equal',...
                                     'Random', 'Min Variance'}
                                 
% Recent comparison (last 10 years)
Recent_results = [MC_recent_results; BV_recent_results; Div_recent_results; Emp_recent_results; Sales_recent_results;...
                    Eq_recent_results;  Rand_recent_results; MinVar_recent_results;];
Recent_table = table(Recent_results(:,1),Recent_results(:,2),Recent_results(:,3),'VariableNames',...
                    {'Average_Annual_Return','Volatility','Sharpe_Ratio'});
Recent_table.Properties.RowNames = {'Market Cap','Dividend','Book Value','Employees','Sales','Equal',...
                                     'Random','Min Variance'}

%% Plots

t = datetime(1981,1,1)+calmonths(0:size(MCInvest,1)-1);

figure(1)
hold on
plot(t,Investments(:,1),'--')
plot(t,Investments(:,2:end))
xlabel("Date");
ylabel("Value in €");
legend("Market Cap Index","Dividend Index","Book Value Index",...
    "Employee Index","Sales Index",...
    "Equal Index", "Random Index", "MinVar Index",...
    "Location","northwest");
title("Growth of a €1 Investment");
%title("Absolute Performance of Alternative Indices ("+numStocks+" Largest German Companies)");
set(gca,'FontSize',20)
hold off

%%

figure(2)

semilogy(t,Investments(:,1),'--',t,Investments(:,2:end))

xlabel("Date");
ylabel("Logarithmic Scale");
legend("Market Cap Index","Dividend Index","Book Value Index",...
    "Employee Index","Sales Index",...
    "Equal Index", "Random Index", "MinVar Index",...
    "Location","northwest");
title("Logarithmic Scale of Index Investments");
  set(gca,'ytick',[1 2 4 8 16])
  set(gca,'FontSize',20)
ylim([0 40])

%% Turnover Bar Graph

turnover = results(:,end);
x = categorical({'Market Cap','Dividend','Book Value',...
    'Employee','Sales','Equal', 'Random', 'MinVar'});
x = reordercats(x,{'Market Cap','Dividend','Book Value',...
    'Employee','Sales','Equal', 'Random', 'MinVar'});
figure(3)
hold on
bar(x,turnover,'FaceColor',[0 0.3294 0.6235]);
  set(gca,'FontSize',17.6);
title("Average Annual Turnover");
ylabel("Turnover in %");
xlabel("Index");



%% Size Bar Graph

sizes = [MCaveragesize2018, Divaveragesize2018, BVaveragesize2018, Empaveragesize2018,...
    Salesaveragesize2018, Eqaveragesize2018, Randaveragesize2018, MinVaraveragesize2018];
x = categorical({'Market Cap','Dividend','Book Value',...
    'Employee','Sales','Equal', 'Random', 'MinVar'});
x = reordercats(x,{'Market Cap','Dividend','Book Value',...
    'Employee','Sales','Equal', 'Random', 'MinVar'});
figure(4)
hold on
bar(x,sizes,'FaceColor',[0 0.3294 0.6235]);
  set(gca,'FontSize',17.6);
title("Weighted Average Size of Component Stocks");
ylabel("Weighted Size in €bns");
xlabel("Index");





