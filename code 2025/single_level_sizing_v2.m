%% TEP is added compared to version 0. Also VoLL is added.
clear,clc,close all
VoLL = 1000;

%% Input Data (excel) Path
TNDataPath = 'C:\Users\nikol\OneDrive\Dokumente\Σχολή\Thesis\models\'; % Transmission Network Data PC location
DNDataPath = 'C:\Users\nikol\OneDrive\Dokumente\Σχολή\Thesis\models\'; % Distribution Network Data excel file path

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INVESTMENT MODEL ASSUMP  
MVA = 100;

%%%%%%% Investment Blocks %%%%%%%%%%%%%%
wind_blocks = 0.05/MVA; % ModelFD-50000 https://inkpv.com/product/solar-and-wind-power-system/50kw-wind-power-plant-cost/
pv_blocks = 0.05/MVA; % 50KW (an industrial area of 250 m^2 is considered as a block)

%%%%%%% Minimum ROI Constraint
% Select if a minimum return on annualized investment cost is added to the model
% Select 1: DERs' clear profits
% Select 0: No
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minimum_return_choice = 0; % Value: 0 or 1!!

%%%%%%% Value-at-Risk %%%%%%%%%%%%%%%%%%%%%%
alpha_cvar = 0.90;
beta_cvar = 0.10;

%%%%%%% Distribution Network Lines' Capacity
% Select the form of the line capacity constraint:
% Select 1: Box Constraints
% Select 2: Quadratic Constraint
line_capacity_constraint_option = 1;
% Number of polygon sizes for the circle linearization
MM = 8;

%%%%%%% Select if DERs Operating Costs are included
% Select 0: No
% Select 1: Yes
DERSOPEX = 1;

%%%%%%% Select if the Users' Bills or the aggregate DN welfare is minimized
% Select 1: Users' Bills
% Select 2: Aggregate DN Welfare
ModelObj = 2; % IN THE SINGLE-LEVEL MODEL, THE CHOICE SHOULD BE 2!!!

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUT PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Budget for TEP
TotalLineBudget = 1000000000;

%%%%%%%%%%% Budget for DER investments
TotalDERBudget = 1000000000;

%%%%%%%%%%% Minimum ROI Constraint
% Minimum Return Ratio
return_ratio = 1.20;

%%%%%%%%%%% Selection of Transmission and Distribution Systems
% List of Test Systems
%%%%%%% 6-Bus Test System:   30 MVA power base %%%%%%%%%%%%%%%%%%
%%%%%%% 24-Bus Test System:   100 MVA power base %%%%%%%%%%%%%%%%%%
%%%%%%% 118-Bus Test System:  100 MVA power base %%%%%%%%%%%%%%%%%%
Nasrolahpour_6BusTestSystem = 'Sizing_6BusTestSystem_v312';                        % 1
IEEE_Reliability_Test_System = 'Sizing_24BusTestSystem_v312';                      % 2
IEEE_118_Test_System = 'Sizing_118BusTestSystem_v312';                             % 3
Zone8_NE_Test_System = 'Sizing_8Zone_v312';                                        % 4
TestSystemList = {Nasrolahpour_6BusTestSystem,IEEE_Reliability_Test_System,IEEE_118_Test_System,Zone8_NE_Test_System};

% Select The Test System from the List Above
TestSystemSelection = 1;

% List of Microgrid Test Systems
DN_15Bus = '15-Node DTS';                                        % 1
DN_33Bus = '33-Node DTS';                                        % 2
DN_69Bus = '69-Node DTS';                                        % 3
Test_MG = 'Testing';                                             % 4
DN_33Bus_Alternative = '33-Node Alter';                          % 5
MGTestSystemList = {DN_15Bus,DN_33Bus, DN_69Bus, Test_MG, DN_33Bus_Alternative};

% Select The MG Test System from the List Above
MGTestSystemSelection = 1;

% Transmission Grid Reference Bus
RefBus = 1;

% Number of Timeslots
T = 24;

% Distribution Network Voltage Base
KV = 11;

% Distribution Network Impedance Base
Zb = (KV^2)/MVA;

%%%%%%%%%%%%%%% Scenarios %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Number of Representative Days
NumberOfScenarios = 4;

% Data Year (2015 or 2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DataYear = 2015; % Should NOT BE CHANGED!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% GenCos' price offers Scenarios
GenOffersScenarios = 1;
GenOffersPossibilities = 1;

%%%%%%%%%%%%%%%%% DERs' operating parameters %%%%%%%%%%%%%%%%%%%%%%%
% Maximum DERs' Capacities
KPVMAX = 0;
KWMAX = 0;
KEMAX = 10/MVA;
KPMAX = 10/MVA;

% Energy-to-Power Ratio Constants
P_CONST = 6;

% PV output efficiency factor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%DERSIncome
PVEfficiency = 0.95;  % Value<1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input Charging and Discharging Efficiencies
Charge_Efficiency = 0.93;
Discharge_Efficiency = 0.93;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Investment Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Batteries %%%%%%%%%%%
%%% Interest rate
% 0.05 -> Hassan2018, Dvorkin2018 (DRvsES), Pandzic2018

%%% Yearlife
% 10 -> Hassan2018, Dvorkin2018 (DRvsES), Pandzic2018
% 15 -> Pandzic_batter_pv_investment

%%%%%%% Solar PV %%%%%%%%%%%%
%%% Interest rate
% 0.03 -> Pandzic_batter_pv_investment

%%% Yearlife
% 20 -> Pandzic_batter_pv_investment, Bhattacharjee2021

%%% Interest rate
rate_of_return = [0.05;0.05;0.05;0.05];
lines_rate_of_return = 0.05;
%%% Yearlife
yearlife = [15;15;15;15];
lines_yearlife = 20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DERs' Buying Costs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Solar %%%%%%%%%%
% $1750/kW  Pandzic_batter_pv_investment

%%%%%%%%% Wind %%%%%%%%%%
% $1000/kW  Baringo2012 & BMGLoaringo2014
% $1300/kW  https://weatherguardwind.com/how-much-does-wind-turbine-cost-worth-it/

%%%%%%%%% Storage %%%%%%%%%%
% $400/kWh | $400/kW Pandzic_batter_pv_investment
% $20/kWh  | $500/kW Hassan2018 & Dvorkin2018(DRvsES) - Low Cost Scenario
% $75/kWh  | $1300/kW Hassan2018 & Dvorkin2018(DRvsES) - Medium Cost Scenario
% $150/kWh | $2600/kW Hassan2018 & Dvorkin2018(DRvsES) - High Cost Scenario
% $291/kWh | $616/kW Bhattacharjee2021

% DERs' Buying Costs
TC_Solar = 890000;%830000;
TC_Wind = 1350000;%1300000;
TC_Energy = 25000;%20000;
TC_Power = 550000;%500000;

% DERS OPEX
CW = 3.5; % Mortaz2020
CPV = 2.5;  % Mortaz2020
CB = 0.5; % Nasrolahpour2020

%%%%%%%%%%%%%%%%%%%%%%%%%% Big-M Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = 2000;
Mline = 100000;
BIGX = 100000;

% DAM Price Cap
Price_Cap = 1000;

%% Generate Scenarios
GenerateScenarios

%% Input from Excel Files
InputFromExcel_v2

%% Piece-Wise Linearization
%PieceWiseCircleLinearization

%% Annualized DERs' Investment Costs
AC_Solar = annualize_investment_cost(rate_of_return(1),yearlife(1),TC_Solar);
AC_Wind = annualize_investment_cost(rate_of_return(2),yearlife(2),TC_Wind);
AC_Energy = annualize_investment_cost(rate_of_return(3),yearlife(3),TC_Energy);
AC_Power = annualize_investment_cost(rate_of_return(4),yearlife(4),TC_Power);
AC_Lines = annualize_investment_cost(lines_rate_of_return,lines_yearlife, CLinesCostMat);
%% Calculate Original LMPs (Prices without investments) that will be used in the single-level model
% Calculate_Original_LMPs_v1
Calculate_Original_LMPs_v2
%% Investment Model
Create_Single_Level_TEP_Investment_Model_new

%% Solve Model
Solve_Single_Level_Investment_Model

%% Incorporate the new investments in the Market
Run_Market_v2

%% Results
% Investments in DERs
Solar_Size = MVA*x_opt.Kpv*pv_blocks;
Wind_Size = MVA*x_opt.Kw*wind_blocks;
Energy_Size = MVA*x_opt.Ke;
Power_Size = MVA*x_opt.Kp;

Solar_Size_pu = x_opt.Kpv*pv_blocks;
Wind_Size_pu = x_opt.Kw*wind_blocks;
Energy_Size_pu = x_opt.Ke;
Power_Size_pu = x_opt.Kp;

% Investment in Lines
NewLines = nonzeros(x_opt.uline.*(1:sum(ncbr))');
NetworkConfiguration = cell(1,TotalNumberOfScenarios);
for ss=1:TotalNumberOfScenarios
    NetworkConfiguration{ss} = nonzeros(x_opt.yline(:,ss).*[1:sum(nbr)]');
end

% Annualized DERs Investment Cost
AIC = AC_Solar*sum(Solar_Size) + AC_Wind*sum(Wind_Size)+AC_Energy*sum(Energy_Size)+AC_Power*sum(Power_Size);
AIC_pu = AC_Solar*sum(Solar_Size_pu) + AC_Wind*sum(Wind_Size_pu)+AC_Energy*sum(Energy_Size_pu)+AC_Power*sum(Power_Size_pu);

% Annualized Lines Investment Cost
ALIC = sum((AC_Lines.*CLinesLenMat).*x_opt.uline);

% Total DERs Investment Cost
TIC = sum(Solar_Size)*TC_Solar+sum(Wind_Size)*TC_Wind + sum(Energy_Size)*TC_Energy + sum(Power_Size)*TC_Power;
TIC_pu = sum(Solar_Size_pu)*TC_Solar+sum(Wind_Size_pu)*TC_Wind + sum(Energy_Size_pu)*TC_Energy + sum(Power_Size_pu)*TC_Power;

% Total Lines Investment Cost
TLIC = sum((CLinesCostMat.*CLinesLenMat).*x_opt.uline);
%% DSO
% Consumers' Bills
FinalLoad = InfLoadMat;
DERSIncome = zeros(nmg,TotalNumberOfScenarios);
DERSIncome_pu = zeros(nmg,TotalNumberOfScenarios);
TotalLoadBill = zeros(nmg,TotalNumberOfScenarios);
TotalLoadBill_pu = zeros(nmg,TotalNumberOfScenarios);
CLoadBill = zeros(nmg,TotalNumberOfScenarios);
CLoadBill_pu = zeros(nmg,TotalNumberOfScenarios);
DERSIncome2 = zeros(nmg,TotalNumberOfScenarios);
DERSIncome2_pu = zeros(nmg,TotalNumberOfScenarios);
TotalLoadBill2 = zeros(nmg,TotalNumberOfScenarios);
TotalLoadBill2_pu = zeros(nmg,TotalNumberOfScenarios);
CLoadBill2 = zeros(nmg,TotalNumberOfScenarios);
CLoadBill2_pu = zeros(nmg,TotalNumberOfScenarios);
Prices_s = [];
for ss=1:TotalNumberOfScenarios
    Prices_s = Original_LMPs(:,:,ss);
    for ii = 1:nmg
        nodal_prices = Prices_s(ii,:);
        nodal_prices2 = MGLoc(:,ii)'*duals2{ss}.Constraints.system_power_balance;
        if ii==1
            DERSIncome(ii,ss) = Possibilities(ss)*sum(nodal_prices.*(sum(x_opt.zw(1:nw(ii),:,ss)*MVA,1)+sum((MVA*x_opt.dis(1:nbn(ii),:,ss)-MVA*x_opt.ch(1:nbn(ii),:,ss)),1)+sum(x_opt.zpv(1:npv(ii),:,ss)*MVA,1)));
            DERSIncome2(ii,ss) = Possibilities(ss)*sum(nodal_prices2.*(sum(x_opt.zw(1:nw(ii),:,ss)*MVA,1)+sum((MVA*x_opt.dis(1:nbn(ii),:,ss)-MVA*x_opt.ch(1:nbn(ii),:,ss)),1)+sum(x_opt.zpv(1:npv(ii),:,ss)*MVA,1)));
            DERSIncome_pu(ii,ss) = Possibilities(ss)*sum(nodal_prices.*(sum(x_opt.zw(1:nw(ii),:,ss),1)+sum((x_opt.dis(1:nbn(ii),:,ss)-x_opt.ch(1:nbn(ii),:,ss)),1)+sum(x_opt.zpv(1:npv(ii),:,ss),1)));
            DERSIncome2_pu(ii,ss) = Possibilities(ss)*sum(nodal_prices2.*(sum(x_opt.zw(1:nw(ii),:,ss),1)+sum((x_opt.dis(1:nbn(ii),:,ss)-x_opt.ch(1:nbn(ii),:,ss)),1)+sum(x_opt.zpv(1:npv(ii),:,ss),1)));
            TotalLoadBill(ii,ss) = Possibilities(ss)*sum(nodal_prices.*sum(MVA*FinalLoad(1:nn(ii),:,ss),1));
            CLoadBill(ii,ss) = Possibilities(ss)*sum(nodal_prices.*sum(MVA*ConnectedLoadMat(1:nn(ii),:,ss),1));
            TotalLoadBill2(ii,ss) = Possibilities(ss)*sum(nodal_prices2.*sum(MVA*FinalLoad(1:nn(ii),:,ss),1));
            CLoadBill2(ii,ss) = Possibilities(ss)*sum(nodal_prices2.*sum(MVA*ConnectedLoadMat(1:nn(ii),:,ss),1));
            TotalLoadBill_pu(ii,ss) = Possibilities(ss)*sum(nodal_prices.*sum(FinalLoad(1:nn(ii),:,ss),1));
            CLoadBill_pu(ii,ss) = Possibilities(ss)*sum(nodal_prices.*sum(ConnectedLoadMat(1:nn(ii),:,ss),1));
            TotalLoadBill2_pu(ii,ss) = Possibilities(ss)*sum(nodal_prices2.*sum(FinalLoad(1:nn(ii),:,ss),1));
            CLoadBill2_pu(ii,ss) = Possibilities(ss)*sum(nodal_prices2.*sum(ConnectedLoadMat(1:nn(ii),:,ss),1));
        else
            DERSIncome(ii,ss) = Possibilities(ss)*sum(nodal_prices.*(sum(x_opt.zw(sum(nw(1:ii-1))+1:sum(nw(1:ii)),:,ss)*MVA,1)+sum((MVA*x_opt.dis(sum(nbn(1:ii-1))+1:sum(nbn(1:ii)),:,ss)-MVA*x_opt.ch(sum(nbn(1:ii-1))+1:sum(nbn(1:ii)),:,ss)),1)+sum(x_opt.zpv(sum(npv(1:ii-1))+1:sum(npv(1:ii)),:,ss)*MVA,1)));
            DERSIncome2(ii,ss) = Possibilities(ss)*sum(nodal_prices2.*(sum(x_opt.zw(sum(nw(1:ii-1))+1:sum(nw(1:ii)),:,ss)*MVA,1)+sum((MVA*x_opt.dis(sum(nbn(1:ii-1))+1:sum(nbn(1:ii)),:,ss)-MVA*x_opt.ch(sum(nbn(1:ii-1))+1:sum(nbn(1:ii)),:,ss)),1)+sum(x_opt.zpv(sum(npv(1:ii-1))+1:sum(npv(1:ii)),:,ss)*MVA,1)));
            DERSIncome_pu(ii,ss) = Possibilities(ss)*sum(nodal_prices.*(sum(x_opt.zw(sum(nw(1:ii-1))+1:sum(nw(1:ii)),:,ss),1)+sum((x_opt.dis(sum(nbn(1:ii-1))+1:sum(nbn(1:ii)),:,ss)-x_opt.ch(sum(nbn(1:ii-1))+1:sum(nbn(1:ii)),:,ss)),1)+sum(x_opt.zpv(sum(npv(1:ii-1))+1:sum(npv(1:ii)),:,ss),1)));
            DERSIncome2_pu(ii,ss) = Possibilities(ss)*sum(nodal_prices2.*(sum(x_opt.zw(sum(nw(1:ii-1))+1:sum(nw(1:ii)),:,ss),1)+sum((x_opt.dis(sum(nbn(1:ii-1))+1:sum(nbn(1:ii)),:,ss)-x_opt.ch(sum(nbn(1:ii-1))+1:sum(nbn(1:ii)),:,ss)),1)+sum(x_opt.zpv(sum(npv(1:ii-1))+1:sum(npv(1:ii)),:,ss),1)));
            TotalLoadBill(ii,ss) = Possibilities(ss)*sum(nodal_prices.*sum(MVA*FinalLoad(sum(nn(1:ii-1))+1:sum(nn(1:ii)),:,ss),1));
            TotalLoadBill2(ii,ss) = Possibilities(ss)*sum(nodal_prices2.*sum(MVA*FinalLoad(sum(nn(1:ii-1))+1:sum(nn(1:ii)),:,ss),1));
            TotalLoadBill_pu(ii,ss) = Possibilities(ss)*sum(nodal_prices.*sum(FinalLoad(sum(nn(1:ii-1))+1:sum(nn(1:ii)),:,ss),1));
            TotalLoadBill2_pu(ii,ss) = Possibilities(ss)*sum(nodal_prices2.*sum(FinalLoad(sum(nn(1:ii-1))+1:sum(nn(1:ii)),:,ss),1));
            
            CLoadBill(ii,ss) = Possibilities(ss)*sum(nodal_prices.*sum(MVA*ConnectedLoadMat(sum(nn(1:ii-1))+1:sum(nn(1:ii)),:,ss),1));
            CLoadBill2(ii,ss) = Possibilities(ss)*sum(nodal_prices2.*sum(MVA*ConnectedLoadMat(sum(nn(1:ii-1))+1:sum(nn(1:ii)),:,ss),1));
            CLoadBill_pu(ii,ss) = Possibilities(ss)*sum(nodal_prices.*sum(ConnectedLoadMat(sum(nn(1:ii-1))+1:sum(nn(1:ii)),:,ss),1));
            CLoadBill2_pu(ii,ss) = Possibilities(ss)*sum(nodal_prices2.*sum(ConnectedLoadMat(sum(nn(1:ii-1))+1:sum(nn(1:ii)),:,ss),1));            
        end
    end
end

% Expected OPEX
if DERSOPEX==1
    TotalWindProd = sum(sum(x_opt.zw*MVA,1),2);
    TotalWindProd_pu = sum(sum(x_opt.zw,1),2);
    TotalSolarProd = sum(sum(x_opt.zpv*MVA,1),2);
    TotalSolarProd_pu = sum(sum(x_opt.zpv,1),2);
    RenewablesOPEX = CW*TotalWindProd+CPV*TotalSolarProd;
    RenewablesOPEX_pu = CW*TotalWindProd_pu+CPV*TotalSolarProd_pu;
    BatteriesOPEX = sum(sum(CB*(MVA*x_opt.dis+MVA*x_opt.ch),1),2);
    BatteriesOPEX_pu = sum(sum(CB*(x_opt.dis+x_opt.ch),1),2);
    ExpectedOPEX = sum(reshape(Possibilities,1,1,TotalNumberOfScenarios).*(RenewablesOPEX+BatteriesOPEX));
    ExpectedOPEX_pu = sum(reshape(Possibilities,1,1,TotalNumberOfScenarios).*(RenewablesOPEX_pu+BatteriesOPEX_pu));
end

% DSO Yearly Market Profits
DSO_Market_Profits_per_Scenario = zeros(1,TotalNumberOfScenarios);
DSO_Market_Profits_per_Scenario2 = zeros(1,TotalNumberOfScenarios);
DSO_Market_Profits_per_Scenario_pu = zeros(1,TotalNumberOfScenarios);
DSO_Market_Profits_per_Scenario2_pu = zeros(1,TotalNumberOfScenarios);
for ss=1:TotalNumberOfScenarios
    DSO_Market_Profits_per_Scenario(ss) = sum(sum((MVA*x_opt.pup(:,:,ss)-MVA*x_opt.pdn(:,:,ss)).*(Original_LMPs(:,:,ss))));
    
    DSO_Market_Profits_per_Scenario2(ss) = sum(sum((MVA*x_TSO{ss}.pup-MVA*x_TSO{ss}.pdn).*(MGLoc'*duals2{ww}.Constraints.system_power_balance)));
    DSO_Market_Profits_per_Scenario2_pu(ss) = sum(sum((x_TSO{ss}.pup-x_TSO{ss}.pdn).*(MGLoc'*duals2{ww}.Constraints.system_power_balance)));
    
    DSO_Market_Profits_per_Scenario_pu(ss) = sum(sum((x_opt.pup(:,:,ss)-x_opt.pdn(:,:,ss)).*(Original_LMPs(:,:,ss))));
    
end
Expected_DSO_Market_Profits = Possibilities*DSO_Market_Profits_per_Scenario';
Expected_DSO_Market_Profits2 = Possibilities*DSO_Market_Profits_per_Scenario2';
Expected_DSO_Market_Profits_pu = Possibilities*DSO_Market_Profits_per_Scenario_pu';
Expected_DSO_Market_Profits2_pu = Possibilities*DSO_Market_Profits_per_Scenario2_pu';

% DSO Yearly Market Net Profits (minus OPEX)
if DERSOPEX==1
    DSO_Net_Expected_Market_Profits = Expected_DSO_Market_Profits - ExpectedOPEX;
    DSO_Net_Expected_Market_Profits2 = Expected_DSO_Market_Profits2 - ExpectedOPEX;
    DSO_Net_Expected_Market_Profits_pu = Expected_DSO_Market_Profits_pu - ExpectedOPEX_pu;
    DSO_Net_Expected_Market_Profits2_pu = Expected_DSO_Market_Profits2_pu - ExpectedOPEX_pu;
elseif DERSOPEX==0
    DSO_Net_Expected_Market_Profits = Expected_DSO_Market_Profits;
    DSO_Net_Expected_Market_Profits2 = Expected_DSO_Market_Profits2;
    DSO_Net_Expected_Market_Profits_pu = Expected_DSO_Market_Profits_pu;
    DSO_Net_Expected_Market_Profits2_pu = Expected_DSO_Market_Profits2_pu;
end

% DSO Yearly Total Profits (Market Net Profits minus Annualized Investment Cost)
Yearly_Expected_DSO_Clear_Profits = DSO_Net_Expected_Market_Profits - AIC-ALIC;
Yearly_Expected_DSO_Clear_Profits2 = DSO_Net_Expected_Market_Profits2 - AIC-ALIC;
Yearly_Expected_DSO_Clear_Profits_pu = DSO_Net_Expected_Market_Profits_pu - AIC_pu-ALIC;
Yearly_Expected_DSO_Clear_Profits2_pu = DSO_Net_Expected_Market_Profits2_pu - AIC_pu-ALIC;

% Consumers' Bills
Yearly_Expected_Consumers_Bills = sum(TotalLoadBill(:));
Yearly_Expected_Consumers_Bills2 = sum(TotalLoadBill2(:));
Yearly_Expected_Consumers_Bills_pu = sum(TotalLoadBill_pu(:));
Yearly_Expected_Consumers_Bills2_pu = sum(TotalLoadBill2_pu(:));

Yearly_Expected_Con_Consumers_Bills = sum(CLoadBill(:));
Yearly_Expected_Con_Consumers_Bills2 = sum(CLoadBill2(:));
Yearly_Expected_Con_Consumers_Bills_pu = sum(CLoadBill_pu(:));
Yearly_Expected_Con_Consumers_Bills2_pu = sum(CLoadBill2_pu(:));

% Check if Objective Function differs from DSO Net Profits
if ((-Yearly_Expected_DSO_Clear_Profits_pu-neg_profits)/neg_profits)>0.00001
    fprintf('----------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------\n')
    fprintf('Objective Function Value and DSO yearly net profits do not match!!\n');
    fprintf('----------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------\n')
end

% Value of Lost Load
VOLL =0;
%% ESP
% ESP Yearly Market Profits
Expected_ESP_Market_Profits = sum(DERSIncome(:));
Expected_ESP_Market_Profits_pu = sum(DERSIncome_pu(:));
Expected_ESP_Market_Profits2 = sum(DERSIncome2(:));
Expected_ESP_Market_Profits2_pu = sum(DERSIncome2_pu(:));
% ESP Yearly Market Net Profits (minus OPEX)
if DERSOPEX==1
    Expected_ESP_Net_Market_Profits = Expected_ESP_Market_Profits-ExpectedOPEX;
    Expected_ESP_Net_Market_Profits_pu = Expected_ESP_Market_Profits_pu-ExpectedOPEX_pu;
    Expected_ESP_Net_Market_Profits2 = Expected_ESP_Market_Profits2-ExpectedOPEX;
    Expected_ESP_Net_Market_Profits2_pu = Expected_ESP_Market_Profits2_pu-ExpectedOPEX_pu;
elseif DERSOPEX==0
    Expected_ESP_Net_Market_Profits = Expected_ESP_Market_Profits;
    Expected_ESP_Net_Market_Profits_pu = Expected_ESP_Market_Profits_pu;
    Expected_ESP_Net_Market_Profits2 = Expected_ESP_Market_Profits2;
    Expected_ESP_Net_Market_Profits2_pu = Expected_ESP_Market_Profits2_pu;    
end

% ESP Yearly Total Profits (Market Net Profits minus Annualized Investment Cost)
Yearly_Expected_ESP_Clear_Profits = Expected_ESP_Net_Market_Profits - AIC;
Yearly_Expected_ESP_Clear_Profits_pu = Expected_ESP_Net_Market_Profits_pu - AIC_pu;
Yearly_Expected_ESP_Clear_Profits2 = Expected_ESP_Net_Market_Profits2 - AIC;
Yearly_Expected_ESP_Clear_Profits2_pu = Expected_ESP_Net_Market_Profits2_pu - AIC_pu;

% Check ESP total profits vs (DSO total profits + Users' Bills)
if abs(Expected_ESP_Net_Market_Profits_pu - (DSO_Net_Expected_Market_Profits_pu+Yearly_Expected_Consumers_Bills_pu))/Expected_ESP_Net_Market_Profits_pu>0.00001
    fprintf('----------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------\n')
    fprintf('ESP Profits and DSO Profits Plus Bills do not match!!!\n')
    fprintf('----------------------------------------------------------------\n')
    fprintf('----------------------------------------------------------------\n')
end
% Check Minimum Recovery Factor
if minimum_return_choice == 2
    if Expected_ESP_Net_Market_Profits_pu-return_ratio*AIC_pu<-0.000001
        fprintf('----------------------------------------------------------------\n')
        fprintf('----------------------------------------------------------------\n')
        fprintf('ESP Profits do not respect the minimum return constraint!!!\n')
        fprintf('----------------------------------------------------------------\n')
        fprintf('----------------------------------------------------------------\n')
    end
end

%% TSO
% Expected Social Welfare & Generation Cost
TSOSocialWelfare = zeros(1,TotalNumberOfScenarios);
TSOSocialWelfare_pu = zeros(1,TotalNumberOfScenarios);
GenerationCost = zeros(1,TotalNumberOfScenarios);
GenerationCost_pu = zeros(1,TotalNumberOfScenarios);
for ss=1:TotalNumberOfScenarios
    TSOSocialWelfare_pu(ss) = -BF(ss);
    TSOSocialWelfare(ss) = -(sum(sum(GenBids(:,:,ss).*x_TSO{ss}.g*MVA))-sum(sum(DemBids(:,:,ss).*x_TSO{ss}.d*MVA))+sum(sum(cpup(:,:,ss).*x_TSO{ss}.pup*MVA))-sum(sum(cpdn(:,:,ss).*x_TSO{ss}.pdn*MVA)));
    GenerationCost(ss) = sum(sum(GenBids(:,:,ss).*x_TSO{ss}.g*MVA))+sum(sum(cpup(:,:,ss).*x_TSO{ss}.pup*MVA));
    GenerationCost_pu(ss) = sum(sum(GenBids(:,:,ss).*x_TSO{ss}.g))+sum(sum(cpup(:,:,ss).*x_TSO{ss}.pup));

end

ExpectedTSOSocialWelfare = Possibilities*TSOSocialWelfare';
ExpectedTSOSocialWelfare_pu = Possibilities*TSOSocialWelfare_pu';


ExpectedGenerationCost = Possibilities*GenerationCost';
ExpectedGenerationCost_pu = Possibilities*GenerationCost_pu';

%% RES Spillage
PV_Spillage = zeros(nmg,T,TotalNumberOfScenarios);
Wind_Spillage = zeros(nmg,T,TotalNumberOfScenarios);
for ii = 1:nmg
    for tt=1:T
        for ww=1:TotalNumberOfScenarios
            if ii==1
                PV_Spillage(1:npv(ii),tt,ww) = Possibilities(ww)*(PVEfficiency*SolarTimeSeries(tt,ii,ww).*Solar_Size(1:npv(ii))-x_opt.zpv(1:npv(ii),tt,ww)*MVA);
                Wind_Spillage(1:nw(ii),tt,ww) = Possibilities(ww)*(WindTimeSeries(tt,ii,ww).*Wind_Size(1:nw(ii))-x_opt.zw(1:nw(ii),tt,ww)*MVA);
            elseif ii>1
                PV_Spillage(sum(npv(1:ii-1))+1:sum(npv(1:ii)),tt,ww) = Possibilities(ww)*(PVEfficiency*SolarTimeSeries(tt,ii,ww).*Solar_Size(sum(npv(1:ii-1))+1:sum(npv(1:ii)))-x_opt.zpv(sum(npv(1:ii-1))+1:sum(npv(1:ii)),tt,ww)*MVA);
                Wind_Spillage(sum(nw(1:ii-1))+1:sum(nw(1:ii)),tt,ww) = Possibilities(ww)*(WindTimeSeries(tt,ii,ww).*Wind_Size(sum(nw(1:ii-1))+1:sum(nw(1:ii)))-x_opt.zw(sum(nw(1:ii-1))+1:sum(nw(1:ii)),tt,ww)*MVA);
            end
        end
    end
end

RES_Spillage = sum(PV_Spillage(:))+sum(Wind_Spillage(:));
RES_Spillage

%% Print Results
SocialWelfareIncreasePercentage = ((ExpectedTSOSocialWelfare_pu-SW0)/SW0)*100;
UserBillsDecreasePercentage = ((Original_Bills-Yearly_Expected_Con_Consumers_Bills_pu)/Original_Bills)*100;
UserBillsDecreasePercentage2 = ((Original_Bills-Yearly_Expected_Con_Consumers_Bills2_pu)/Original_Bills)*100;
HypoUserBillsDecreasePercentage = ((Original_Hypo_Bills-Yearly_Expected_Consumers_Bills_pu)/Original_Hypo_Bills)*100;
HypoUserBillsDecreasePercentage2 = ((Original_Hypo_Bills-Yearly_Expected_Consumers_Bills2_pu)/Original_Hypo_Bills)*100;
ESPNetProfitPercentage = ((Expected_ESP_Net_Market_Profits-AIC)/AIC)*100;
ESPNetProfitPercentage2 = ((Expected_ESP_Net_Market_Profits2-AIC)/AIC)*100;
GenerationCostReductionPercentage = ((-ExpectedGenerationCost_pu+GC0)/GC0)*100;
ValueOfLostLoadReductionPercentage = ((VOLL0-VOLL)/VOLL0)*100;
fprintf('Changes in DSO, ESP and TSO objectives\n')
fprintf('--------------------------------------\n')
fprintf('RES Spillage: %2.2f% MWh\n')
fprintf('ESP: %2.2f%% investors'' profit percentage vs %2.2f%% desired profit percentage\n',ESPNetProfitPercentage,100*(return_ratio-1))
fprintf('DSO: %2.2f%% decrease in Consumers'' Bills\n',UserBillsDecreasePercentage)
fprintf('DSO: Hypothetical %2.2f%% decrease in Consumers'' Bills\n',HypoUserBillsDecreasePercentage)
fprintf('DSO: %2.2f%% decrease in Value of Lost Load\n',ValueOfLostLoadReductionPercentage)
fprintf('TSO: %2.2f%% increaase in Social Welfare\n',SocialWelfareIncreasePercentage)
fprintf('TSO: %2.2f%% decrease in Generation Cost\n',GenerationCostReductionPercentage)
fprintf('--------------------------------------\n')
fprintf('Version2 ESP: %2.2f%% investors'' profit percentage vs %2.2f%% desired profit percentage\n',ESPNetProfitPercentage2,100*(return_ratio-1))
fprintf('Version2 DSO: %2.2f%% decrease in Consumers'' Bills\n',UserBillsDecreasePercentage2)
fprintf('Version2 DSO: Hypothetical %2.2f%% decrease in Consumers'' Bills\n',HypoUserBillsDecreasePercentage2)

%% Print DER and Line Investment Totals
fprintf('\n----- Investment Totals -----\n')
fprintf('Total DER Investments:\n')
fprintf('  Solar PV: %2.2f MW across %d locations\n', sum(Solar_Size), sum(Solar_Size > 0))
fprintf('  Wind: %2.2f MW across %d locations\n', sum(Wind_Size), sum(Wind_Size > 0))
fprintf('  Energy Storage Capacity: %2.2f MWh across %d locations\n', sum(Energy_Size), sum(Energy_Size > 0))
fprintf('  Energy Storage Power: %2.2f MW across %d locations\n', sum(Power_Size), sum(Power_Size > 0))
fprintf('  Total DER Investment Cost: $%2.2f\n', TIC)
fprintf('  Annualized DER Investment Cost: $%2.2f\n', AIC)

fprintf('\nTotal Line Investments:\n')
fprintf('  Number of New Lines: %d\n', length(NewLines))
fprintf('  Total Line Investment Cost: $%2.2f\n', TLIC)
fprintf('  Annualized Line Investment Cost: $%2.2f\n', ALIC)

fprintf('\nTotal Investment (DER + Lines):\n')
fprintf('  Total Investment Cost: $%2.2f\n', TIC + TLIC)
fprintf('  Total Annualized Investment Cost: $%2.2f\n', AIC + ALIC)
fprintf('--------------------------\n\n')

% Generate investment report
generate_investment_report(Solar_Size, Wind_Size, Energy_Size, Power_Size, NewLines, TIC, TLIC, AIC, ALIC, minimum_return_choice, nmg, PVLoc, ResLoc, BatLoc);

% Generate power flow report
generate_power_flow_report(ppf, x_opt);

function result = annualize_investment_cost(r,y,TC)
capital_recovery_factor = (r*((1+r)^y))/((1+r)^y-1);
result = TC*capital_recovery_factor;
end


function generate_investment_report(Solar_Size, Wind_Size, Energy_Size, Power_Size, NewLines, TIC, TLIC, AIC, ALIC, minimum_return_choice, nmg, PVLoc, ResLoc, BatLoc)
    % Create timestamp for unique filename
    timestamp = datestr(now, 'yyyy-mm-dd_HH-MM');
    filename = sprintf('investment_report_%s.txt', timestamp);
    
    % Open file for writing
    fid = fopen(filename, 'w');
    
    % Write header
    fprintf(fid, '=================================================\n');
    fprintf(fid, 'DER and Line Investment Report\n');
    fprintf(fid, 'Generated on: %s\n', datestr(now));
    if minimum_return_choice == 1
        fprintf(fid, 'Model Version: With ESP Profit Constraint\n');
    else
        fprintf(fid, 'Model Version: Without ESP Profit Constraint\n');
    end
    fprintf(fid, '=================================================\n\n');
    
    % Define eligible nodes for each microgrid
    mg1_pv_nodes = [1, 2, 5, 6, 9, 16, 17, 19];
    mg2_pv_nodes = [2, 3, 6, 11, 12, 14, 17, 18];
    mg1_wind_nodes = [1, 3, 4, 11, 12, 16, 18, 19];
    mg2_wind_nodes = [1, 6, 7, 8, 9, 10, 16, 18];
    
    % Track current index for each type of DER
    pv_start_idx = 1;
    wind_start_idx = 1;
    battery_start_idx = 1;
    
    % Iterate through each microgrid
    for mg = 1:nmg
        fprintf(fid, '\nMICROGRID %d\n', mg);
        fprintf(fid, '=================\n\n');
        
        % Solar PV placements
        fprintf(fid, 'Solar PV Installations:\n');
        if mg == 1
            pv_nodes = mg1_pv_nodes;
        else
            pv_nodes = mg2_pv_nodes;
        end
        pv_end_idx = pv_start_idx + length(pv_nodes) - 1;
        mg_solar_size = Solar_Size(pv_start_idx:min(pv_end_idx, length(Solar_Size)));
        
        for i = 1:length(pv_nodes)
            if i <= length(mg_solar_size) && mg_solar_size(i) > 0
                fprintf(fid, '  Node %d: %.2f MW\n', pv_nodes(i), mg_solar_size(i));
            end
        end
        fprintf(fid, 'Microgrid %d Total Solar PV: %.2f MW\n\n', mg, sum(mg_solar_size));
        pv_start_idx = pv_end_idx + 1;
        
        % Wind placements
        fprintf(fid, 'Wind Installations:\n');
        if mg == 1
            wind_nodes = mg1_wind_nodes;
        else
            wind_nodes = mg2_wind_nodes;
        end
        wind_end_idx = wind_start_idx + length(wind_nodes) - 1;
        mg_wind_size = Wind_Size(wind_start_idx:min(wind_end_idx, length(Wind_Size)));
        
        for i = 1:length(wind_nodes)
            if i <= length(mg_wind_size) && mg_wind_size(i) > 0
                fprintf(fid, '  Node %d: %.2f MW\n', wind_nodes(i), mg_wind_size(i));
            end
        end
        fprintf(fid, 'Microgrid %d Total Wind: %.2f MW\n\n', mg, sum(mg_wind_size));
        wind_start_idx = wind_end_idx + 1;
        
        % Energy Storage placements
        fprintf(fid, 'Energy Storage Installations:\n');
        battery_nodes = 1:19; % Both microgrids have nodes 1-19 eligible for batteries
        battery_end_idx = battery_start_idx + length(battery_nodes) - 1;
        mg_energy_size = Energy_Size(battery_start_idx:min(battery_end_idx, length(Energy_Size)));
        mg_power_size = Power_Size(battery_start_idx:min(battery_end_idx, length(Power_Size)));
        
        for i = 1:length(battery_nodes)
            if i <= length(mg_energy_size) && (mg_energy_size(i) > 0 || mg_power_size(i) > 0)
                fprintf(fid, '  Node %d: %.2f MWh capacity, %.2f MW power\n', ...
                    battery_nodes(i), mg_energy_size(i), mg_power_size(i));
            end
        end
        fprintf(fid, 'Microgrid %d Total Storage: %.2f MWh capacity, %.2f MW power\n\n', ...
            mg, sum(mg_energy_size), sum(mg_power_size));
        battery_start_idx = battery_end_idx + 1;
    end
    
    % New Lines
    fprintf(fid, '\nNew Transmission Lines:\n');
    fprintf(fid, '-------------------------\n');
    fprintf(fid, 'Number of New Lines: %d\n', length(NewLines));
    for i = 1:length(NewLines)
        fprintf(fid, 'Line %d\n', NewLines(i));
    end
    fprintf(fid, '\n');
    
    % Investment Summary
    fprintf(fid, '\nInvestment Summary:\n');
    fprintf(fid, 'Total DER Investment: $%.2f\n', TIC);
    fprintf(fid, 'Total Line Investment: $%.2f\n', TLIC);
    fprintf(fid, 'Total Investment: $%.2f\n', (TIC + TLIC));
    fprintf(fid, 'Annualized DER Investment: $%.2f/year\n', AIC);
    fprintf(fid, 'Annualized Line Investment: $%.2f/year\n', ALIC);
    fprintf(fid, 'Total Annualized Investment: $%.2f/year\n', (AIC + ALIC));

    % Add new metrics
    fprintf(fid, '\nInvestment Ratios:\n');
    fprintf(fid, 'Lines to DER Ratio (LtDR): %.2f%%\n', (TLIC/TIC)*100);
    fprintf(fid, 'Lines to Total Ratio (LtTR): %.2f%%\n', (TLIC/(TIC + TLIC))*100);
    
    % Close file
    fclose(fid);
    
    fprintf('Investment report generated: %s\n', filename);
end

% Generate detailed investment report
generate_investment_report(Solar_Size, Wind_Size, Energy_Size, Power_Size, NewLines, TIC, TLIC, AIC, ALIC, minimum_return_choice, nmg, PVLoc, ResLoc, BatLoc);