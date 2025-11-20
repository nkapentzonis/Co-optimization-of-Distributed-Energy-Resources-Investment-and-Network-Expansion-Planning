rng(1)
RepDays = 1;
%% Number Of Year Days
if DataYear==2015
    NumberOfDays = 365;
elseif DataYear == 2016
    NumberOfDays = 366;
end
%% System-Wide and Load Consumption (scaling factors from ENTSO-E)

Yearly_Demand_Vec = xlsread(['ENTSOE_Load_Data_Greece_',num2str(DataYear),'.xlsx'],'Sheet1','A1:A8784');

Yearly_Demand_Matrix = reshape(Yearly_Demand_Vec,24,NumberOfDays)';
% Adjust to T timeslots
Final_Yearly_Demand_Matrix = reshape(sum((reshape(Yearly_Demand_Matrix',24/T,T,NumberOfDays)),1),T,NumberOfDays)';

% Scaling Factors
Maximum_Demand_Per_Load_Scenario = max(Final_Yearly_Demand_Matrix,[],2)';
YDSF = Final_Yearly_Demand_Matrix./repmat(Maximum_Demand_Per_Load_Scenario',1,T);
YDSF(isnan(YDSF)) = 0;
Yearly_Demand_Scaling_Factors =  YDSF;

%% Solar Power (scaling factors from VIMSEN data)
%%%%%%%%%%%%%%%%% 2015 Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PVList2015 = [1 6 17 25 33 48];
% PVList2016 = [1 6 17 25 33 38];
PVList2015 = [17 1 33];
PVList2016 = [17 1 33];
if DataYear==2015
    Yearly_Solar_Matrix_96 = cell(length(PVList2015));
    Yearly_Solar_Matrix = zeros(365,24,length(PVList2015));
    Final_Yearly_Solar_Matrix = zeros(365,T,length(PVList2015));
    it = 0;
    Maximum_Solar_Per_Scenario = zeros(length(PVList2016),NumberOfDays);
    Yearly_Solar_Scaling_Factors = zeros(NumberOfDays,T,length(PVList2016));
    
    for ii  = PVList2015
        it = it + 1;
        Yearly_Solar_Matrix_96{it} = xlsread(['VIMSEN_PV_Data_Elefsina',num2str(PVList2015(it)),'_2015.xlsx'],'Sheet1','A1:CR365');
        % From 96 to 24 timeslots
        for tt=1:24
            YSM96 = Yearly_Solar_Matrix_96{it};
            Yearly_Solar_Matrix(:,tt,it) = YSM96(:,4*tt)+YSM96(:,4*tt-1) + YSM96(:,4*tt-2) + YSM96(:,4*tt-3);
        end
        Final_Yearly_Solar_Matrix(:,:,it) = reshape(sum((reshape(Yearly_Solar_Matrix(:,:,it)',24/T,T,365)),1),T,365)';
        Maximum_Solar_Per_Scenario(it,:) = max(Final_Yearly_Solar_Matrix(:,:,it),[],2)';
        YSSF = Final_Yearly_Solar_Matrix(:,:,it)./repmat(Maximum_Solar_Per_Scenario(it,:)',1,T);
        YSSF(isnan(YSSF))=0;
        Yearly_Solar_Scaling_Factors(:,:,it) = YSSF;
    end
elseif DataYear==2016
%%%%%%%%%%%%%%%%%% 2016 Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Yearly_Solar_Vec_96 = cell(length(PVList2016));
    Yearly_Solar_Matrix_96 = cell(length(PVList2016));
    it = 0;
    Yearly_Solar_Matrix = zeros(335,24,length(PVList2016));
    Final_Yearly_Solar_Matrix = zeros(335,T,length(PVList2016));
    Maximum_Solar_Per_Scenario = zeros(length(PVList2016),335);
    Yearly_Solar_Scaling_Factors = zeros(335,T,length(PVList2016));
    for ii = PVList2016
        it = it + 1;
        Yearly_Solar_Vec_96{it} = xlsread('PVs_20kW_VIMSEN_2016.xlsx','PVs_20kW_VIMSEN_2016',xlRC2A1(3,PVList2016(it)+1,32162,PVList2016(it)+1));
        Yearly_Solar_Matrix_96{it} = reshape(Yearly_Solar_Vec_96{it},4*24,335)';
        for tt=1:24
            YSM96 = Yearly_Solar_Matrix_96{it};
            Yearly_Solar_Matrix(:,tt,it) = YSM96(:,4*tt)+YSM96(:,4*tt-1) + YSM96(:,4*tt-2) + YSM96(:,4*tt-3);
        end
        Final_Yearly_Solar_Matrix(:,:,it) = reshape(sum((reshape(Yearly_Solar_Matrix(:,:,it)',24/T,T,335)),1),T,335)';
        Maximum_Solar_Per_Scenario(it,:) = max(Final_Yearly_Solar_Matrix(:,:,it),[],2)';
        YSSF =Final_Yearly_Solar_Matrix(:,:,it)./repmat(Maximum_Solar_Per_Scenario(it,:)',1,T);
        YSSF(isnan(YSSF)) = 0;
        Yearly_Solar_Scaling_Factors(:,:,it) =  YSSF;
    end
end


%% Wind Power (scaling factors from VIMSEN data)
WindList2015 = [17 1 33];
WindList2016 = [17 1 33];
if DataYear == 2015
    it = 0;
    Yearly_Wind_Matrix_96 = cell(length(WindList2015));
    Yearly_Wind_Matrix = zeros(NumberOfDays,24,length(WindList2015));
    Final_Yearly_Wind_Matrix = zeros(NumberOfDays,T,length(WindList2015));
    Maximum_Wind_Per_Scenario = zeros(length(WindList2015),NumberOfDays);
    Yearly_Wind_Scaling_Factors = zeros(NumberOfDays,T,length(WindList2015));
    for ii = WindList2015
        it = it + 1;
        Yearly_Wind_Matrix_96{it} = xlsread('aiolika_MV_2015.xlsx','Sheet1',xlRC2A1((WindList2015(it)-1)*NumberOfDays+1, 4, WindList2015(it)*NumberOfDays, 99));
        for tt=1:24
            YWM96 = Yearly_Wind_Matrix_96{it};
            Yearly_Wind_Matrix(:,tt,it) = YWM96(:,4*tt)+YWM96(:,4*tt-1) + YWM96(:,4*tt-2) + YWM96(:,4*tt-3);
        end
        Final_Yearly_Wind_Matrix(:,:,it) = reshape(sum((reshape(Yearly_Wind_Matrix(:,:,it)',24/T,T,365)),1),T,365)'; 
        Maximum_Wind_Per_Scenario(it,:) = max(Final_Yearly_Wind_Matrix(:,:,it),[],2)';
        YWSF = Final_Yearly_Wind_Matrix(:,:,it)./repmat(Maximum_Wind_Per_Scenario(it,:)',1,T);
        YWSF(isnan(YWSF)) = 0;
        Yearly_Wind_Scaling_Factors(:,:,it) = YWSF;
    end
elseif DataYear==2016
    it = 0;
    Yearly_Wind_Matrix_96 = cell(length(WindList2016));
    Yearly_Wind_Matrix = zeros(274,24,length(WindList2016));
    Final_Yearly_Wind_Matrix = zeros(274,T,length(WindList2016));
    Maximum_Wind_Per_Scenario = zeros(length(WindList2016),274);
    Yearly_Wind_Scaling_Factors = zeros(274,T,length(WindList2016));
    for ii = WindList2016
        it = it + 1;
        Yearly_Wind_Matrix_96{it} = xlsread('aiolika_MV_2016.xlsx','Sheet2',xlRC2A1((WindList2015(it)-1)*274+1,4,WindList2015(it)*274,99));
        for tt=1:24
            YWM96 = Yearly_Wind_Matrix_96{it};
            Yearly_Wind_Matrix(:,tt,it) = YWM96(:,4*tt)+YWM96(:,4*tt-1) + YWM96(:,4*tt-2) + YWM96(:,4*tt-3);
        end
        Final_Yearly_Wind_Matrix(:,:,it) = reshape(sum((reshape(Yearly_Wind_Matrix(:,:,it)',24/T,T,274)),1),T,274)';
        Maximum_Wind_Per_Scenario(it,:) = max(Final_Yearly_Wind_Matrix(:,:,it),[],2)';
        YWSF = Final_Yearly_Wind_Matrix(:,:,it)./repmat(Maximum_Wind_Per_Scenario(it,:)',1,T);
        YWSF(isnan(YWSF))=0;
        Yearly_Wind_Scaling_Factors(:,:,it) = YWSF;
    end
end

%% Clustering
% SolarLocation = 1;
% WindLocation = 1;
% NumberOfScenarios = 1;
if RepDays
    NumberOfRepresentativeDays = NumberOfScenarios;
    switch DataYear
        case 2015
            SizingScenarios = [Yearly_Demand_Scaling_Factors];
            for kk=1:length(WindList2015)
                SizingScenarios = [SizingScenarios,Yearly_Wind_Scaling_Factors(:,:,kk),Yearly_Solar_Scaling_Factors(:,:,kk)];
            end
        case 2016
            SizingScenarios = [Yearly_Demand_Scaling_Factors];
            for kk=1:length(WindList2015)
                SizingScenarios = [SizingScenarios,Yearly_Wind_Scaling_Factors(:,:,kk),Yearly_Solar_Scaling_Factors(:,:,kk)];
            end
        otherwise
            error('Wrong Data Year!!\n')
    end
    [ids,RepresentativeDays,s,d] = kmeans(SizingScenarios,NumberOfRepresentativeDays);
    RepDaysPoss = zeros(1,NumberOfRepresentativeDays);
    for ss=1:NumberOfRepresentativeDays
        RepDaysPoss(ss) = sum(ids == ss);
    end
    
    SizingLoadScenarios = RepresentativeDays(:,1:T);
    SizingWindScenarios = cell(1,length(WindList2015));
    SizingSolarScenarios = cell(1,length(WindList2015));
    for gg = 1:length(WindList2015)
        SizingWindScenarios{gg} = RepresentativeDays(:,(2*gg-1)*T+1:2*gg*T);
        SizingSolarScenarios{gg} = RepresentativeDays(:,2*gg*T+1:2*gg*T+T);
    end
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Load Scenarios %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     NumberOfLoadScenarios = NumberOfScenarios;
    [idx_d,Clusters_d,s_d,d_d] = kmeans(Yearly_Demand_Scaling_Factors,NumberOfLoadScenarios);
    LoadScenPoss = zeros(1,NumberOfLoadScenarios);
    for ss=1:NumberOfLoadScenarios
        LoadScenPoss(ss) = sum(idx_d==ss);
    end
    
    NumberOfWindScenarios = NumberOfScenarios;
    NumberOfSolarScenarios = NumberOfScenarios;
    switch DataYear
        case 2015
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PV Scenarios %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [idx_pv,Clusters_pv,s_pv,d_pv] = kmeans(Yearly_Solar_Scaling_Factors(:,:,PVList2015(SolarLocation)),NumberOfSolarScenarios);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Wind Scenarios %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [idx_w,Clusters_w,s_w,d_w] = kmeans(Yearly_Wind_Scaling_Factors(:,:,WindList2015(WindLocation)),NumberOfWindScenarios);
        case 2016
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PV Scenarios %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [idx_pv,Clusters_pv,s_pv,d_pv] = kmeans(Yearly_Solar_Scaling_Factors(:,:,PVList2016(SolarLocation)),NumberOfSolarScenarios);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Wind Scenarios %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [idx_w,Clusters_w,s_w,d_w] = kmeans(Yearly_Wind_Scaling_Factors(:,:,WindList2016(WindLocation)),NumberOfWindScenarios);            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Wind Scenarios %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        otherwise
            error('Wrong Year!!\n');
    end
    SolarScenPoss = zeros(1,NumberOfSolarScenarios);
    for ss=1:NumberOfSolarScenarios
        SolarScenPoss(ss) = sum(idx_pv==ss);
    end
    WindScenPoss = zeros(1,NumberOfWindScenarios);
    for ss=1:NumberOfWindScenarios
        WindScenPoss(ss) = sum(idx_w==ss);
    end  

end
