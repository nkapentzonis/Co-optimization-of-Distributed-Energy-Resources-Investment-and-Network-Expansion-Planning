%% Model
sizing_op = optimproblem('ObjectiveSense','min');
%%%%%%%%% Variables  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
KPV = optimvar('Kpv',sum(npv),1,'Type','integer','LowerBound',0);
KW = optimvar('Kw',sum(nw),1,'Type','integer','LowerBound',0);
KE = optimvar('Ke',sum(nbn),1,'LowerBound',0);
KP = optimvar('Kp',sum(nbn),1,'LowerBound',0);
x = optimvar('x',sum(nbn),T,TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1); % (a.16)
soc = optimvar('soc',sum(nbn),T,TotalNumberOfScenarios,'LowerBound',0);
dis = optimvar('dis',sum(nbn),T,TotalNumberOfScenarios,'LowerBound',0);
ch = optimvar('ch',sum(nbn),T,TotalNumberOfScenarios,'LowerBound',0);
m = optimvar('m',nmg,T,TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1); % (a.37)
ppf = optimvar('Ppf',sum(nbr),T,TotalNumberOfScenarios);
qpf = optimvar('Qpf',sum(nbr),T,TotalNumberOfScenarios);
v = optimvar('V',sum(nn),T,TotalNumberOfScenarios,'LowerBound',repmat(NVMlMat,1,T,TotalNumberOfScenarios),'UpperBound',repmat(NVMuMat,1,T,TotalNumberOfScenarios)); % (a.32)
qmg = optimvar('Q',nmg,T,TotalNumberOfScenarios);
zpv = optimvar('zpv',sum(npv),T,TotalNumberOfScenarios,'LowerBound',0);
zw = optimvar('zw',sum(nw),T,TotalNumberOfScenarios,'LowerBound',0);
pup = optimvar('pup',nmg,T,TotalNumberOfScenarios,'LowerBound',0);
pdn = optimvar('pdn',nmg,T,TotalNumberOfScenarios,'LowerBound',0);
uline = optimvar('uline',sum(ncbr),'Type','integer','LowerBound',0,'UpperBound',1); % (a.2)
yline = optimvar('yline',sum(nbr),TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1); % (a.24)
% bline = optimvar('bline',sum(nbr),TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1);
bline = optimvar('bline',sum(nn)+nmg,sum(nn)+nmg,TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1); % (a.30)

%% CVaR Variables
zeta_obj = optimvar('eta_obj',1);
eta_w_obj = optimvar('sw_obj',1,NumberOfScenarios,'LowerBound',0);

%% Constraints
% Total Line Budget Constraint (a.3)
sizing_op.Constraints.total_line_budget = sum((CLinesCostMat.*CLinesLenMat).*uline)<=TotalLineBudget;

% Total Budget Constraint (a.4)
storage_investment_cost = sum(TC_Energy*KE+TC_Power*KP);
wind_investment_cost = sum(TC_Wind*(KW*wind_blocks));
pv_investment_cost = sum(TC_Solar*(KPV*pv_blocks));
sizing_op.Constraints.total_budget = (storage_investment_cost+wind_investment_cost+pv_investment_cost)<=TotalDERBudget;

% Bounds on Size variables
KW.UpperBound = ceil(KWMAX/wind_blocks); % (a.5)
KPV.UpperBound = ceil(KPVMAX/pv_blocks); % (a.6)
KE.UpperBound = KEMAX; % (a.7)
KP.UpperBound = KPMAX; % (a.8)
% Storage Power-to-Energy ratio Constraint (a.9)
sizing_op.Constraints.energy_to_power_ratio = P_CONST*KP==KE;
% RES Production Bounds
wind_upper = optimconstr(sum(nw),T,TotalNumberOfScenarios);
solar_upper = optimconstr(sum(npv),T,TotalNumberOfScenarios);
for ii = 1:nmg
    for tt=1:T
        for ww=1:TotalNumberOfScenarios
            if ii==1
                solar_upper(1:npv(ii),tt,ww) = zpv(1:npv(ii),tt,ww)<=PVEfficiency*SolarTimeSeries(tt,ii,ww).*(KPV(1:npv(ii))*pv_blocks);
                wind_upper(1:nw(ii),tt,ww) = zw(1:nw(ii),tt,ww)<=WindTimeSeries(tt,ii,ww).*(KW(1:nw(ii))*wind_blocks);
            elseif ii>1
                solar_upper(sum(npv(1:ii-1))+1:sum(npv(1:ii)),tt,ww) = zpv(sum(npv(1:ii-1))+1:sum(npv(1:ii)),tt,ww)<=PVEfficiency*SolarTimeSeries(tt,ii,ww).*(KPV(sum(npv(1:ii-1))+1:sum(npv(1:ii)))*pv_blocks);
                wind_upper(sum(nw(1:ii-1))+1:sum(nw(1:ii)),tt,ww) = zw(sum(nw(1:ii-1))+1:sum(nw(1:ii)),tt,ww)<=WindTimeSeries(tt,ii,ww).*(KW(sum(nw(1:ii-1))+1:sum(nw(1:ii)))*wind_blocks);
            end
        end
    end
end
sizing_op.Constraints.wind_upper = wind_upper; % (a.10)
sizing_op.Constraints.solar_upper = solar_upper; % (a.11)
% Minimum Return Constraint (a.12)
switch minimum_return_choice
    case 1
        WPROD = optimexpr(nmg,T,TotalNumberOfScenarios);
        PVPROD = optimexpr(nmg,T,TotalNumberOfScenarios);
        DISCHARGE = optimexpr(nmg,T,TotalNumberOfScenarios);
        CHARGE = optimexpr(nmg,T,TotalNumberOfScenarios);
        for ii=1:nmg
            for tt=1:T
                for ww=1:TotalNumberOfScenarios
                    if ii==1
                        WPROD(ii,tt,ww) = sum(zw(1:nw(ii),tt,ww));
                        PVPROD(ii,tt,ww) = sum(zpv(1:npv(ii),tt,ww));
                        DISCHARGE(ii,tt,ww) = sum(dis(1:nbn(ii),tt,ww));
                        CHARGE(ii,tt,ww) = sum(ch(1:nbn(ii),tt,ww));
                    elseif ii>1
                        WPROD(ii,tt,ww) = sum(zw(sum(nw(1:ii-1))+1:sum(nw(1:ii)),tt,ww));
                        PVPROD(ii,tt,ww) = sum(zpv(sum(npv(1:ii-1))+1:sum(npv(1:ii)),tt,ww));
                        DISCHARGE(ii,tt,ww) = sum(dis(sum(nbn(1:ii-1))+1:sum(nbn(1:ii)),tt,ww));
                        CHARGE(ii,tt,ww) = sum(ch(sum(nbn(1:ii-1))+1:sum(nbn(1:ii)),tt,ww));
                    end
                end
            end
        end
        TotalWPROD = sum(sum(WPROD,1),2);
        TotalPVPROD = sum(sum(PVPROD,1),2);
        TotalRESOPEX = CW*TotalWPROD+CPV*TotalPVPROD;
        TotalBSUOPEX = sum(sum(CB*(dis+ch),1),2);
        TotalExpectedOPEX = sum(reshape(Possibilities,1,1,TotalNumberOfScenarios).*(TotalRESOPEX+TotalBSUOPEX));
        
        ESP_profits_per_scen = optimexpr(nmg,TotalNumberOfScenarios);
        for gg=1:nmg
            for ss=1:TotalNumberOfScenarios
                ESP_profits_per_scen(gg,ss)= Possibilities(ss)*(sum(Original_LMPs(gg,:,ss).*(WPROD(gg,:,ss)+PVPROD(gg,:,ss)+DISCHARGE(gg,:,ss)-CHARGE(gg,:,ss))));
            end
        end
        if DERSOPEX==1
            sizing_op.Constraints.minimum_return = (-TotalExpectedOPEX+sum(ESP_profits_per_scen(:)))>=return_ratio*(sum(AC_Energy*KE+AC_Power*KP)+sum(AC_Wind*KW*wind_blocks)+sum(AC_Solar*KPV*pv_blocks));
        elseif DERSOPEX==0
            sizing_op.Constraints.minimum_return = (sum(ESP_profits_per_scen(:)))>=return_ratio*(sum(AC_Energy*KE+AC_Power*KP)+sum(AC_Wind*KW*wind_blocks)+sum(AC_Solar*KPV*pv_blocks));
        end
    case 0
        WPROD = optimexpr(nmg,T,TotalNumberOfScenarios);
        PVPROD = optimexpr(nmg,T,TotalNumberOfScenarios);
        DISCHARGE = optimexpr(nmg,T,TotalNumberOfScenarios);
        CHARGE = optimexpr(nmg,T,TotalNumberOfScenarios);
        for ii=1:nmg
            for tt=1:T
                for ww=1:TotalNumberOfScenarios
                    if ii==1
                        WPROD(ii,tt,ww) = sum(zw(1:nw(ii),tt,ww));
                        PVPROD(ii,tt,ww) = sum(zpv(1:npv(ii),tt,ww));
                        DISCHARGE(ii,tt,ww) = sum(dis(1:nbn(ii),tt,ww));
                        CHARGE(ii,tt,ww) = sum(ch(1:nbn(ii),tt,ww));
                    elseif ii>1
                        WPROD(ii,tt,ww) = sum(zw(sum(nw(1:ii-1))+1:sum(nw(1:ii)),tt,ww));
                        PVPROD(ii,tt,ww) = sum(zpv(sum(npv(1:ii-1))+1:sum(npv(1:ii)),tt,ww));
                        DISCHARGE(ii,tt,ww) = sum(dis(sum(nbn(1:ii-1))+1:sum(nbn(1:ii)),tt,ww));
                        CHARGE(ii,tt,ww) = sum(ch(sum(nbn(1:ii-1))+1:sum(nbn(1:ii)),tt,ww));
                    end
                end
            end
        end
        TotalWPROD = sum(sum(WPROD,1),2);
        TotalPVPROD = sum(sum(PVPROD,1),2);
        TotalRESOPEX = CW*TotalWPROD+CPV*TotalPVPROD;
        TotalBSUOPEX = sum(sum(CB*(dis+ch),1),2);
        TotalExpectedOPEX = sum(reshape(Possibilities,1,1,TotalNumberOfScenarios).*(TotalRESOPEX+TotalBSUOPEX));
        
        ESP_profits_per_scen = optimexpr(nmg,TotalNumberOfScenarios);
        for gg=1:nmg
            for ss=1:TotalNumberOfScenarios
                ESP_profits_per_scen(gg,ss)= Possibilities(ss)*(sum(Original_LMPs(gg,:,ss).*(WPROD(gg,:,ss)+PVPROD(gg,:,ss)+DISCHARGE(gg,:,ss)-CHARGE(gg,:,ss))));
            end
        end
        
    otherwise
        error('Wrong choice for the incorporation of a minimum return on annualized investment cost!!')
end
% Power Discharge/Charge Upper Bound (a.13)
sizing_op.Constraints.dis_upper = dis<=repmat(KP,1,T,TotalNumberOfScenarios);
sizing_op.Constraints.ch_upper = ch<=repmat(KP,1,T,TotalNumberOfScenarios);
% Storage Binary Variable Constraints
sizing_op.Constraints.dis_x = dis<=BIGX.*x; % (a.14)
sizing_op.Constraints.ch_x = ch<=BIGX.*(1-x); % (a.15)
% State of Energy Dynamics (a.17 - a.18)
sizing_op.Constraints.soc_dynamicst2 = soc(:,2:T,:) == soc(:,1:T-1,:)-dis(:,2:T,:)./repmat(DischargeEff,1,T-1,TotalNumberOfScenarios)+ch(:,2:T,:).*repmat(ChargeEff,1,T-1,TotalNumberOfScenarios);
sizing_op.Constraints.soc_dynamicst1 = soc(:,1,:) == repmat(InitialSOCPerMat,1,1,TotalNumberOfScenarios).*repmat(KE,1,1,TotalNumberOfScenarios)-dis(:,1,:)./repmat(DischargeEff,1,1,TotalNumberOfScenarios)+ch(:,1,:).*repmat(ChargeEff,1,1,TotalNumberOfScenarios);
% State of Energy Upper Bound (a.19)
sizing_op.Constraints.soc_upper_bound = soc<=repmat(KE,1,T,TotalNumberOfScenarios);
% State of Energy Last Timeslot (a.20)
sizing_op.Constraints.soc_last_timeslot = soc(:,T,:)>=repmat(FinalSOCPerMat,1,1,TotalNumberOfScenarios).*repmat(KE,1,1,TotalNumberOfScenarios);
% DN Active and Reactive Power Balance (a.21), (a.22)
CoefActivePF = cell(1,nmg);
CoefBatteryLoc = cell(1,nmg);
CoefReactivePF = cell(1,nmg);
for mm=1:nmg
    Line = Lines{mm};
    MGCoefActivePF = zeros(nn(mm),nbr(mm));
    MGCoefReactivePF = zeros(nn(mm),nbr(mm));
    for ii=1:nn(mm)
        MGCoefActivePF(ii,Line((ToBusMG{mm}==ii))) = -1;
        MGCoefActivePF(ii,Line((FromBusMG{mm}==ii))) = 1;
        MGCoefReactivePF(ii,Line((ToBusMG{mm}==ii))) = -1;
        MGCoefReactivePF(ii,Line((FromBusMG{mm}==ii))) = 1;
    end
    CoefActivePF{mm} = MGCoefActivePF;
    CoefReactivePF{mm} = MGCoefReactivePF;
end

active_power_flow = optimconstr(sum(nn),T,TotalNumberOfScenarios);
reactive_power_flow = optimconstr(sum(nn),T,TotalNumberOfScenarios);
for ii=1:nmg
    for tt=1:T
        for ww=1:TotalNumberOfScenarios
            if ii==1
                active_power_flow(1:nn(ii),tt,ww) = CoefActivePF{ii}*ppf(1:nbr(ii),tt,ww) - ResLoc{ii}*zw(1:nw(ii),tt,ww) - PVLoc{ii}*zpv(1:npv(ii),tt,ww)- BatLoc{ii}*dis(1:nbn(ii),tt,ww) + BatLoc{ii}*ch(1:nbn(ii),tt,ww) == -InfLoadMat(1:nn(ii),tt,ww);
                reactive_power_flow(1:nn(ii),tt,ww) = CoefReactivePF{ii}*qpf(1:nbr(ii),tt,ww) - (ResLoc{ii}*zw(1:nw(ii),tt,ww).*WindPF(1)) - ((PVLoc{ii}*zpv(1:npv(ii),tt,ww)).*SolarPF(1)) == -aInfMat(1:nn(ii),tt).*(InfLoadMat(1:nn(ii),tt,ww));
            elseif ii>1
                active_power_flow(sum(nn(1:ii-1))+1:sum(nn(1:ii)),tt,ww) = CoefActivePF{ii}*ppf(sum(nbr(1:ii-1))+1:sum(nbr(1:ii)),tt,ww) - ResLoc{ii}*zw(sum(nw(1:ii-1))+1:sum(nw(1:ii)),tt,ww) - PVLoc{ii}*zpv(sum(npv(1:ii-1))+1:sum(npv(1:ii)),tt,ww) - BatLoc{ii}*dis(sum(nbn(1:ii-1))+1:sum(nbn(1:ii)),tt,ww) + BatLoc{ii}*ch(sum(nbn(1:ii-1))+1:sum(nbn(1:ii)),tt,ww) == -InfLoadMat(sum(nn(1:ii-1))+1:sum(nn(1:ii)),tt,ww);
                reactive_power_flow(sum(nn(1:ii-1))+1:sum(nn(1:ii)),tt,ww) = CoefReactivePF{ii}*qpf(sum(nbr(1:ii-1))+1:sum(nbr(1:ii)),tt,ww) - (ResLoc{ii}*(zw(sum(nw(1:ii-1))+1:sum(nw(1:ii)),tt,ww))).*WindPF(1) - (PVLoc{ii}*(zpv(sum(npv(1:ii-1))+1:sum(npv(1:ii)),tt,ww))).*SolarPF(1)  == -aInfMat(sum(nn(1:ii-1))+1:sum(nn(1:ii)),tt).*(InfLoadMat(sum(nn(1:ii-1))+1:sum(nn(1:ii)),tt,ww));
            end
        end
    end
end


sizing_op.Constraints.active_power_flow_power_balance = active_power_flow;
sizing_op.Constraints.reactive_power_flow_power_balance = reactive_power_flow;

% Voltage-Branch Flow (a.23)
CoefNodalV = cell(1,nmg);
CoefppfNodalV = cell(1,nmg);
CoefqpfNodalV = cell(1,nmg);
Vparam = cell(1,nmg);
for mm=1:nmg
    Line = Lines{mm};
    LineResistance = LinesR{mm};
    LineReactance = LinesX{mm};
    DestBus = ToBusMG{mm};
    SourceBus = FromBusMG{mm};
    MGppfNodalV = diag(2*LineResistance);
    MGqpfNodalV = diag(2*LineReactance);
    MGNodalV = zeros(nbr(mm),nn(mm));
    for ii=1:nbr(mm)
        MGNodalV(ii,DestBus(ii)) = 1;
        if SourceBus(ii)>0
            MGNodalV(ii,SourceBus(ii)) = -1;
        end
    end
    CoefNodalV{mm} = MGNodalV;
    CoefppfNodalV{mm} = MGppfNodalV;
    CoefqpfNodalV{mm} = MGqpfNodalV;
    Vparam{mm} = SourceBus==0;
end

voltage_branch_lower = optimconstr(sum(nbr),T,TotalNumberOfScenarios);
voltage_branch_upper = optimconstr(sum(nbr),T,TotalNumberOfScenarios);
for ii=1:nmg
    for tt=1:T
        for ww=1:TotalNumberOfScenarios
            if ii==1
                voltage_branch_upper(1:nbr(ii),tt,ww) = (CoefNodalV{ii}*v(1:nn(ii),tt,ww)+ CoefppfNodalV{ii}*ppf(1:nbr(ii),tt,ww) + CoefqpfNodalV{ii}*qpf(1:nbr(ii),tt,ww) - double(Vparam{ii}))<=Mline*(1-yline(1:nbr(ii),ww));
                voltage_branch_lower(1:nbr(ii),tt,ww) = -(CoefNodalV{ii}*v(1:nn(ii),tt,ww)+ CoefppfNodalV{ii}*ppf(1:nbr(ii),tt,ww) + CoefqpfNodalV{ii}*qpf(1:nbr(ii),tt,ww) - double(Vparam{ii}))<=Mline*(1-yline(1:nbr(ii),ww));
            elseif ii>1
                voltage_branch_upper(sum(nbr(1:ii-1))+1:sum(nbr(1:ii)),tt,ww) = (CoefNodalV{ii}*v(sum(nn(1:ii-1))+1:sum(nn(1:ii)),tt,ww)+ CoefppfNodalV{ii}*ppf(sum(nbr(1:ii-1))+1:sum(nbr(1:ii)),tt,ww) + CoefqpfNodalV{ii}*qpf(sum(nbr(1:ii-1))+1:sum(nbr(1:ii)),tt,ww) - double(Vparam{ii}))<=Mline*(1-yline(sum(nbr(1:ii-1))+1:sum(nbr(1:ii)),ww));
                voltage_branch_lower(sum(nbr(1:ii-1))+1:sum(nbr(1:ii)),tt,ww) = -(CoefNodalV{ii}*v(sum(nn(1:ii-1))+1:sum(nn(1:ii)),tt,ww)+ CoefppfNodalV{ii}*ppf(sum(nbr(1:ii-1))+1:sum(nbr(1:ii)),tt,ww) + CoefqpfNodalV{ii}*qpf(sum(nbr(1:ii-1))+1:sum(nbr(1:ii)),tt,ww) - double(Vparam{ii}))<=Mline*(1-yline(sum(nbr(1:ii-1))+1:sum(nbr(1:ii)),ww));
            end
        end
    end
end
sizing_op.Constraints.voltage_branch_upper = voltage_branch_upper;
sizing_op.Constraints.voltage_branch_lower = voltage_branch_lower;

% Constraints (a.25 - a.26) that forbid the usage of a non-installed line
a25 = optimconstr(sum(nebr),TotalNumberOfScenarios);
a26 = optimconstr(sum(ncbr),TotalNumberOfScenarios);
for ii=1:nmg
    for ss=1:TotalNumberOfScenarios
        if ii==1
            a25(1:nebr(ii),ss) = yline(1:nebr(ii),ss)<=1;
            a26(1:ncbr(ii),ss) = yline(nebr(ii)+1:nbr(ii),ss)<=uline(1:ncbr(ii));
        elseif ii>1
            a25(sum(nebr(1:ii-1))+1:sum(nebr(1:ii)),ss) = yline(sum(nbr(1:ii-1))+1:sum(nbr(1:ii-1))+sum(nebr(ii)),ss)<=1;
            a26(sum(ncbr(1:ii-1))+1:sum(ncbr(1:ii)),ss) = yline(sum(nbr(1:ii-1))+sum(nebr(ii))+1:sum(nbr(1:ii)),ss)<=uline(sum(ncbr(1:ii-1))+1:sum(ncbr(1:ii)));
        end
    end
end
sizing_op.Constraints.a25 = a25;% Existing lines can always be utilized (is it necessary????? I think NO!)
sizing_op.Constraints.a26 = a26;

% Constraint (a.27)
a27 = optimconstr(sum(nbr),TotalNumberOfScenarios);
for ss=1:TotalNumberOfScenarios
    %     for ii=1:nmg
    for ll=1:sum(nbr)
        a27(ll,ss) = yline(ll,ss) == bline(FromNodesMat(ll),ToNodesMat(ll),ss)+bline(ToNodesMat(ll),FromNodesMat(ll),ss);
    end
    %     end
end
sizing_op.Constraints.a27 = a27;
%%%%%%% SHOULD A CONSTRAINT THAT NULLIFIES ? NO-LINE ELEMENTS??? %%%%%%%
% (a.28) is not necessary (?)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constraint (a.29)
SourceNodes = zeros(1,nmg);
for ii=1:nmg
    if ii==1
        temp = FromBusMG{ii}+1;
        SourceNodes(ii) = temp(1);
    elseif ii>1
        temp = FromBusMG{ii}+sum(nn(ii-1))+ii;
        SourceNodes(ii) = temp(1);
    end
end
BCoef = zeros(sum(nn)+nmg);
a29 = optimconstr(sum(nn),TotalNumberOfScenarios);
for ii=1:sum(nn)+nmg
    if ~ismember(ii,SourceNodes)
        BCoef(FromNodesMat(LinesMat(ToNodesMat==ii)),ii) = 1;
    end
end
AllBs = sum(bline(:,:,1).*BCoef,1);
AllBs(SourceNodes) = [];
for ss=1:TotalNumberOfScenarios
    a29(:,ss) = AllBs<=ones(1,sum(nn)); % == OR <= ?????
end
sizing_op.Constraints.a29 = a29;

% Line Flow - Apparent Power Upper Bound (a.31)
switch line_capacity_constraint_option
    case 1
        ppf_upper = optimconstr(sum(nbr),T,TotalNumberOfScenarios);
        ppf_lower = optimconstr(sum(nbr),T,TotalNumberOfScenarios);
        qpf_upper = optimconstr(sum(nbr),T,TotalNumberOfScenarios);
        qpf_lower = optimconstr(sum(nbr),T,TotalNumberOfScenarios);
        for ss=1:TotalNumberOfScenarios
            for tt=1:T
                ppf_upper(:,tt,ss) = ppf(:,tt,ss)<=LinesPuMat.*(yline(:,ss));
                ppf_lower(:,tt,ss) = -ppf(:,tt,ss)<=-LinesPlMat.*(yline(:,ss));
                qpf_upper(:,tt,ss) = qpf(:,tt,ss)<=LinesQuMat.*(yline(:,ss));
                qpf_lower(:,tt,ss) = -qpf(:,tt,ss)<=-LinesQlMat.*(yline(:,ss));
            end
        end
        sizing_op.Constraints.ppf_upper = ppf_upper;
        sizing_op.Constraints.ppf_lower = ppf_lower;
        sizing_op.Constraints.qpf_upper = qpf_upper;
        sizing_op.Constraints.qpf_lower = qpf_lower;
    case 2
        apparent_bound = optimconstr(sum(nbr),T,MM,TotalNumberOfScenarios);
        for ii=1:sum(nbr)
            for tt=1:T
                for ww=1:TotalNumberOfScenarios
                    %                     apparent_bound(ii,tt,:,ww) = A_theta{ii}.*repmat(ppf(ii,tt,ww),MM,1) + B_theta{ii}.*repmat(qpf(ii,tt,ww),MM,1)<=C_theta((ii-1)*MM+1:ii*MM);
                    apparent_bound(ii,tt,:,ww) = A_theta{ii}.*repmat(ppf(ii,tt,ww),MM,1) + B_theta{ii}.*repmat(qpf(ii,tt,ww),MM,1)<=C_theta((ii-1)*MM+1:ii*MM).*repmat(yline(ii,ww),MM,1);
                end
            end
        end
        sizing_op.Constraints.dn_flows_constraint = apparent_bound;
    otherwise
        error('There is no such an option for line capacity constraint!');
end

% Constraint (a.33) - SEEMS REDUNDANT IF CONSTRAINT (a.29) IS AN EQUALITY CONSTRAINT!!!!!!!!!!
% %%%%%%%%%%%%%%%%%%%%% NEEDS FIXING!!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % nn(ii) is wrong. The number of the connected nodes should be there and
% % not the number of the total existing and candidate nodes!!!!!!!!!
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a33 = optimconstr(nmg,TotalNumberOfScenarios);
% for ss=1:TotalNumberOfScenarios
%     for ii=1:nmg
%         if ii==1
%             a33(ii,ss) = sum(yline(1:nbr(ii),ss))==nn(ii);
%         elseif ii>1
%             a33(ii,ss) = sum(yline(sum(nbr(1:ii-1))+1:sum(nbr(1:ii)),ss))==nn(ii);
%         end
%     end
% end
% % % % % % sizing_op.Constraints.a33 = a33;
% Active Power Balance at substation (a.34) and Reactive Power Balance at substation
substation_ppb = optimconstr(nmg,T,TotalNumberOfScenarios);
substation_qpb = optimconstr(nmg,T,TotalNumberOfScenarios);
for ii=1:nmg
    if ii==1
        substation_ppb(ii,:,:)= (sum(repmat(double(Vparam{ii}),1,T,TotalNumberOfScenarios).*ppf(1:nbr(ii),:,:))) + pup(ii,:,:) - pdn(ii,:,:)==0;
        substation_qpb(ii,:,:)= (sum(repmat(double(Vparam{ii}),1,T,TotalNumberOfScenarios).*qpf(1:nbr(ii),:,:))) - qmg(ii,:,:)==0;
    else
        substation_ppb(ii,:,:)= sum(repmat(double(Vparam{ii}),1,T,TotalNumberOfScenarios).*ppf(sum(nbr(1:ii-1))+1:sum(nbr(1:ii)),:,:))+pup(ii,:,:) - pdn(ii,:,:) ==0;
        substation_qpb(ii,:,:)= sum(repmat(double(Vparam{ii}),1,T,TotalNumberOfScenarios).*qpf(sum(nbr(1:ii-1))+1:sum(nbr(1:ii)),:,:))-qmg(ii,:,:)==0;
    end
end
sizing_op.Constraints.substation_active_power_balance = substation_ppb;
sizing_op.Constraints.substation_reactive_power_balance = substation_qpb;

% Forbidding simultaneous draw/upply power from/to the transmission grid
sizing_op.Constraints.offer_binary_bound = pup<=m.*repmat(SubPuMat,1,T,TotalNumberOfScenarios); % (a.35)
sizing_op.Constraints.bid_binary_bound = pdn<=(1-m).*repmat(SubPuMat,1,T,TotalNumberOfScenarios); % (a.36)

%% Objective
DERS_CAPEX = sum(AC_Wind*KW*wind_blocks)+sum(AC_Solar*KPV*pv_blocks)+sum(AC_Energy*KE)+sum(AC_Power*KP);
LINES_CAPEX = sum((AC_Lines.*CLinesLenMat).*uline);
if DERSOPEX==1
    OPEX = TotalExpectedOPEX - sum(reshape(Possibilities,1,1,TotalNumberOfScenarios).*sum(sum(Original_LMPs.*(pup-pdn))));
    OPEX_per_scenario = optimexpr(1,NumberOfScenarios);
    for ss=1:TotalNumberOfScenarios
        OPEX_per_scenario(ss) = TotalRESOPEX(ss)+TotalBSUOPEX(ss)+sum(sum(Original_LMPs(:,:,ss).*(pup(:,:,ss)-pdn(:,:,ss))));
    end
elseif DERSOPEX==0
    OPEX = -sum(reshape(Possibilities,1,1,TotalNumberOfScenarios).*sum(sum(Original_LMPs.*(pup-pdn))));
    OPEX_per_scenario = optimexpr(1,NumberOfScenarios);
    for ss=1:TotalNumberOfScenarios
        OPEX_per_scenario(ss) = sum(sum(Original_LMPs(:,:,ss).*(pup(:,:,ss)-pdn(:,:,ss))));
    end    
end

CVAR_OBJ = -(zeta_obj - (1/(1-alpha_cvar))*(Possibilities*eta_w_obj'));
sizing_op.Constraints.cvar_cons1 = zeta_obj-eta_w_obj<=OPEX_per_scenario;
% sizing_op.Constraints.cvar_cons2 = (1/(1-alpha_cvar))*(Possibilities*sw_obj')<=eta_obj;
sizing_op.Objective = DERS_CAPEX+LINES_CAPEX+(1-beta_cvar)*OPEX+beta_cvar*CVAR_OBJ;


