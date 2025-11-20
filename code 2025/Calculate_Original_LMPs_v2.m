%% Model
original_bilevel = optimproblem('ObjectiveSense','min');
%% Variables
KPV = optimvar('Kpv',sum(npv),1,'LowerBound',0); % PV Size Variable
KW = optimvar('Kw',sum(nw),1,'LowerBound',0); % Wind Size Variable
KE = optimvar('Ke',sum(nbn),1,'LowerBound',0); % Energy Size Variable (batteries)
KP = optimvar('Kp',sum(nbn),1,'LowerBound',0); % Power Size Variable (batteries)
x = optimvar('x',sum(nbn),T,TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1); % Binary variable indicating charging or discharging power (batteries) (a.19)
soc = optimvar('soc',sum(nbn),T,TotalNumberOfScenarios,'LowerBound',0); % State of energy variable (batteries)
dis = optimvar('dis',sum(nbn),T,TotalNumberOfScenarios,'LowerBound',0); % Batteries discharge variable
ch = optimvar('ch',sum(nbn),T,TotalNumberOfScenarios,'LowerBound',0); % Batterues charge variable
m = optimvar('m',nmg,T,TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1); % Binary variable indicating DN selling or buying power (a.15)
ppf = optimvar('Ppf',sum(nbr),T,TotalNumberOfScenarios); % DN active power flow variable
qpf = optimvar('Qpf',sum(nbr),T,TotalNumberOfScenarios); % DN reactive power flow variable
v = optimvar('V',sum(nn),T,TotalNumberOfScenarios,'LowerBound',repmat(NVMlMat,1,T,TotalNumberOfScenarios),'UpperBound',repmat(NVMuMat,1,T,TotalNumberOfScenarios)); % DN voltage variable - voltage bounds (a.35)
qmg = optimvar('Q',nmg,T,TotalNumberOfScenarios); % Reactive power traded between TN and DN  
o = optimvar('o',nmg,T,TotalNumberOfScenarios,'LowerBound',0); % DN offer variable 
b = optimvar('b',nmg,T,TotalNumberOfScenarios,'LowerBound',0); % DN bid variable
g = optimvar('g',ng,T,TotalNumberOfScenarios,'LowerBound',repmat(GMIN,1,1,TotalNumberOfScenarios),'UpperBound',repmat(GMAX,1,1,TotalNumberOfScenarios)); % (TN) Generation variable - generation bounds (b.3)
d = optimvar('d',nd,T,TotalNumberOfScenarios,'LowerBound',DMIN(:,1:T,:),'UpperBound',DMAX(:,1:T,:)); % (TN) demand variable - demand bounds (b.6)
pup = optimvar('pup',nmg,T,TotalNumberOfScenarios,'LowerBound',0); % DN power sold to market
pdn = optimvar('pdn',nmg,T,TotalNumberOfScenarios,'LowerBound',0); % DN power bought from market
theta = optimvar('theta',n-1,T,TotalNumberOfScenarios); % TN voltage angle variable
lamda = optimvar('l',n,T,TotalNumberOfScenarios); % Market prices
fgmin = optimvar('fgmin',ng,T,TotalNumberOfScenarios,'LowerBound',0); % Lower-level dual variable 
fgmax = optimvar('fgmax',ng,T,TotalNumberOfScenarios,'LowerBound',0); % Lower-level dual variable 
fgrd = optimvar('fgrd',ng,T,TotalNumberOfScenarios,'LowerBound',0); % Lower-level dual variable 
fgru = optimvar('fgru',ng,T,TotalNumberOfScenarios,'LowerBound',0); % Lower-level dual variable 
fdmin = optimvar('fdmin',nd,T,TotalNumberOfScenarios,'LowerBound',0); % Lower-level dual variable 
fdmax = optimvar('fdmax',nd,T,TotalNumberOfScenarios,'LowerBound',0); % Lower-level dual variable 
fpupmin = optimvar('fpupmin',nmg,T,TotalNumberOfScenarios,'LowerBound',0); % Lower-level dual variable 
fpdnmin = optimvar('fpdnmin',nmg,T,TotalNumberOfScenarios,'LowerBound',0); % Lower-level dual variable 
fpupmax = optimvar('fpupmax',nmg,T,TotalNumberOfScenarios,'LowerBound',0); % Lower-level dual variable 
fpdnmax = optimvar('fpdnmax',nmg,T,TotalNumberOfScenarios,'LowerBound',0); % Lower-level dual variable 
flmin = optimvar('flmin',nl,T,TotalNumberOfScenarios,'LowerBound',0); % Lower-level dual variable 
flmax = optimvar('flmax',nl,T,TotalNumberOfScenarios,'LowerBound',0); % Lower-level dual variable 
ugmin = optimvar('ugmin',ng,T,TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1); % Big-M auxiliary binary variable (complementarity conditions in lower-level problem KKTs)
ugmax = optimvar('ugmax',ng,T,TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1); % Big-M auxiliary binary variable (complementarity conditions in lower-level problem KKTs)
ugrd = optimvar('ugrd',ng,T,TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1); % Big-M auxiliary binary variable (complementarity conditions in lower-level problem KKTs)
ugru = optimvar('ugru',ng,T,TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1); % Big-M auxiliary binary variable (complementarity conditions in lower-level problem KKTs)
udmin = optimvar('udmin',nd,T,TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1); % Big-M auxiliary binary variable (complementarity conditions in lower-level problem KKTs)
udmax = optimvar('udmax',nd,T,TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1); % Big-M auxiliary binary variable (complementarity conditions in lower-level problem KKTs)
upupmin = optimvar('upupmin',nmg,T,TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1); % Big-M auxiliary binary variable (complementarity conditions in lower-level problem KKTs)
updnmin = optimvar('updnmin',nmg,T,TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1); % Big-M auxiliary binary variable (complementarity conditions in lower-level problem KKTs)
upupmax = optimvar('upupmax',nmg,T,TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1); % Big-M auxiliary binary variable (complementarity conditions in lower-level problem KKTs)
updnmax = optimvar('updnmax',nmg,T,TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1); % Big-M auxiliary binary variable (complementarity conditions in lower-level problem KKTs)
ulmin = optimvar('ulmin',nl,T,TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1); % Big-M auxiliary binary variable (complementarity conditions in lower-level problem KKTs)
ulmax = optimvar('ulmax',nl,T,TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1); % Big-M auxiliary binary variable (complementarity conditions in lower-level problem KKTs)
zpv = optimvar('zpv',sum(npv),T,TotalNumberOfScenarios,'LowerBound',0); % PV production variable
zw = optimvar('zw',sum(nw),T,TotalNumberOfScenarios,'LowerBound',0); % Wind production variable
uline = optimvar('uline',sum(ncbr),'Type','integer','LowerBound',0,'UpperBound',1); % (a.2)
yline = optimvar('yline',sum(nbr),TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1); % (a.27)
bline = optimvar('bline',sum(nn)+nmg,sum(nn)+nmg,TotalNumberOfScenarios,'Type','integer','LowerBound',0,'UpperBound',1); % (a.33)

%% Constraints
% % Total Line Budget Constraint (a.3)
sizing_op.Constraints.total_line_budget = sum((CLinesCostMat.*CLinesLenMat).*uline)<=TotalLineBudget;

% Total Budget Constraint (a.4)
original_bilevel.Constraints.total_budget = sum(TC_Energy*KE+TC_Power*KP)+sum(TC_Wind*KW)+sum(TC_Solar*KPV)<=TotalDERBudget;

% Bounds on Size variables
KW.UpperBound = 0; % (a.5)
KPV.UpperBound = 0; % (a.6)
KE.UpperBound = 0; % (a.7)
KP.UpperBound = 0; % (a.8)
% Storage Power-to-Energy ratio Constraint (a.9)
original_bilevel.Constraints.energy_to_power_ratio = P_CONST*KP==KE;
% RES Production Bounds
wind_upper = optimconstr(sum(nw),T,TotalNumberOfScenarios);
solar_upper = optimconstr(sum(npv),T,TotalNumberOfScenarios);
for ii = 1:nmg
    for tt=1:T
        for ww=1:TotalNumberOfScenarios
            if ii==1
                solar_upper(1:npv(ii),tt,ww) = zpv(1:npv(ii),tt,ww)<=PVEfficiency*SolarTimeSeries(tt,ii,ww).*KPV(1:npv(ii));
                wind_upper(1:nw(ii),tt,ww) = zw(1:nw(ii),tt,ww)<=WindTimeSeries(tt,ii,ww).*KW(1:nw(ii));
            elseif ii>1
                solar_upper(sum(npv(1:ii-1))+1:sum(npv(1:ii)),tt,ww) = zpv(sum(npv(1:ii-1))+1:sum(npv(1:ii)),tt,ww)<=PVEfficiency*SolarTimeSeries(tt,ii,ww).*KPV(sum(npv(1:ii-1))+1:sum(npv(1:ii)));
                wind_upper(sum(nw(1:ii-1))+1:sum(nw(1:ii)),tt,ww) = zw(sum(nw(1:ii-1))+1:sum(nw(1:ii)),tt,ww)<=WindTimeSeries(tt,ii,ww).*KW(sum(nw(1:ii-1))+1:sum(nw(1:ii)));
            end
        end
    end
end
original_bilevel.Constraints.wind_upper = wind_upper; % (a.10)
original_bilevel.Constraints.solar_upper = solar_upper; % (a.11)
% Minimum Return Constraint (a.12)
switch minimum_return_choice
    case 1
        WPROD = optimexpr(nmg,T,TotalNumberOfScenarios);
        PVPROD = optimexpr(nmg,T,TotalNumberOfScenarios);
%         DISCHARGE = optimexpr(nmg,T,TotalNumberOfScenarios);
%         CHARGE = optimexpr(nmg,T,TotalNumberOfScenarios);
        for ii=1:nmg
            for tt=1:T
                for ww=1:TotalNumberOfScenarios
                    if ii==1
                        WPROD(ii,tt,ww) = sum(zw(1:nw(ii),tt,ww));
                        PVPROD(ii,tt,ww) = sum(zpv(1:npv(ii),tt,ww));
%                         DISCHARGE(ii,tt,ww) = sum(dis(1:nbn(ii),tt,ww));
%                         CHARGE(ii,tt,ww) = sum(ch(1:nbn(ii),tt,ww));
                    elseif ii>1
                        WPROD(ii,tt,ww) = sum(zw(sum(nw(1:ii-1))+1:sum(nw(1:ii)),tt,ww));
                        PVPROD(ii,tt,ww) = sum(zpv(sum(npv(1:ii-1))+1:sum(npv(1:ii)),tt,ww));
%                         DISCHARGE(ii,tt,ww) = sum(dis(sum(nbn(1:ii-1))+1:sum(nbn(1:ii)),tt,ww));
%                         CHARGE(ii,tt,ww) = sum(ch(sum(nbn(1:ii-1))+1:sum(nbn(1:ii)),tt,ww));
                    end
                end
            end
        end
        TotalWPROD = sum(sum(WPROD,1),2);
        TotalPVPROD = sum(sum(PVPROD,1),2);
        TotalRESOPEX = CW*TotalWPROD+CPV*TotalPVPROD;
        TotalBSUOPEX = sum(sum(CB*(dis+ch),1),2);
        TotalExpectedOPEX = sum(reshape(Possibilities,1,1,TotalNumberOfScenarios).*(TotalRESOPEX+TotalBSUOPEX));
        
        PPP = -(sum(reshape(Possibilities,1,1,TotalNumberOfScenarios).*((sum(GenBids(:,1,:).*g(:,1,:))-sum(DemBids(:,1,:).*d(:,1,:)))-sum(repmat(GMIN(:,1,:),1,1,TotalNumberOfScenarios).*fgmin(:,1,:))+sum(repmat(GMAX(:,1,:),1,1,TotalNumberOfScenarios).*fgmax(:,1,:))+sum(repmat(RD-G0,1,1,TotalNumberOfScenarios).*fgrd(:,1,:))+sum(repmat(RU+G0,1,1,TotalNumberOfScenarios).*fgru(:,1,:))-sum(fdmin(:,1,:).*DMIN(:,1,:))+sum(fdmax(:,1,:).*DMAX(:,1,:))+sum(repmat(FMAX,1,1,TotalNumberOfScenarios).*flmin(:,1,:))+sum(repmat(FMAX,1,1,TotalNumberOfScenarios).*flmax(:,1,:))))+sum(reshape(Possibilities,1,1,TotalNumberOfScenarios).*((sum(sum(GenBids(:,2:T,:).*g(:,2:T,:)))-sum(sum(DemBids(:,2:T,:).*d(:,2:T,:))))-sum(sum(repmat(GMIN(:,2:T,:),1,1,TotalNumberOfScenarios).*fgmin(:,2:T,:)))+sum(sum(repmat(GMAX(:,2:T,:),1,1,TotalNumberOfScenarios).*fgmax(:,2:T,:)))+sum(sum(repmat(RD,1,T-1,TotalNumberOfScenarios).*fgrd(:,2:T,:)))+sum(sum(repmat(RU,1,T-1,TotalNumberOfScenarios).*fgru(:,2:T,:)))-sum(sum(fdmin(:,2:T,:).*DMIN(:,2:T,:)))+sum(sum(fdmax(:,2:T,:).*DMAX(:,2:T,:)))+sum(sum(repmat(FMAX,1,T-1,TotalNumberOfScenarios).*flmin(:,2:T,:)))+sum(sum(repmat(FMAX,1,T-1,TotalNumberOfScenarios).*flmax(:,2:T,:))))));
        Bill_per_scen = optimexpr(nmg,TotalNumberOfScenarios);
        for gg=1:nmg
            if gg==1
                for ss=1:TotalNumberOfScenarios
                    Bill_per_scen(gg,ss)=Possibilities(ss)*(sum(MGLoc(:,gg)'*lamda(:,:,ss).*sum(ConnectedLoadMat(1:nn(gg),:,ss),1)));
                end
            else
                for ss=1:TotalNumberOfScenarios
                    Bill_per_scen(gg,ss)=Possibilities(ss)*(sum(MGLoc(:,gg)'*lamda(:,:,ss).*sum(ConnectedLoadMat(sum(nn(1:gg-1))+1:sum(nn(1:gg)),:,ss),1)));
                end
            end
        end
        if DERSOPEX==1
            original_bilevel.Constraints.minimum_return = (-TotalExpectedOPEX+PPP+sum(Bill_per_scen(:)))>=return_ratio*(sum(AC_Energy*KE+AC_Power*KP)+sum(AC_Wind*KW)+sum(AC_Solar*KPV));
        elseif DERSOPEX==0
            original_bilevel.Constraints.minimum_return = (PPP+sum(Bill_per_scen(:)))>=return_ratio*(sum(AC_Energy*KE+AC_Power*KP)+sum(AC_Wind*KW)+sum(AC_Solar*KPV));
        end
    case 0
        WPROD = optimexpr(nmg,T,TotalNumberOfScenarios);
        PVPROD = optimexpr(nmg,T,TotalNumberOfScenarios);
%         DISCHARGE = optimexpr(nmg,T,TotalNumberOfScenarios);
%         CHARGE = optimexpr(nmg,T,TotalNumberOfScenarios);
        for ii=1:nmg
            for tt=1:T
                for ww=1:TotalNumberOfScenarios
                    if ii==1
                        WPROD(ii,tt,ww) = sum(zw(1:nw(ii),tt,ww));
                        PVPROD(ii,tt,ww) = sum(zpv(1:npv(ii),tt,ww));
%                         DISCHARGE(ii,tt,ww) = sum(dis(1:nbn(ii),tt,ww));
%                         CHARGE(ii,tt,ww) = sum(ch(1:nbn(ii),tt,ww));
                    elseif ii>1
                        WPROD(ii,tt,ww) = sum(zw(sum(nw(1:ii-1))+1:sum(nw(1:ii)),tt,ww));
                        PVPROD(ii,tt,ww) = sum(zpv(sum(npv(1:ii-1))+1:sum(npv(1:ii)),tt,ww));
%                         DISCHARGE(ii,tt,ww) = sum(dis(sum(nbn(1:ii-1))+1:sum(nbn(1:ii)),tt,ww));
%                         CHARGE(ii,tt,ww) = sum(ch(sum(nbn(1:ii-1))+1:sum(nbn(1:ii)),tt,ww));
                    end
                end
            end
        end
        TotalWPROD = sum(sum(WPROD,1),2);
        TotalPVPROD = sum(sum(PVPROD,1),2);
        TotalRESOPEX = CW*TotalWPROD+CPV*TotalPVPROD;
        TotalBSUOPEX = sum(sum(CB*(dis+ch),1),2);
        TotalExpectedOPEX = sum(reshape(Possibilities,1,1,TotalNumberOfScenarios).*(TotalRESOPEX+TotalBSUOPEX));
               
    otherwise
        error('Wrong choice for the incorporation of a minimum return on annualized investment cost!!')
end

% Offer/Bid Binary Upper Bound (a.13)-(a.14)
original_bilevel.Constraints.offer_binary_bound = o<=m.*repmat(SubPuMat,1,T,TotalNumberOfScenarios);
original_bilevel.Constraints.bid_binary_bound = b<=(1-m).*repmat(SubPuMat,1,T,TotalNumberOfScenarios);

% Power Discharge/Charge Upper Bound (a.16)
original_bilevel.Constraints.dis_upper = dis<=repmat(KP,1,T,TotalNumberOfScenarios);
original_bilevel.Constraints.ch_upper = ch<=repmat(KP,1,T,TotalNumberOfScenarios);

% Storage Binary Variable Constraints
original_bilevel.Constraints.dis_x = dis<=BIGX.*x; % (a.17)
original_bilevel.Constraints.ch_x = ch<=BIGX.*(1-x); % (a.18)

% State of Energy Dynamics (a.20 - a.21)
original_bilevel.Constraints.soc_dynamicst2 = soc(:,2:T,:) == soc(:,1:T-1,:)-dis(:,2:T,:)./repmat(DischargeEff,1,T-1,TotalNumberOfScenarios)+ch(:,2:T,:).*repmat(ChargeEff,1,T-1,TotalNumberOfScenarios);
original_bilevel.Constraints.soc_dynamicst1 = soc(:,1,:) == repmat(InitialSOCPerMat,1,1,TotalNumberOfScenarios).*repmat(KE,1,1,TotalNumberOfScenarios)-dis(:,1,:)./repmat(DischargeEff,1,1,TotalNumberOfScenarios)+ch(:,1,:).*repmat(ChargeEff,1,1,TotalNumberOfScenarios);
% State of Energy Upper Bound (a.22)
original_bilevel.Constraints.soc_upper_bound = soc<=repmat(KE,1,T,TotalNumberOfScenarios);
% State of Energy Last Timeslot (a.23)
original_bilevel.Constraints.soc_last_timeslot = soc(:,T,:)>=repmat(FinalSOCPerMat,1,1,TotalNumberOfScenarios).*repmat(KE,1,1,TotalNumberOfScenarios);

% DN Active and Reactive Power Balance (a.24), (a.25)
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
                active_power_flow(1:nn(ii),tt,ww) = CoefActivePF{ii}*ppf(1:nbr(ii),tt,ww) - ResLoc{ii}*zw(1:nw(ii),tt,ww) - PVLoc{ii}*zpv(1:npv(ii),tt,ww)- BatLoc{ii}*dis(1:nbn(ii),tt,ww) + BatLoc{ii}*ch(1:nbn(ii),tt,ww) == -ConnectedLoadMat(1:nn(ii),tt,ww);
                reactive_power_flow(1:nn(ii),tt,ww) = CoefReactivePF{ii}*qpf(1:nbr(ii),tt,ww) - (ResLoc{ii}*zw(1:nw(ii),tt,ww).*WindPF(1)) - ((PVLoc{ii}*zpv(1:npv(ii),tt,ww)).*SolarPF(1)) == -aInfMat(1:nn(ii),tt).*(ConnectedLoadMat(1:nn(ii),tt,ww));
            elseif ii>1
                active_power_flow(sum(nn(1:ii-1))+1:sum(nn(1:ii)),tt,ww) = CoefActivePF{ii}*ppf(sum(nbr(1:ii-1))+1:sum(nbr(1:ii)),tt,ww) - ResLoc{ii}*zw(sum(nw(1:ii-1))+1:sum(nw(1:ii)),tt,ww) - PVLoc{ii}*zpv(sum(npv(1:ii-1))+1:sum(npv(1:ii)),tt,ww) - BatLoc{ii}*dis(sum(nbn(1:ii-1))+1:sum(nbn(1:ii)),tt,ww) + BatLoc{ii}*ch(sum(nbn(1:ii-1))+1:sum(nbn(1:ii)),tt,ww) == -ConnectedLoadMat(sum(nn(1:ii-1))+1:sum(nn(1:ii)),tt,ww);
                reactive_power_flow(sum(nn(1:ii-1))+1:sum(nn(1:ii)),tt,ww) = CoefReactivePF{ii}*qpf(sum(nbr(1:ii-1))+1:sum(nbr(1:ii)),tt,ww) - (ResLoc{ii}*(zw(sum(nw(1:ii-1))+1:sum(nw(1:ii)),tt,ww))).*WindPF(1) - (PVLoc{ii}*(zpv(sum(npv(1:ii-1))+1:sum(npv(1:ii)),tt,ww))).*SolarPF(1)  == -aInfMat(sum(nn(1:ii-1))+1:sum(nn(1:ii)),tt).*(ConnectedLoadMat(sum(nn(1:ii-1))+1:sum(nn(1:ii)),tt,ww));
            end
        end
    end
end

original_bilevel.Constraints.active_power_flow_power_balance = active_power_flow;
original_bilevel.Constraints.reactive_power_flow_power_balance = reactive_power_flow;

% Voltage-Branch Flow (a.26)
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
original_bilevel.Constraints.voltage_branch_upper = voltage_branch_upper;
original_bilevel.Constraints.voltage_branch_lower = voltage_branch_lower;


% Constraints (a.28 - a.29) that forbid the usage of a non-installed line
a28 = optimconstr(sum(nebr),TotalNumberOfScenarios);
a32 = optimconstr(sum(ncbr),TotalNumberOfScenarios);
for ii=1:nmg
    for ss=1:TotalNumberOfScenarios
        if ii==1
            a28(1:nebr(ii),ss) = yline(1:nebr(ii),ss)<=1;
            a32(1:ncbr(ii),ss) = yline(nebr(ii)+1:nbr(ii),ss)<=uline(1:ncbr(ii));
        elseif ii>1
            a28(sum(nebr(1:ii-1))+1:sum(nebr(1:ii)),ss) = yline(sum(nbr(1:ii-1))+1:sum(nbr(1:ii-1))+sum(nebr(ii)),ss)<=1;
            a32(sum(ncbr(1:ii-1))+1:sum(ncbr(1:ii)),ss) = yline(sum(nbr(1:ii-1))+sum(nebr(ii))+1:sum(nbr(1:ii)),ss)<=uline(sum(ncbr(1:ii-1))+1:sum(ncbr(1:ii)));
        end
    end
end
original_bilevel.Constraints.a28 = a28;% Existing lines can always be utilized (is it necessary????? I think NO!)
original_bilevel.Constraints.a29 = a32;

% Constraint (a.30)
a30 = optimconstr(sum(nbr),TotalNumberOfScenarios);
for ss=1:TotalNumberOfScenarios
    %     for ii=1:nmg
    for ll=1:sum(nbr)
        a30(ll,ss) = yline(ll,ss) == bline(FromNodesMat(ll),ToNodesMat(ll),ss)+bline(ToNodesMat(ll),FromNodesMat(ll),ss);
    end
    %     end
end
original_bilevel.Constraints.a30 = a30;
%%%%%%% SHOULD A CONSTRAINT THAT NULLIFIES ? NO-LINE ELEMENTS??? %%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Constraint (a.32)
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
a32 = optimconstr(sum(nn),TotalNumberOfScenarios);
for ii=1:sum(nn)+nmg
    if ~ismember(ii,SourceNodes)
        BCoef(FromNodesMat(LinesMat(ToNodesMat==ii)),ii) = 1;
    end
end
AllBs = sum(bline(:,:,1).*BCoef,1);
AllBs(SourceNodes) = [];
for ss=1:TotalNumberOfScenarios
    a32(:,ss) = AllBs<=ones(1,sum(nn)); % == OR <= ?????
end
original_bilevel.Constraints.a32 = a32;

% Line Flow - Apparent Power Upper Bound (a.34)
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
        original_bilevel.Constraints.ppf_upper = ppf_upper;
        original_bilevel.Constraints.ppf_lower = ppf_lower;
        original_bilevel.Constraints.qpf_upper = qpf_upper;
        original_bilevel.Constraints.qpf_lower = qpf_lower;
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
        original_bilevel.Constraints.dn_flows_constraint = apparent_bound;
    otherwise
        error('There is no such an option for line capacity constraint!');
end

% Constraint (a.36) - SEEMS REDUNDANT IF CONSTRAINT (a.32) IS AN EQUALITY CONSTRAINT!!!!!!!!!!
% %%%%%%%%%%%%%%%%%%%%% NEEDS FIXING!!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % nn(ii) is wrong. The number of the connected nodes should be there and
% % not the number of the total existing and candidate nodes!!!!!!!!!
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a36 = optimconstr(nmg,TotalNumberOfScenarios);
% for ss=1:TotalNumberOfScenarios
%     for ii=1:nmg
%         if ii==1
%             a36(ii,ss) = sum(yline(1:nbr(ii),ss))==nn(ii);
%         elseif ii>1
%             a36(ii,ss) = sum(yline(sum(nbr(1:ii-1))+1:sum(nbr(1:ii)),ss))==nn(ii);
%         end
%     end
% end
% % % % % % sizing_op.Constraints.a36 = a36;

% Active Power Balance at substation (a.37) and Reactive Power Balance at substation
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
original_bilevel.Constraints.substation_active_power_balance = substation_ppb;
original_bilevel.Constraints.substation_reactive_power_balance = substation_qpb;

%%
%%%%%%%%%%%%%%%% Lower Level KKT Conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Primal Constraints
% Power Balance (b.2)/(c.1)
TSOPowerBalance = optimconstr(n,T,TotalNumberOfScenarios);
for tt=1:T
    for ww=1:TotalNumberOfScenarios
        TSOPowerBalance(:,tt,ww) = -GenLoc*g(:,tt,ww)+DemLoc*d(:,tt,ww)-MGLoc*(pup(:,tt,ww)-pdn(:,tt,ww))+B*theta(:,tt,ww)==0;
    end
end
original_bilevel.Constraints.system_power_balance = TSOPowerBalance;
% Ramp Up/Down Constraints (b.4 - b.5)
original_bilevel.Constraints.ramp_up_t1 = g(:,1,:)<=repmat(G0+RU,1,1,TotalNumberOfScenarios);
original_bilevel.Constraints.ramp_up_t2 = g(:,2:T,:)-g(:,1:T-1,:)<=repmat(RU,1,T-1,TotalNumberOfScenarios);

original_bilevel.Constraints.ramp_dn_t1 = -g(:,1,:)<=repmat(-G0+RD,1,1,TotalNumberOfScenarios);
original_bilevel.Constraints.ramp_dn_t2 = g(:,1:T-1,:)-g(:,2:T,:)<=repmat(RD,1,T-1,TotalNumberOfScenarios);
% DN Dispatch Lower/Upper Bounds 
original_bilevel.Constraints.upper_bound_on_ESP = pup<= o; % (b.7)
original_bilevel.Constraints.lower_bound_on_ESP = pdn<= b; % (b.8)
% Line Flows Bounds (b.9)
LineFlowsUpperBound = optimconstr(nl,T,TotalNumberOfScenarios);
LineFlowsLowerBound = optimconstr(nl,T,TotalNumberOfScenarios);
for tt=1:T
    for ww=1:TotalNumberOfScenarios
        LineFlowsUpperBound(:,tt,ww) = Yline*theta(:,tt,ww)<=FMAX;
        LineFlowsLowerBound(:,tt,ww) = -Yline*theta(:,tt,ww)<=FMAX;
    end
end

original_bilevel.Constraints.upper_bound_on_line_flows = LineFlowsUpperBound;
original_bilevel.Constraints.lower_bound_on_line_flows = LineFlowsLowerBound;
%%%%%%%%%%%%%%%%% Stationarity Conditions
% Partial Derivative of Lagrangian w.r.t. g
gen_derivative_T = optimconstr(ng,1,TotalNumberOfScenarios);
gen_derivative = optimconstr(ng,T-1,TotalNumberOfScenarios);
for ww = 1:TotalNumberOfScenarios
    for tt=1:T-1
        gen_derivative(:,tt,ww) = GenBids(:,tt,ww) - GenLoc'*lamda(:,tt,ww)-fgmin(:,tt,ww)+fgmax(:,tt,ww)-fgrd(:,tt,ww)+fgrd(:,tt+1,ww)+fgru(:,tt,ww)-fgru(:,tt+1,ww)==0;
    end
    gen_derivative_T(:,:,ww) = GenBids(:,T,ww) - GenLoc'*lamda(:,T,ww)-fgmin(:,T,ww)+fgmax(:,T,ww)-fgrd(:,T,ww)+fgru(:,T,ww)==0;
end

original_bilevel.Constraints.gen_derivative_T = gen_derivative_T; % (c.3)
original_bilevel.Constraints.gen_derivative = gen_derivative; % (c.2)

% Partial Derivatives of Lagrangian w.r.t. d, pup, pdn, theta
demand_derivative = optimconstr(nd,T,TotalNumberOfScenarios);
pup_derivative = optimconstr(nmg,T,TotalNumberOfScenarios);
pdn_derivative = optimconstr(nmg,T,TotalNumberOfScenarios);
theta_derivative = optimconstr(n-1,T,TotalNumberOfScenarios);

for tt=1:T
    for ww=1:TotalNumberOfScenarios
        demand_derivative(:,tt,ww) = -DemBids(:,tt,ww)+ DemLoc'*lamda(:,tt,ww) - fdmin(:,tt,ww)+fdmax(:,tt,ww) ==0;
        pup_derivative(:,tt,ww) = cpup(:,tt,ww) - MGLoc'*lamda(:,tt,ww)-fpupmin(:,tt,ww)+fpupmax(:,tt,ww) == 0;
        pdn_derivative(:,tt,ww) = -cpdn(:,tt,ww) + MGLoc'*lamda(:,tt,ww)-fpdnmin(:,tt,ww)+fpdnmax(:,tt,ww) == 0;
        theta_derivative(:,tt,ww) = B'*lamda(:,tt,ww)-Yline'*flmin(:,tt,ww)+Yline'*flmax(:,tt,ww) == 0;
    end
end

original_bilevel.Constraints.demand_derivative = demand_derivative; % (c.4)
original_bilevel.Constraints.pup_derivative = pup_derivative; % (c.5)
original_bilevel.Constraints.pdn_derivative = pdn_derivative; % (c.6)
original_bilevel.Constraints.theta_derivative = theta_derivative; % (c.7)
%%%%%%%%%%%%% Complementarity Slackness
original_bilevel.Constraints.comp_slack_dual1 = fgmin<=M*(1-ugmin); % (c.8)
original_bilevel.Constraints.comp_slack_dual2 = fgmax<=M*(1-ugmax); % (c.9)
original_bilevel.Constraints.comp_slack_dual3 = fgrd<=M*(1-ugrd); % (c.10)
original_bilevel.Constraints.comp_slack_dual4 = fgru<=M*(1-ugru); % (c.11)
original_bilevel.Constraints.comp_slack_dual5 = fdmin<=M*(1-udmin); % (c.12)
original_bilevel.Constraints.comp_slack_dual6 = fdmax<=M*(1-udmax); % (c.13)
original_bilevel.Constraints.comp_slack_dual7 = fpupmin<=M*(1-upupmin); % (c.14)
original_bilevel.Constraints.comp_slack_dual8 = fpupmax<=M*(1-upupmax); % (c.15)
original_bilevel.Constraints.comp_slack_dual9 = fpdnmin<=M*(1-updnmin); % (c.16)
original_bilevel.Constraints.comp_slack_dual10 = fpdnmax<=M*(1-updnmax); % (c.17)
original_bilevel.Constraints.comp_slack_dual11 = flmin<=M*(1-ulmin); % (c.18)
original_bilevel.Constraints.comp_slack_dual12 = flmax<=M*(1-ulmax); % (c.19)

original_bilevel.Constraints.comp_slack_primal1 = g-repmat(GMIN,1,1,TotalNumberOfScenarios)<=M*ugmin; % (c.8)
original_bilevel.Constraints.comp_slack_primal2 = -g+repmat(GMAX,1,1,TotalNumberOfScenarios)<=M*ugmax; % (c.9)
original_bilevel.Constraints.comp_slack_primal3_t1 = g(:,1,:)-repmat(G0,1,1,TotalNumberOfScenarios)+repmat(RD,1,1,TotalNumberOfScenarios)<=M*ugrd(:,1,:); % (c.10)
original_bilevel.Constraints.comp_slack_primal3_t2 = g(:,2:T,:)-g(:,1:T-1,:)+repmat(RD,1,T-1,TotalNumberOfScenarios)<=M*ugrd(:,2:T,:); % (c.10)
original_bilevel.Constraints.comp_slack_primal4_t1 = -g(:,1,:)+repmat(G0,1,1,TotalNumberOfScenarios)+repmat(RU,1,1,TotalNumberOfScenarios)<=M*ugru(:,1,:); % (c.11)
original_bilevel.Constraints.comp_slack_primal4_t2 = -g(:,2:T,:)+g(:,1:T-1,:)+repmat(RU,1,T-1,TotalNumberOfScenarios)<=M*ugru(:,2:T,:); % (c.11)
original_bilevel.Constraints.comp_slack_primal5 = d-DMIN(:,1:T,:)<=M*udmin; % (c.12)
original_bilevel.Constraints.comp_slack_primal6 = -d+DMAX(:,1:T,:)<=M*udmax; % (c.13)
original_bilevel.Constraints.comp_slack_primal10 = pup<=M*upupmin; % (c.14)
original_bilevel.Constraints.comp_slack_primal9 = -pup+o<=M*upupmax; % (c.15)
original_bilevel.Constraints.comp_slack_primal8 = pdn<=M*updnmin; % (c.16)
original_bilevel.Constraints.comp_slack_primal7 = -pdn+b<=M*updnmax; % (c.17)



comp_slack_primal11 = optimconstr(nl,T,TotalNumberOfScenarios);
comp_slack_primal12 = optimconstr(nl,T,TotalNumberOfScenarios);
for ww=1:TotalNumberOfScenarios
    for tt=1:T
        comp_slack_primal11(:,tt,ww) = Yline*theta(:,tt,ww)+FMAX<=M*ulmin(:,tt,ww);
        comp_slack_primal12(:,tt,ww) = -Yline*theta(:,tt,ww)+FMAX<=M*ulmax(:,tt,ww);
    end
end

original_bilevel.Constraints.comp_slack_primal11 = comp_slack_primal11; % (c.18)
original_bilevel.Constraints.comp_slack_primal12 = comp_slack_primal12; % (c.19)

% Guarantee No Lines Investments
original_bilevel.Constraints.fix_u = uline==0;

%% Objective
DERS_CAPEX = sum(AC_Wind*KW)+sum(AC_Solar*KPV)+sum(AC_Energy*KE)+sum(AC_Power*KP);
LINES_CAPEX = sum((AC_Lines.*CLinesLenMat).*uline);
UsersBills = optimexpr(nmg,TotalNumberOfScenarios);
HypotheticUserBills = optimexpr(nmg,TotalNumberOfScenarios);
for ss=1:TotalNumberOfScenarios
    Prices_s = MGLoc'*lamda(:,:,ss);
    for ii=1:nmg
        nodal_prices = Prices_s(ii,:);
        if ii==1
            UsersBills(ii,ss) = Possibilities(ss)*sum(nodal_prices.*sum(ConnectedLoadMat(1:nn(ii),:,ss),1));
            HypotheticUserBills(ii,ss) = Possibilities(ss)*sum(nodal_prices.*sum(InfLoadMat(1:nn(ii),:,ss),1));
        elseif ii>1
            UsersBills(ii,ss) = Possibilities(ss)*sum(nodal_prices.*sum(ConnectedLoadMat(sum(nn(1:ii-1))+1:sum(nn(1:ii)),:,ss),1)); 
            HypotheticUserBills(ii,ss) = Possibilities(ss)*sum(nodal_prices.*sum(InfLoadMat(sum(nn(1:ii-1))+1:sum(nn(1:ii)),:,ss),1)); 
        end
    end
end
TotalBills = sum(sum(UsersBills));
TotalHypotheticBills=  sum(sum(HypotheticUserBills));

switch ModelObj
    case 1
        Operating_Objective = TotalBills;
    case 2
        Operating_Objective = sum(reshape(Possibilities,1,1,TotalNumberOfScenarios).*((sum(GenBids(:,1,:).*g(:,1,:))-sum(DemBids(:,1,:).*d(:,1,:)))-sum(fgmin(:,1,:).*repmat(GMIN(:,1,:),1,1,TotalNumberOfScenarios))+sum(fgmax(:,1,:).*repmat(GMAX(:,1,:),1,1,TotalNumberOfScenarios))+sum(repmat(RD-G0,1,1,TotalNumberOfScenarios).*fgrd(:,1,:))+sum(repmat(RU+G0,1,1,TotalNumberOfScenarios).*fgru(:,1,:))-sum(fdmin(:,1,:).*DMIN(:,1,:))+sum(fdmax(:,1,:).*DMAX(:,1,:))+sum(repmat(FMAX,1,1,TotalNumberOfScenarios).*flmin(:,1,:))+sum(repmat(FMAX,1,1,TotalNumberOfScenarios).*flmax(:,1,:))))+sum(reshape(Possibilities,1,1,TotalNumberOfScenarios).*(sum(sum(GenBids(:,2:T,:).*g(:,2:T,:))) - sum(sum(DemBids(:,2:T,:).*d(:,2:T,:))) - sum(sum(fgmin(:,2:T,:).*repmat(GMIN(:,2:T,:),1,1,TotalNumberOfScenarios))) + sum(sum(fgmax(:,2:T,:).*repmat(GMAX(:,2:T,:),1,1,TotalNumberOfScenarios)))++sum(sum(repmat(RD,1,T - 1,TotalNumberOfScenarios).*fgrd(:,2:T,:))) + sum(sum(repmat(RU,1,T - 1,TotalNumberOfScenarios).*fgru(:,2:T,:))) - sum(sum(fdmin(:,2:T,:).*DMIN(:,2:T,:))) + sum(sum(fdmax(:,2:T,:).*DMAX(:,2:T,:))) + sum(sum(repmat(FMAX,1,T - 1,TotalNumberOfScenarios).*flmin(:,2:T,:))) + sum(sum(repmat(FMAX,1,T - 1,TotalNumberOfScenarios).*flmax(:,2:T,:)))));
    otherwise
        error('Wrong choice for the Operating Objective!!')
end

if DERSOPEX==1
    OPEX = TotalExpectedOPEX + Operating_Objective;
elseif DERSOPEX==0
    OPEX = Operating_Objective;                    
end

original_bilevel.Objective = DERS_CAPEX+LINES_CAPEX+OPEX;
%% Solve Model
options = optimoptions('intlinprog','Display','iter');
[init_x_opt,neg_profits,exitflag,output]=solve(original_bilevel,'solver','intlinprog');

%% Original Results
TSOSocialWelfare_pu = zeros(1,TotalNumberOfScenarios);
GenerationCost_pu = zeros(1,TotalNumberOfScenarios);
ExpectedVOLL_pu = zeros(1,TotalNumberOfScenarios);
for ss=1:TotalNumberOfScenarios
    TSOSocialWelfare_pu(ss) = -(sum(sum(GenBids(:,:,ss).*init_x_opt.g(:,:,ss)))-sum(sum(DemBids(:,:,ss).*init_x_opt.d(:,:,ss)))+sum(sum(cpup(:,:,ss).*init_x_opt.pup(:,:,ss)))-sum(sum(cpdn(:,:,ss).*init_x_opt.pdn(:,:,ss))));
    GenerationCost_pu(ss) = sum(sum(GenBids(:,:,ss).*init_x_opt.g(:,:,ss)))+sum(sum(cpup(:,:,ss).*init_x_opt.pup(:,:,ss)));
    ExpectedVOLL_pu(ss) = sum(sum(IslandedLoadMat(:,:,ss)*VoLL));
end
ExpectedGenerationCost_pu = Possibilities*GenerationCost_pu';
ExpectedTSOSocialWelfare_pu = Possibilities*TSOSocialWelfare_pu';
ExpectedVOLL_pu = Possibilities*ExpectedVOLL_pu';

Original_LMPs = zeros(nmg,T,TotalNumberOfScenarios);
for ss=1:TotalNumberOfScenarios
    Original_LMPs(:,:,ss) = MGLoc'*init_x_opt.l(:,:,ss);
end

Original_Bills = evaluate(TotalBills,init_x_opt);
Original_Hypo_Bills = evaluate(TotalHypotheticBills,init_x_opt);
GC0 = ExpectedGenerationCost_pu;
SW0 = ExpectedTSOSocialWelfare_pu;
VOLL0 = ExpectedVOLL_pu;
VOLL0(VOLL0<0)=0;
Original_OPEX = Original_Bills+VOLL0;