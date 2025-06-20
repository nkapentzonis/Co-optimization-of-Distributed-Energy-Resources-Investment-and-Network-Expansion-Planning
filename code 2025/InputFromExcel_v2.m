TestSystemXLFile = [TNDataPath,TestSystemList{TestSystemSelection}];
MGTestSystemXLSheet = MGTestSystemList{MGTestSystemSelection};

%% Grid Structure
% Number of Buses
n = xlsread([TestSystemXLFile '.xlsx'],'Sheet1','A2');

% Number of Lines
nl = xlsread([TestSystemXLFile '.xlsx'],'Sheet1','B2');

% Number of Generators
ng = xlsread([TestSystemXLFile '.xlsx'],'Sheet1','C2');

% Number of Loads
nd = xlsread([TestSystemXLFile '.xlsx'],'Sheet1','D2');

% Number of MGs
nmg = xlsread([TestSystemXLFile '.xlsx'],'Sheet1','E2');


% GenLoc defines the position (bus index) of each generator 
LastRow = sprintf('%d',2+ng-1);
GenBus = xlsread([TestSystemXLFile '.xlsx'],'Sheet3',strcat('B2:B',LastRow));
GenLoc = zeros(n,ng);
for gg = 1:ng
    GenLoc(GenBus(gg),gg) = 1;
end

% DemLoc defines the position (bus index) of each FSP's load
LastRow = sprintf('%d',2+nd-1);
DemBus = xlsread([TestSystemXLFile '.xlsx'],'Sheet4',strcat('A2:A',LastRow));
DemLoc = zeros(n,nd);
for dd=1:nd
    DemLoc(DemBus(dd),dd) = 1;
end

% MGLoc defines the position (bus index) of each DN
LastRow = sprintf('%d',2+nmg-1);
MGBus = xlsread([TestSystemXLFile '.xlsx'],'Sheet1',strcat('F2:F',LastRow));
MGLoc = zeros(n,nmg);
for mg=1:nmg
    MGLoc(MGBus(mg),mg) = 1;
end

% Power Lines Susceptance
LastRow = sprintf('%d',2+nl-1);
LinesSusc = transpose(xlsread([TestSystemXLFile '.xlsx'],'Sheet2',strcat('E2:E',LastRow)));

% Source Buses
FromBus = xlsread([TestSystemXLFile '.xlsx'],'Sheet2',strcat('B2:B',LastRow));
% Destination Buses
ToBus = xlsread([TestSystemXLFile '.xlsx'],'Sheet2',strcat('C2:C',LastRow));

% Susceptance Matrix and LineFlows Matrix
% LineFlows indicates the starting and ending nodes of each line
B = zeros(n);
LineFlows = zeros(nl,n);
for kk=1:nl
    % Off Diagonal elements of B Matrix
    B(FromBus(kk),ToBus(kk)) = -LinesSusc(kk);
    B(ToBus(kk),FromBus(kk)) = B(FromBus(kk),ToBus(kk));
    % LineFlows Matrix
    LineFlows(kk,FromBus(kk)) = 1;
    LineFlows(kk,ToBus(kk)) = -1;
end
% Diagonal Elements of B Matrix
for ii=1:n
    B(ii,ii) = -sum(B(ii,:));
end

% Power Lines' Capacities
LastRow = sprintf('%d',2+nl-1);
FMAX = xlsread([TestSystemXLFile '.xlsx'],'Sheet2',strcat('F2:F',LastRow))./MVA;

%% Thermal Generators
LastRow = sprintf('%d',2+ng-1);
% Generators' Price Offers
GenBidsPerTimeslot = xlsread([TestSystemXLFile '.xlsx'],'Sheet3',strcat('H2:','H',LastRow));
GenBids = zeros(ng,T,TotalNumberOfScenarios);
for ss=1:TotalNumberOfScenarios
    if RepDays==1
        GenBids(:,:,ss) = GenOffersScenarios(floor((ss-1)/NumberOfScenarios)+1)*GenBidsPerTimeslot*ones(1,T);
    elseif RepDays==0
        GenBids(:,:,ss) = GenOffersScenarios(floor((ss-1)/(NumberOfLoadScenarios*NumberOfSolarScenarios*NumberOfWindScenarios))+1)*GenBidsPerTimeslot*ones(1,T);
    end
end

% Generators' Production Lower Bounds
GMINPerTimeslot = xlsread([TestSystemXLFile '.xlsx'],'Sheet3',strcat('D2:D',LastRow))./MVA;
GMIN = GMINPerTimeslot*ones(1,T);

% Generators' Production Upper Bounds
GMAXPerTimeslot = xlsread([TestSystemXLFile '.xlsx'],'Sheet3',strcat('C2:C',LastRow))./MVA;
GMAX = GMAXPerTimeslot*ones(1,T);

% Generators' Ramping Up Limits
RU = xlsread([TestSystemXLFile '.xlsx'],'Sheet3',strcat('E2:E',LastRow))./MVA;

% Generators' Ramping Down Limits
RD = xlsread([TestSystemXLFile '.xlsx'],'Sheet3',strcat('F2:F',LastRow))./MVA;

% Generators' Initial Production Level
G0 = xlsread([TestSystemXLFile '.xlsx'],'Sheet3',strcat('G2:G',LastRow))./MVA;

%% Demand
LastRow = sprintf('%d',2+nd-1);
DemBidsPerTimeslot = xlsread([TestSystemXLFile '.xlsx'],'Sheet4',strcat('B2:B',LastRow));
DemBids = repmat(DemBidsPerTimeslot*ones(1,T),1,1,TotalNumberOfScenarios);
PeakDemand = xlsread([TestSystemXLFile '.xlsx'],'Sheet4',strcat('C2:C',LastRow))./MVA;
ValleyDemand = xlsread([TestSystemXLFile '.xlsx'],'Sheet4',strcat('D2:D',LastRow))./MVA;
DMAX = zeros(nd,T,TotalNumberOfScenarios);
DMIN = repmat(ValleyDemand,1,T,TotalNumberOfScenarios);
for dd=1:nd
    for ss=1:TotalNumberOfScenarios
        DMAX(dd,:,ss) = TotalSizingLoadScenarios(ss,:).*PeakDemand(dd);
    end
end

%% Distribution Networks
% DN Structure Data
MaxRow = 100;
MGData = cell(1,nmg);
for mm = 1:nmg
    MGData{mm} = xlsread([DNDataPath,'\TEPSizingInputData'],MGTestSystemXLSheet,xlRC2A1(3,(mm-1)*26+1,MaxRow+2,mm*26));
end

nen = zeros(nmg,1);
ncn = zeros(nmg,1);
nn = zeros(nmg,1);
nebr = zeros(nmg,1);
ncbr = zeros(nmg,1);
nbr = zeros(nmg,1);
nfl = zeros(nmg,1);
ninfl = zeros(nmg,1);
nbn = zeros(nmg,1);
nw = zeros(nmg,1);
npv = zeros(nmg,1);
ExFromBusMG = cell(1,nmg);
ExToBusMG = cell(1,nmg);
CFromBusMG = cell(1,nmg);
CToBusMG = cell(1,nmg);
FromBusMG = cell(1,nmg);
ToBusMG = cell(1,nmg);
ExLines = cell(1,nmg);
CLines = cell(1,nmg);
Lines = cell(1,nmg);
ExLinesR = cell(1,nmg);
ExLinesX = cell(1,nmg);
CLinesR = cell(1,nmg);
CLinesX = cell(1,nmg);
LinesR = cell(1,nmg);
LinesX = cell(1,nmg);
CLinesLen = cell(1,nmg);
CLinesCost = cell(1,nmg);
ExLinesPu = cell(1,nmg);
ExLinesPl = cell(1,nmg);
ExLinesQu = cell(1,nmg);
ExLinesQl = cell(1,nmg);
CLinesPu = cell(1,nmg);
CLinesPl = cell(1,nmg);
CLinesQu = cell(1,nmg);
CLinesQl = cell(1,nmg);
NVMu = cell(1,nmg);
NVMl = cell(1,nmg);
initialSOC = cell(1,nmg);
finalSOC = cell(1,nmg);
InfLoc = cell(1,nmg);
InfLoad = cell(1,nmg);
TotalInfLoadPerMG = cell(1,nmg);
alphaInflex = cell(1,nmg);
alphaRes = cell(1,nmg);
alphaPV = cell(1,nmg);
SubPu = cell(1,nmg);
SubPl = cell(1,nmg);
WindInt = cell(1,nmg);
PVInt = cell(1,nmg);
BatLoc = cell(1,nmg);
ResLoc = cell(1,nmg);
PVLoc = cell(1,nmg);
for ii=1:nmg
    mgd = MGData{ii};  
    fbmg = mgd(:,2);tbmg = mgd(:,3);
    fbmg(isnan(fbmg)) = [];tbmg(isnan(tbmg)) = [];
    ExFromBusMG{ii} = fbmg;ExToBusMG{ii} = tbmg;
    mgdd = mgd(:,1);
    ExLines{ii} = mgdd(mgdd>0);    
    
    % Number of Existing Nodes per DN
    nen(ii) = max(max(ExFromBusMG{ii}),max(ExToBusMG{ii}));
 
    % Number of Existing Branches per DN
    nebr(ii) = max(mgd(:,1));
    
    lr = mgd(:,4)./Zb;
    lr(isnan(lr)) = [];
    % Existing Lines Resistance
    ExLinesR{ii} = lr;
    
    lx = mgd(:,5)./Zb;
    lx(isnan(lx)) = [];
    
    % Existing Lines Reactance
    ExLinesX{ii} = lx;
    
    fbmg = mgd(:,7);tbmg = mgd(:,8);
    fbmg(isnan(fbmg)) = [];tbmg(isnan(tbmg)) = [];
    CFromBusMG{ii} = fbmg;CToBusMG{ii} = tbmg;
    FromBusMG{ii} = [ExFromBusMG{ii};CFromBusMG{ii}];
    ToBusMG{ii} = [ExToBusMG{ii};CToBusMG{ii}];
    mgdd = mgd(:,6);
    CLines{ii} = mgdd(mgdd>0);
    Lines{ii} = [ExLines{ii};CLines{ii}];
    
    % Number of Candidate Nodes per DN
    ncn(ii) = max(max(max(CFromBusMG{ii}),max(CToBusMG{ii}))-nen(ii), 0);
    
    % Number of Candidate Branches per DN
    ncbr(ii) = max(mgd(:,6))-nebr(ii);
    
    % Total Number of Nodes per DN
    nn(ii) = nen(ii) + ncn(ii);
    
    % Total Number of Branches per DN
    nbr(ii) = nebr(ii) + ncbr(ii);
    
    lr = mgd(:,9)./Zb;
    lr(isnan(lr)) = [];
    % Candidate Lines Resistance
    CLinesR{ii} = lr;
    
    lx = mgd(:,10)./Zb;
    lx(isnan(lx)) = [];
    
    % Candidate Lines Reactance
    CLinesX{ii} = lx;    
    
    LinesR{ii} = [ExLinesR{ii};CLinesR{ii}];
    LinesX{ii} = [ExLinesX{ii};CLinesX{ii}];
    
    cll = mgd(:,11);
    cll(isnan(cll)) = [];
    % Candidate Lines Length
    CLinesLen{ii} = cll;
    
    clinc = mgd(:,12);
    clinc(isnan(clinc)) = [];
    
    % Candidate Lines Building Cost
    CLinesCost{ii} = clinc;
    
    lpu = mgd(:,13)./MVA;
    lpu(isnan(lpu)) = [];    
    
    % Lines Active Power Capacity Upper Limit
    ExLinesPu{ii} = lpu(1:nebr(ii));
    CLinesPu{ii} = lpu(nebr(ii)+1:end);
    
    % Substation Power Capacity Upper Limit
    SubPu{ii} = lpu(1);
    
    lpl = -lpu;
    lpl(isnan(lpl)) = [];
     
    % Lines Active Power Capacity Lower Lomit
    ExLinesPl{ii} = lpl(1:nebr(ii));
    CLinesPl{ii} = lpl(nebr(ii)+1:end);
    SubPl{ii} = lpl(1);
    
    lqu = lpu;
    lqu(isnan(lqu)) = [];
    
    % Lines Reactive Power Capacity Upper Limit
    ExLinesQu{ii} = lqu(1:nebr(ii));
    CLinesQu{ii} = lqu(nebr(ii)+1:end);
    
    lql = -lqu;
    lql(isnan(lql)) = [];
    
    % Lines Reactive Power Capacity Lower Limit
    ExLinesQl{ii} = lql(1:nebr(ii));
    CLinesQl{ii} = lql(nebr(ii)+1:end);
    
    vu = mgd(:,15);
    vu(isnan(vu)) = [];
    
    % Nodal Voltage Magnitude Upper Limit
    NVMu{ii} = vu;
    
    vl = mgd(:,16);
    vl(isnan(vl)) = [];
    
    % Nodal Voltage Magnitude Lower Limit
    NVMl{ii} = vl;  
    
    
    % Battery Eligible Nodes
    % Eligible nodes for Wind Capacity Installation
    batnodes = mgd(:,17);
    batnodes(isnan(batnodes)) = [];
    nbn(ii) = length(batnodes);
    
    BATL = zeros(nen(ii),nbn(ii));
    for bb=1:nbn(ii)
        BATL(batnodes(bb),bb) = 1;
    end
    BatLoc{ii} = BATL;
    
    Sinit = (mgd(:,18));
    Sinit(isnan(Sinit)) = [];
    
    % Battery Initial SOC per Distribution Network
    initialSOC{ii} = Sinit;   
  
    Sfinal = mgd(:,19);
    Sfinal(isnan(Sfinal)) = [];
    
    % Battery Final SOC per Distribution Network
    finalSOC{ii} = Sfinal;

    InfBus = mgd(:,20);
    InfBus(isnan(InfBus)) = [];
    
    % Number of Inflexible Loads per Distribution Network
    ninfl(ii) = length(InfBus);    
    
    % Inflexible Loads Locations per Distribution Network
    InfLoc{ii} = InfBus;    
    
    aInflex = tan(acos(mgd(:,21)));
    aInflex(isnan(aInflex)) = [];
    
    % Active to Reactive Power Conversion Parameter (tan(arcos(power factor))) of flexible Loads per Microgrid
    alphaInflex{ii} = aInflex;  
    
    % Local Load Consumption
    peakload = mgd(:,22)./MVA;
    peakload(isnan(peakload)) = [];
    InfLoadPerMG = zeros(ninfl(ii),T,TotalNumberOfScenarios);
    totalinflpermg = zeros(1,T,TotalNumberOfScenarios);
    for dd=1:ninfl(ii)
        for ss=1:TotalNumberOfScenarios
            InfLoadPerMG(dd,:,ss) = TotalSizingLoadScenarios(ss,:).*peakload(dd);
            totalinflpermg(:,:,ss) = sum(InfLoadPerMG(:,:,ss),1);
        end
    end
    InfLoad{ii} = InfLoadPerMG;
    
    TotalInfLoadPerMG{ii} = totalinflpermg;
    
    % Eligible nodes for Wind Capacity Installation
    wnodes = mgd(:,23);
    wnodes(isnan(wnodes)) = [];
    nw(ii) = length(wnodes);
    
    RESL = zeros(nn(ii),nw(ii));
    for rr=1:nw(ii)
        RESL(wnodes(rr),rr) = 1;
    end
    ResLoc{ii} = RESL;
    
    WindInt{ii} = TotalSizingWindScenarios{ii}';
    
    ARes = tan(acos(mgd(:,24)));
    ARes(isnan(ARes)) = [];
    
    % Active to Reactive Power Conversion Parameter (tan(arcos(power factor))) of Renewable Generators per Microgrid
    alphaRes{ii} = ARes;
    
    % Eligible nodes for Solar Capacity Installation
    pvnodes = mgd(:,25);
    pvnodes(isnan(pvnodes)) = [];
    npv(ii) = length(pvnodes);
    
    PVL = zeros(nn(ii),npv(ii));
    for rr=1:npv(ii)
        PVL(pvnodes(rr),rr) = 1;
    end    
    PVLoc{ii} = PVL; 
    
    PVInt{ii} = TotalSizingSolarScenarios{ii}';
    
    APV = tan(acos(mgd(:,26)));
    APV(isnan(APV)) = [];    
    
    % Active to Reactive Power Conversion Parameter (tan(arcos(power factor))) of Renewable Generators per Microgrid
    alphaRes{ii} = ARes;
    alphaPV{ii} = APV;    
    
end

InitialSOCPerMat = [];
ExLinesPuMat = [];
CLinesPuMat = [];
SubPuMat = [];
ExLinesPlMat = [];
CLinesPlMat = [];
SubPlMat = [];
ExLinesQuMat = [];
ExLinesQlMat = [];
CLinesQuMat = [];
CLinesQlMat = [];
NVMuMat = [];
NVMlMat = [];
WindIntMat = [];
PVIntMat = [];
TotalInfLoadMatPerMG = [];
CLinesCostMat = [];
CLinesLenMat = [];
LinesPuMat = [];
LinesPlMat = [];
LinesQuMat = [];
LinesQlMat = [];
% FromBusMGMat = [];
% ToBusMGMat = [];
FromNodesMat = [];
ToNodesMat = [];
LinesMat = [];
for ii=1:nmg
    TotalInfLoadMatPerMG = [TotalInfLoadMatPerMG;TotalInfLoadPerMG{ii}];
    WindIntMat(:,:,ii) = WindInt{ii};
    PVIntMat(:,:,ii) = PVInt{ii};
    InitialSOCPerMat = [InitialSOCPerMat;initialSOC{ii}*ones(nbn(ii),1)];
    ExLinesPuMat = [ExLinesPuMat;ExLinesPu{ii}];
    CLinesPuMat = [CLinesPuMat;CLinesPu{ii}];
    SubPuMat = [SubPuMat;SubPu{ii}];
    ExLinesPlMat = [ExLinesPlMat;ExLinesPl{ii}];
    CLinesPlMat = [CLinesPlMat;CLinesPl{ii}];
    SubPlMat = [SubPlMat;SubPl{ii}];
    ExLinesQuMat = [ExLinesQuMat;ExLinesQu{ii}];
    ExLinesQlMat = [ExLinesQlMat;ExLinesQl{ii}];
    CLinesQuMat = [CLinesQuMat;CLinesQu{ii}];
    CLinesQlMat = [CLinesQlMat;CLinesQl{ii}];    
    NVMuMat = [NVMuMat;NVMu{ii}];
    NVMlMat = [NVMlMat;NVMl{ii}];   
    CLinesCostMat = [CLinesCostMat;CLinesCost{ii}];
    LinesPuMat=[LinesPuMat;ExLinesPu{ii};CLinesPu{ii}];
    LinesPlMat=[LinesPlMat;ExLinesPl{ii};CLinesPl{ii}];
    LinesQuMat=[LinesQuMat;ExLinesQu{ii};CLinesQu{ii}];
    LinesQlMat=[LinesQlMat;ExLinesQl{ii};CLinesQl{ii}];   
    CLinesLenMat = [CLinesLenMat;CLinesLen{ii}];
%     FromBusMGMat = [FromBusMGMat;FromBusMG{ii}];
%     ToBusMGMat = [ToBusMGMat;ToBusMG{ii}];
    if ii==1
        FromNodesMat = [FromNodesMat;FromBusMG{ii}+1];
        ToNodesMat = [ToNodesMat;ToBusMG{ii}+1];
        LinesMat = [LinesMat;Lines{ii}];
    elseif ii>1
        FromNodesMat = [FromNodesMat;FromBusMG{ii}+sum(nn(ii-1))+ii];
        ToNodesMat = [ToNodesMat;ToBusMG{ii}+sum(nn(ii-1))+ii];
        LinesMat = [LinesMat;Lines{ii}+sum(nbr(ii-1))];
    end
end

WindTimeSeries = permute(WindIntMat,[1 3 2]);
SolarTimeSeries = permute(PVIntMat,[1 3 2]);
InfLoadMat = zeros(sum(nn),T,TotalNumberOfScenarios);
for ss=1:TotalNumberOfScenarios
    infld = InfLoad{1};
    InfLoadMat(InfLoc{1},:,ss) = infld(:,:,ss);
end
aInfMat = zeros(sum(nn),T);
aInfMat(InfLoc{1},:) = alphaInflex{1}*ones(1,T);

for ii=2:nmg
    for ss=1:TotalNumberOfScenarios
        infld = InfLoad{ii};
        InfLoadMat(InfLoc{ii}+nn(ii-1),:,ss) = infld(:,:,ss);
    end
    aInfMat(InfLoc{ii}+nn(ii-1),:) = alphaInflex{ii}*ones(1,T);
    
end

WindPF = alphaRes{1};
SolarPF = alphaPV{1};
KPMAT = [];
for mm=1:nmg
    KPMAT = blkdiag(KPMAT,ones(1,nen(mm)));
end

WMAT = [];
WWMAT = [];
PVMAT = [];
PPVMAT = [];
for ss = 1:TotalNumberOfScenarios
    for tt=1:T
        WW = [];
        WWWW = [];
        PV = [];
        PPV = [];
        for gg=1:nmg
            WW = blkdiag(WW,(WindTimeSeries(tt,gg,ss)*ones(1,nw(gg))));
            PV = blkdiag(PV,(SolarTimeSeries(tt,gg,ss)*ones(1,npv(gg))));
            WWWW = blkdiag(WWWW,(WindTimeSeries(tt,gg,ss)*ResLoc{gg}));
            PPV = blkdiag(PPV,(SolarTimeSeries(tt,gg,ss)*PVLoc{gg}));
        end
        WMAT = [WMAT;WW];
        WWMAT = [WWMAT;WWWW];
        PVMAT = [PVMAT;PV];
        PPVMAT = [PPVMAT;PPV];
    end
end

InitialSOCMat = zeros(sum(nen),1);

DischargeEff = Discharge_Efficiency*ones(sum(nbn),1);
ChargeEff = Charge_Efficiency*ones(sum(nbn),1);
FinalSOCPerMat = [];
for ii=1:nmg
    FinalSOCPerMat = [FinalSOCPerMat;finalSOC{ii}*ones(nbn(ii),1)];
end

B(:,RefBus) = [];
Yline = (LineFlows'.*LinesSusc)';
Yline(:,RefBus) = [];

cpup = zeros(nmg,T,TotalNumberOfScenarios);
cpdn = Price_Cap*ones(nmg,T,TotalNumberOfScenarios);

ConnectedLoadMat = zeros(sum(nn),T,TotalNumberOfScenarios);
IslandedLoadMat = zeros(sum(nn),T,TotalNumberOfScenarios);
for ii=1:nmg
    if ii==1
        ConnectedLoadMat(1:nen(ii),:,:) = InfLoadMat(1:nen(ii),:,:);
        IslandedLoadMat(nen(ii)+1:nen(ii)+ncn(ii),:,:) = InfLoadMat(nen(ii)+1:nen(ii)+ncn(ii),:,:);
    else
        ConnectedLoadMat(sum(nn(1:ii-1))+1:sum(nn(1:ii-1))+nen(ii),:,:) = InfLoadMat(sum(nn(1:ii-1))+1:sum(nn(1:ii-1))+nen(ii),:,:);
        IslandedLoadMat(sum(nn(1:ii-1))+nen(ii)+1:sum(nn(1:ii)),:,:) = InfLoadMat(sum(nn(1:ii-1))+nen(ii)+1:sum(nn(1:ii)),:,:);
    end
end
