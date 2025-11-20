%% Scenario Generation
ScenariosInputAndReduction

% Total Number Of Scenarios
TotalNumberOfScenarios = NumberOfScenarios*length(GenOffersScenarios);

Possibilities = zeros(1,TotalNumberOfScenarios);
pp = 0;
for ii=1:length(GenOffersScenarios)
    for jj=1:NumberOfScenarios
        pp = pp +1;
        Possibilities(pp) = 365*(GenOffersPossibilities(ii)*(RepDaysPoss(jj)./365));
    end
end

TotalSizingLoadScenarios = repmat(SizingLoadScenarios,length(GenOffersScenarios),1);
TotalSizingWindScenarios = cell(1,length(WindList2015));
TotalSizingSolarScenarios = cell(1,length(WindList2015));
for gg = 1:length(WindList2015)
    TotalSizingWindScenarios{gg} = repmat(SizingWindScenarios{gg},length(GenOffersScenarios),1);
    TotalSizingSolarScenarios{gg} = repmat(SizingSolarScenarios{gg},length(GenOffersScenarios),1);
end