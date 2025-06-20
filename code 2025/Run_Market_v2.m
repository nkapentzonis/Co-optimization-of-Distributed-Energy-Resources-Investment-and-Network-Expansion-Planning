final_TSO_market = optimproblem('ObjectiveSense','min');
g = optimvar('g',ng,T,'LowerBound',GMIN,'UpperBound',GMAX); % (TN) Generation variable - generation bounds (b.3)
d = optimvar('d',nd,T); % (TN) demand variable - demand bounds (b.6)
pup = optimvar('pup',nmg,T,'LowerBound',0); % DN power sold to market
pdn = optimvar('pdn',nmg,T,'LowerBound',0); % DN power bought from market
theta = optimvar('theta',n-1,T); % TN voltage angle variable
o = x_opt.pup;
b = x_opt.pdn;
% pup = optimvar('pup',nmg,T,'LowerBound',0); % DN power sold to market
% pdn = optimvar('pdn',nmg,T,'LowerBound',0);
for ww=1:TotalNumberOfScenarios
    d.LowerBound = DMIN(:,:,ww);
    d.UpperBound = DMAX(:,:,ww);
    TSOPowerBalance = optimconstr(n,T);
    for tt=1:T
        TSOPowerBalance(:,tt) = -GenLoc*g(:,tt)+DemLoc*d(:,tt)-MGLoc*(pup(:,tt)-pdn(:,tt))+B*theta(:,tt)==0;
    end
    final_TSO_market.Constraints.system_power_balance = TSOPowerBalance;
    final_TSO_market.Constraints.upper_bound_on_ESP = pup(:,:)== o(:,:,ww);
    final_TSO_market.Constraints.lower_bound_on_ESP = pdn(:,:)== b(:,:,ww);
    
    % Line Flows Bounds (b.6)
    LineFlowsUpperBound = optimconstr(nl,T);
    LineFlowsLowerBound = optimconstr(nl,T);
    for tt=1:T
        LineFlowsUpperBound(:,tt) = Yline*theta(:,tt)<=FMAX;
        LineFlowsLowerBound(:,tt) = -Yline*theta(:,tt)<=FMAX;
    end
    
    final_TSO_market.Constraints.upper_bound_on_line_flows = LineFlowsUpperBound;
    final_TSO_market.Constraints.lower_bound_on_line_flows = LineFlowsLowerBound;
    
    % Ramp Up/Down Constraints (not in manuscript so far!!!!!)
    final_TSO_market.Constraints.ramp_up_t1 = g(:,1)<=repmat(G0+RU,1,1);
    final_TSO_market.Constraints.ramp_up_t2 = g(:,2:T)-g(:,1:T-1)<=repmat(RU,1,T-1);
    
    final_TSO_market.Constraints.ramp_dn_t1 = -g(:,1)<=repmat(-G0+RD,1,1);
    final_TSO_market.Constraints.ramp_dn_t2 = g(:,1:T-1)-g(:,2:T)<=repmat(RD,1,T-1);
    
    final_TSO_market.Objective = sum(sum(GenBids(:,:,ww).*g,1),2)-sum(sum(DemBids(:,:,ww).*d,1),2)+sum(sum(cpup(:,:,ww).*pup))-sum(sum(cpdn(:,:,ww).*pdn));    
    
    [x_TSO{ww},BF(ww),exitflag2{ww},output2{ww},duals2{ww}]=solve(final_TSO_market);
    Prices2{ww} = duals2{ww}.Constraints.system_power_balance;    
end