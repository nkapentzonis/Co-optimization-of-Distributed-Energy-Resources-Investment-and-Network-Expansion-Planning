%% Solve Model
options = optimoptions('intlinprog','Display','iter');

% showproblem(sizing_op); 
[x_opt,neg_profits,exitflag,output] = solve(sizing_op,'solver','intlinprog');
