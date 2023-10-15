% -----------------------------------------------------------------------------
% File: PV_BESS_OPTM.m
% Author: M. Rafey Tahir
% Date: 09th October 2023
% Reference: Kiran (USPCAS-E)
%  "Optimization of PV_BESS" 
% 
% -----------------------------------------------------------------------------
%%
clc
clear all
close all

P_load_E = randi([500 1000], 3, 24);
%%Cost parameters
E_charge=0.89;                      %electricty supply cahrgge
G_charge=0.78;                      % supply charge of
Cap_Cost_PV=  1500;                 %Capital Cost of PV
Rep_Cost_PV= 300;                   %replacement cost of PV
Maint_cost_PV= 1250;                %maintanance cost of PV
Cap_Cost_BESS=350;                  %Capital cost of BESS
Rep_Cost_BESS=200;                  %Replacement cost of BESS
Maint_cost_BESS=0;                  %maintainance cost of BESS
G_Cost=0.0468;                      %Gas price
%%Tarrif and interest
imp_T=48;                          %Import tariff
exp_T= 17;                         %Export Tariff
w= 0.08;                           %interest rate
e= 0.02;                           %escalation rate

%%Number of units
Kbs=3;                     % number of Battery energy storage units
Kpv=3;                      % Number of solar PV
%%time related parameters
Lu=  25;                 % life time of units
n=   20;                 %project life span
D=1;                     % Num of days in year( it is 1 as we are taking just one day demand)
Y=10;                    %replacement year
Q=21;                    %total num of cycle in one year of operation(random assume)
T=24;                    %T=24 as we are taking one day
Ru=10;                   %remaining life of unit   (Assume)
% Power/energy Related
% 1-Battery
n_charg=0.925;                            %charging efficiency
n_disch=0.725;                            %discharging efficeniey
Bess_Cap=17000;                           %Capacity of Battery in Watt-hour(Wh)
PV_Cap=10000;                             %Capacity of PV in Watt
Ebs=1000;                                 %Rated energy of Battery
%SOC indicateds State of charge of battery
SOC_min=0.2;
SOC_max=1;
MAX_BESS_charge=SOC_max*Bess_Cap;
MAX_BESS_discharge=SOC_min*Bess_Cap;

%% Variables
Ppv=optimvar('Power_PV',Kpv,T,'LowerBound',0,'UpperBound',PV_Cap) %W

Pbs=optimvar('Bess_power',Kbs,T,'LowerBound',0,'UpperBound',Bess_Cap)
P_charge=optimvar('BESS_charge',Kbs,T,'type','integer','LowerBound',0,'UpperBound',MAX_BESS_charge)  %Wh
P_dis=optimvar('BESS_discharge',Kbs,T,'type','integer','LowerBound',20,'UpperBound',MAX_BESS_discharge)
Pexp=optimvar('Power_exp',Kpv,T,'LowerBound',0,'UpperBound',5000) %W
Pimp=optimvar('Power_imp',Kpv,T,'LowerBound',0,'UpperBound',3000) %W

SOC = optimvar('Battery_SOC', Kbs, T, 'LowerBound', SOC_min * Bess_Cap, 'UpperBound', SOC_max * Bess_Cap);

%% Equations( formulas or equations are taken from paper)
u=(w-e)/1+e;      %utility interest calculation
CRF= ((1+u)^(n-1)) /u*(1+u)^n;  %Capital Recovery factor
%%
%annual electricity bill
ABe= E_charge*D+((imp_T*Pimp-exp_T*Pexp))
%net present cost of electricity
NPVe=CRF*ABe
show(NPVe)
%%
%%Annual Gas bill
% ABg= G_charge*D+((G_Cost*G_imp))        %annual gas bill
% NPVg=CRF*ABg
% show(NPVg)
%% Costs For PV
Cc_PV=Ppv*Cap_Cost_PV;       %net present discounted value of Capital cost per unit
Cm_PV=(Maint_cost_PV)/((1+w)*Lu);   %net present discounted value of Maintainance cost per unit
Cr_PV=Rep_Cost_PV*CRF;       %net present discounted value of Replacement cost per unit
%% Cost For BESS
Cc_BESS=Pbs*Cap_Cost_BESS;    %net present discounted value of Capital cost per unit
Cm_BESS=Maint_cost_BESS/(1+w)*Lu;  %net present discounted value of Maintannce cost per unit
Cr_BESS=Rep_Cost_BESS*CRF;    %net present discounted value of replacement cost per unit
%% total net present value of PV-BESS units(summation of all above)
NPVu=Kpv*(Cc_PV+Cr_PV+Cm_PV)+Kbs*(Cc_BESS+Cr_BESS+Cm_BESS);
%NPVt=NPVe+(NPVg)+NPVu          %total net present cost
NPVt=NPVe+NPVu;

% Call the function to solve the problem
result = myOptimizationFunction(NPVt, Kpv, T, P_load_E, Ppv, Pimp, Pexp, Pbs,Kbs,P_charge,P_dis,SOC,Bess_Cap,SOC_min,n_disch,n_charg);

% Access the solution and other information
solution = result.Solution;
objectiveValue = result.ObjectiveValue;
exitFlag = result.ExitFlag;

if exitFlag == 1
    disp(['Optimal Solution: ', num2str(solution)])
    disp(['Objective Value: ', num2str(objectiveValue)])
else
    disp('Optimization problem did not converge to a solution.')
end

% Display the results
% disp(['Optimal Solution: ', num2str(solution)])
% disp(['Objective Value: ', num2str(objectiveValue)])
% disp(['Exit Flag: ', num2str(exitFlag)])
Bess_power=solution.Bess_power;
Power_PV=solution.Power_PV;
Power_exp=solution.Power_exp;
Power_imp=solution.Power_imp;
plots
%%
function result = myOptimizationFunction(NPVt, Kpv, T, P_load_E, Ppv, Pimp, Pexp, Pbs,Kbs,P_charge,P_dis,SOC,Bess_Cap,SOC_min,n_disch,n_charg)
% Define optimization problem
prob = optimproblem;
x = optimvar('x', 1, 1, 'Type', 'integer', 'LowerBound', 0);
NPVt_total = sum(NPVt(:) * x');
prob.Objective = NPVt_total;
% Define constraints
econstraints1 = optimconstr(Kpv,T);

for i=1:T
    for j=1:Kpv
        if (i >= 1 && i <= 9) || (i >= 19 && i <= 24)
            econstraints1(j,i) = P_load_E(j,i) == Pimp(j,i) - Pexp(j,i) + Pbs(j,i);
        else
            econstraints1(j,i) = P_load_E(j,i) == Ppv(j,i) + Pimp(j,i) - Pexp(j,i) + Pbs(j,i);
        end
    end
end

prob.Constraints.econstraints1 = econstraints1;

% Battery charging constraint at mid-day
econstraints2 = optimconstr(Kbs, T);

for i = 1:T
    for j = 1:Kbs
        % Charging constraint
        P_charge>=0;
        econstraints2(j, i) = P_charge(j, i)==Ppv(j,i)-P_load_E(j,i);

    end
end

prob.Constraints.econstraints2 = econstraints2;

% Battery discharging constraint (for each unit and time period)
econstraints3 = optimconstr(Kbs, T);

for i = 1:T
    for j = 1:Kbs
        P_dis>=0;
        econstraints3(j, i) = P_dis(j, i) >= P_load_E(j, i) - Ppv(j, i);
    end
end

% Add the battery discharging constraints to the optimization problem
prob.Constraints.econstraints3 = econstraints3;

% Update SOC constraints to ensure feasibility
econstraints4 = optimconstr(Kbs, T);

for t = 1:T
    for j = 1:Kbs
        if t == 1
            % Initial SOC calculation considering discharge in the morning
            econstraints4(j, t) = SOC(j, t) == Bess_Cap * SOC_min + (Pbs(j, t) / n_disch) - P_dis(j, t) / n_disch;
        else
            % SOC calculation for the rest of the day
            econstraints4(j, t) = SOC(j, t) == SOC(j, t - 1) + (P_charge(j, t) / n_charg - P_dis(j, t) / n_disch);
        end
    end
end

% Add SOC constraints to the problem
prob.Constraints.econstraints4 = econstraints4;

% Import power constraint (for each unit and time period)
econstraints5 = optimconstr(Kpv, T);

for i = 1:T
    for j = 1:Kpv
        econstraints5(j, i) = Pimp(j, i) >= P_load_E(j, i) - Ppv(j, i) - P_dis(j, i);
    end
end

% Add the import power constraints to the optimization problem
prob.Constraints.econstraints5 = econstraints5;

% Export power constraint (for each unit and time period)
econstraints6 = optimconstr(Kpv, T);

for i = 1:T
    for j = 1:Kpv
        econstraints6(j, i) = Pexp(j, i) >= Ppv(j, i) - P_load_E(j, i) - P_charge(j, i);
    end
end

% Add the export power constraints to the optimization problem
prob.Constraints.econstraints6 = econstraints6;
% Solve the optimization problem
[sol, fval, exitflag, output] = solve(prob);

% Process the solution
result = struct('Solution', sol, 'ObjectiveValue', fval, 'ExitFlag', exitflag, 'Output', output);
end