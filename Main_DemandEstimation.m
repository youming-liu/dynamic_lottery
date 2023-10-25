%================================================================
% Dynamic Model of Automobile License Allocaton
%================================================================
% Youming Liu
%================================================================

clear;
clc;
%set up cluster environment
%setenv('MATLAB_WORKER_ARGS','-p long --mem-per-cpu=20G --time=7-00:00:00');
%p= parpool('edith',64);

%add path
path(path,'Auxillary') 
path(path,'../Input')
path(path,'Output')
addpath(genpath('Auxillary/OpenOpt'));


%Specifications
FlagSpecs.PrePolicy=0;
FlagSpecs.isNoLast2=1;  %1, drop last two months before policy and first month after policy
FlagSpecs.isGroupSeg=1; %0, disaggreagted sales data; 1, aggregate sales to segment-type level

Cst.FlagSpecs = FlagSpecs;
%===============================================================   
    % Time: t1  = Jan, 2008  (Data Starts)
    %       t25 = Jan, 2011  (Lottery Begins in Beijing)
    %       t60 = Dec, 2012  (Data Ends)

%Indicate which markets to turn on (beijing, tianjin) (=B,T)
Cst.nmarkets = 2;

% Number of years and periods in data
Cst.nyears=5;
Cst.nyears_af = 2; %number of year after lottery
Cst.nyears_bf = Cst.nyears-Cst.nyears_af;
Cst.year0 =2007;
Cst.nmonths = 12;   
Cst.T    = 6;       % Validate Periods of Lottery Winners

Cst.nperiods   = Cst.nyears*Cst.nmonths;
Cst.nperiods_bf= Cst.nyears_bf*Cst.nmonths;
Cst.nperiods_af= Cst.nyears_af*Cst.nmonths;

%Tolerence level: for fixed point algorithm
Cst.tolV            = 1e-12;    % Value Function
Cst.tolT            = 1e-14;     % Overall tolerance

%Maximum Iterations
Cst.iterMaxV        = 1000;
Cst.iterMaxT        = 200;

%Heterogenity
Cst.nconsumers = 150; % Number of Household Types
Cst.nrand      = 5;  % Number of random coefficients
Cst.nincgroups = 4;  % Number of income groups (micro moments)



%===========================================================================%
%===========================Dynamic Parameters==============================%
%===========================================================================%

Cst.ntrans = 100; %Number of Simulations when Compute Transition of States
Cst.beta= 0.99;

%===================Vehicle Purchase Parameters=============================%
Cst.ndeltaf  = 21;          % Number of bins for flow utility of product you hold
Cst.ndeltafShr = 201;       % Number of bins for market shares
Cst.ndelta_it = 20;         % Number of bins for flow utility from purchasing new product:
Cst.min_delta_it  = [-4;-4; -4;-4];    % Min IVS
Cst.max_delta_it  = [1; 1; 1; 1];      % Max IVS

% State-space
%Cst.mid_delta_it  = linspace(Cst.min_delta_it,Cst.max_delta_it,Cst.ndelta_it)';



%===================Lottery Participation Parameters========================%
Cst.min_phi = -3;
Cst.max_phi = 1;
Cst.nphi    = 10;
Cst.mid_phi = linspace(Cst.min_phi,Cst.max_phi,Cst.nphi)';

%% Load and clean data; construct variables
Demand_ImportData 

%% GMM Estimation

%Initialize coefficients
%1st col: coefs for persistent shocks (coefs by 4 income groups)
%2nd col: coefs of price disutility (can vary across income groups) 
%3rd col: random coefs (price, cons, fuel, size, disp)
%4th col: owner coefs
%5th col: lottery coefs (sigma_omega)

%Matrix for parameter structure, 1 indicating the parameter should be
%estimated
Cst.theta2m = [ 1    0    1    1    1
                1    1    1    0    0
                1    1    0    0    0
                1    1    0    0    0
                0    0    1    0    0];                             
                                          
                                   
[Cst.thetai,Cst.thetaj] = find(Cst.theta2m);
theta2_ini = [0.1822    0.7671    1.1854    2.4826    0.00001    0.0364    0.1507    0.1937  -1.1760    0.6229    1    1.0923]';
theta_lb = [ 0, 0, 0, 0,  0, 0,  0,  -3, -3, -3, -5, 0];
theta_ub = [  5, 5,  5,  5,  15, 15,  15,  5,  5,  5, 5,  8];
options = optimoptions('fmincon','Display','iter', ...
        'SpecifyObjectiveGradient',false,'CheckGradients',false, ...
        'FiniteDifferenceType', 'forward', 'TolCon', 1e-4, ...
        'FiniteDifferenceStepSize', 1e-4,'StepTolerance',1e-4, ...
        'MaxFunEval', 1800, 'MaxIter', inf);
%     , 'HessianApproximation',...
%        'finite-difference','SubproblemAlgorithm','cg');
    %'Algorithm','active-set'
% options = optimoptions('fminunc','GradObj', 'on',  'DerivativeCheck', 'off');

fileTheta = 'thetapath.txt';
Cst.fileTheta = fileTheta;
fileID = fopen(fileTheta,'w');
fprintf(fileID, 'First stage \n');
fclose(fileID);

%First-stage
Cst.optw = 0;
Cst.iter = 0;

[theta2hat_stage1, fval_stage1] = fmincon('f_GMMobj', theta2_ini,[],[],[],[],...
     theta_lb,theta_ub,[],options, Data, Cst);

Cst.iter = 0;
[~, ~, Iter] = f_GMMobj(theta2hat_stage1, Data, Cst);
Data.OW_mom1 = Iter.OW_mom1;
Data.OW_mom2 = Iter.OW_mom2;

fileID = fopen(fileTheta,'a');
fprintf(fileID, 'Second stage \n');
fclose(fileID);

%Second-stage
Cst.optw = 1;

options = optimoptions('fmincon','Display','iter', ...
        'SpecifyObjectiveGradient',false,'CheckGradients',false, ...
        'FiniteDifferenceType', 'forward', 'TolCon', 1e-4, ...
        'FiniteDifferenceStepSize', 1e-4,'StepTolerance',1e-8, ...
        'MaxFunEval', 800, 'MaxIter', inf);

[theta2hat_stage2, fval_stage2] = fmincon('f_GMMobj', theta2hat_stage1,[],[],[],[],...
     theta_lb,theta_ub,[],options, Data, Cst);


[fval, paras, Iter] = f_GMMobj(theta2hat_stage2, Data, Cst);

%standard errors
[se, cov] = f_se(theta2hat_stage2, Data, Cst, Iter);

%put estimates to corresponding parameter groups
price_coefs  = -10+paras.lin_paras(end)+[0;theta2hat_stage2(5:7)];
se_price_coefs = paras.std_lin_paras(end)+[0; se(5:7)];
kappas = theta2hat_stage2(1:4);
replace = 1/(1+exp(theta2hat_stage2(end-1)));
se_replace = exp(theta2hat_stage2(end-1))/(1+exp(theta2hat_stage2(end-1)))^2*se(end-1);
sigmas = theta2hat_stage2(8:10);
se_sigmas = se(8:10);
logit_err = theta2hat_stage2(end);
se_err = se(end);
linparas = paras.lin_paras(1:end-1);
se_linparas = paras.std_lin_paras(1:end-1);
save Output/EstResult.mat


