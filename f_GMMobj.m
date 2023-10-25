function [fval, paras, Iter] = f_GMMobj(theta2, Data, Cst)

%Print the current parameters
printm(theta2');

%set up timer
stTime   = tic;

%setup parameter structure
theta2m = Cst.theta2m;
theta2m(Cst.thetai+ size(theta2m,1)*(Cst.thetaj-1))= theta2;


%Initialize mean utilities and value functions
if Cst.iter==0 
%If not to use the mean utility and value functions from previous iteration
    shares       = Data.share(:,1:Cst.nmarkets);
    Tshares      = 1-Data.TimeProduct*shares(:,1:Cst.nmarkets);
    Iter.deltaf_jt = (log(shares)-log(Tshares(Data.Time,:)));
    Iter.eta_t     = log(Data.parate_bj)-log(1-Data.parate_bj);

%Initialize value functions
    %Vehicle Purchases w/o lottery  
    Iter.ev_value = zeros(2,4,Cst.ndelta_it,Cst.nmonths,Cst.nmarkets, Cst.nconsumers);%First state: owner and non-owner; second state income groups  
    %Vehicle Purchases w/ lottery
    Iter.ev0_value = zeros(Cst.T, 4, Cst.ndelta_it, Cst.nmonths, Cst.nconsumers);
    
    Iter.ew_value = zeros(4, Cst.nphi, Cst.nperiods_af,Cst.nconsumers);%First state; income groups
    Iter.ewfwd    = zeros(Cst.nperiods_af, Cst.nconsumers);
else
%Load from previous iteration
load 'Interim/Iter_save.mat';
end

%Compute random Utilities
%Individual preference
Data.mu_ijt = f_mu(theta2m,Data,Cst);

%Price disutility
[Data.rho, Data.rho_ijt]= f_rho(theta2m,Data,Cst);    


%Initialize iteration checkers
Iter.flag_NAN= false;
Iter.flag_maxBLP = false;

%BLP contraction mapping
[Iter,Cst]=f_DynamicMeanVal(theta2m, Data, Cst, Iter);


if ~Iter.flag_NAN && ~Iter.flag_maxBLP
%If the contraction mapping does not return to NaN values or not converge
    %Compute moments
    [moments, paras, Cst]= f_moments(theta2m, Data, Cst, Iter);
    %Cst.Iter = 1;
    if Cst.optw==0 %if optimal weight matrix is not used
    Iter.OW_mom1 = moments.OW_mom1;
    Iter.OW_mom2 = moments.OW_mom2;
    end


%% GMM Value
    gmm1 = (Cst.nprods*Cst.nmarkets)*moments.mom1*(moments.W_mom1\(moments.mom1.'));
    gmm2 = Cst.nperiods_af*moments.mom2*(moments.W_mom2\(moments.mom2.'));
    gmm3 = moments.mom3'*(moments.W_mom3\(moments.mom3));
    fval = gmm1+gmm2+gmm3;
    
    if ~isnan(fval) 
    if ~Iter.flag_maxBLP
    save 'Interim/Iter_save.mat' 'Iter' 'moments';
    end

    theta1 = paras.lin_paras;
    nK1 = length(theta1);
    nK2 = length(theta2);
    fileID = fopen(Cst.fileTheta,'a');
    format{1} = ['%6.2f \r\n'];
    for i=2:nK1
    format{1} = strcat('%6.2f ', format{1});
    end
    format{2} = ['%6.3f \r\n'];
    for i=1:nK2
    format{2} = strcat('%6.3f  ', format{2});
    end
    format{3} = ['%6.2f  %6.2f  %6.2f  %6.2f\r\n'];
    fprintf(fileID, 'linear paras:');
    fprintf(fileID, format{1}, theta1);
    fprintf(fileID, 'nonlinear paras:');
    fprintf(fileID, format{2}, theta2);
    fprintf(fileID, 'gmms:');
    fprintf(fileID, format{3}, [gmm1*1000, gmm2*1000, gmm3, fval]);
    fclose(fileID);
    end
else
% In this case, the program diverged so we output GMM Value = infinity
    fval = 1e6;
end

 elapTime = toc(stTime);
 fprintf('Time elapsed: %8.2f obj = %8.2f\n',elapTime,fval);
end
