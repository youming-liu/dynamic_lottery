function [se, cov] = f_se(theta2, Data, Cst, Iter)

%% Initialization
nprods = Cst.nprods;
nperiods = Cst.nperiods;
nperiods_af = Cst.nperiods_af;

% Size of theta2
ntheta = length(theta2);
%recover theta2 matrix
theta2m = Cst.theta2m;
theta2m(Cst.thetai+ size(theta2m,1)*(Cst.thetaj-1))= theta2;

% Set derivative step
tolx = 10^-6;

%Use Coarse Tolerance Compute the Std
Cst.tolDelta        = 1e-10;    % BLP fixed point
Cst.tolV            = 1e-12;    % Value Function
Cst.tolT            = 1e-10;     % Overall tolerance
Cst.iterMaxT        = 500;
Cst.flag_NAN        = 0;
Cst.flag_maxBLP     = 0;
%% Get base moments
%Random Utilities
    % Get Individual preference
    Data.mu_ijt = f_mu(theta2m,Data,Cst);

    % Get price disutility
    [Data.rho, Data.rho_ijt]= f_rho(theta2m,Data,Cst);

%BLP contraction mapping
   [Iter,Cst]=f_DynamicMeanVal(theta2m, Data, Cst, Iter);

%Compute moments
[m, ~, ~]= f_moments(theta2m, Data, Cst, Iter);   

%% Get adjusted moments  

% Create gradient matrix
G1 = zeros(size(m.mom1,2),ntheta);
G2 = zeros(size(m.mom2,2),ntheta);
G3 = zeros(size(m.mom3,1),ntheta);
Cst.tolT            = 1e-11;     % Overall tolerance
for kk=1:ntheta
tolx = 10^-6;
%Upper
theta2_up = theta2;
theta2_up(kk) = theta2(kk).*(1+tolx);

%recover theta2 matrix
theta2m_up = theta2m;
theta2m_up(Cst.thetai+ size(theta2m,1)*(Cst.thetaj-1))= theta2_up;
    
%Random Utilities
    % Get Individual preference
    Data.mu_ijt = f_mu(theta2m_up,Data,Cst);

    % Get price disutility
    [Data.rho, Data.rho_ijt]= f_rho(theta2m_up,Data,Cst); 

%BLP contraction mapping
   [Iter_up,Cst]=f_DynamicMeanVal(theta2m_up, Data, Cst, Iter);

%Compute moments
[m_up, ~, ~]= f_moments(theta2m_up, Data, Cst, Iter_up);   

%Lower
theta2_lo = theta2;
theta2_lo(kk) = theta2(kk).*(1-tolx);

%recover theta2 matrix
theta2m_lo = theta2m;
theta2m_lo(Cst.thetai+ size(theta2m,1)*(Cst.thetaj-1))= theta2_lo;

%Random Utilities
    % Get Individual preference
    Data.mu_ijt = f_mu(theta2m_lo,Data,Cst);

    % Get price disutility
    [Data.rho, Data.rho_ijt]= f_rho(theta2m_lo,Data,Cst); 
    %Data.rho_ijt= f_rho(theta2m_lo,Data,Cst);  

%BLP contraction mapping
   [Iter_lo,Cst]=f_DynamicMeanVal(theta2m_lo, Data, Cst, Iter);

%Compute moments
[m_lo, ~, ~]= f_moments(theta2m_up, Data, Cst, Iter_lo); 

%Calculate numeric gradient
G1(:,kk) = 0.5*(m_up.mom1-m_lo.mom1)'./(theta2(kk)*tolx);
G2(:,kk) = 0.5*(m_up.mom2-m_lo.mom2)'./(theta2(kk)*tolx);
G3(:,kk) = 0.5*(m_up.mom3-m_lo.mom3)'./(theta2(kk)*tolx);
fprintf('Theta %i complete\n',kk);
end

% Step 3: Create Var_Cov Matrix
% cov_Tehta = (G'AG)^-1G'A*Omega*A*G*(G'AG)^-1
% G is the graident matrix
% A is weighting matrix (inv(Cov))
% Omega is variance matrix of moments divided by number of observations
save Output/SE ;

% Variance of theta
var1 = (G1'/(nprods*2)) * inv(m.OW_mom1) * (G1/(nprods*2));
var2 = (G2'/nperiods_af * inv(m.OW_mom2) * G2/nperiods_af);
var3 = (G3' * inv(m.W_mom3) * (G3))/max(Data.obs3);

cov = inv(var1 + var2+ var3)/(nprods*2);
se = sqrt(diag(cov));

end
