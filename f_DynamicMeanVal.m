function[Iter,Cst]=f_DynamicMeanVal(theta2m, Data, Cst, Iter)

%Step 1: initialize delta, phi and value functions
%Step 2: value function iterations
%Step 3: Update delta, phi
%Step 4: Compute market shares and update deltaf
%Iterating between step 2 to 4 until convergence


%Load mean utilities and stack up the utilities from both stages
deltaf_jt = Iter.deltaf_jt(:,1:Cst.nmarkets);
eta_t  = Iter.eta_t;
Iter.x = [deltaf_jt(:);eta_t];

%Start iteration
try 
[Iter]=fp_squarem(@(Iter)fp(theta2m, Data, Cst, Iter), Iter,  'noisy',2,'con_tol',Cst.tolT,'max_iter',Cst.iterMaxT,'algorithm',1);
Iter.deltaf_jt = reshape(Iter.x(1:Cst.nprods*Cst.nmarkets), Cst.nprods, Cst.nmarkets);
Iter.eta_t     = Iter.x(Cst.nprods*Cst.nmarkets+1:end);
catch ME
 warning('Found a NAN value.');
 Iter.flag_NAN = true;   
end

end


%% Define the fixed point function 
function [Iter] = fp(theta2m, Data, Cst, Iter)

%Step 1. Initialize mean utilities and value functions
x = Iter.x;

%Full value functions
ev_value = Iter.ev_value;
ev0_value = Iter.ev0_value;
ew_value = Iter.ew_value;

%Value functions evaluated at each period
ewfwd    = Iter.ewfwd;
evfwd  = zeros(2,4,Cst.nperiods,Cst.nmarkets, Cst.nconsumers);
ev0  = zeros(Cst.T, 4, Cst.nperiods, Cst.nconsumers);
ev0fwd = zeros(Cst.T, Cst.nperiods, Cst.nconsumers);

%price utility for all income groups
rho     = Data.rho; 
%price utility evaluated by income draws
rho_ijt = Data.rho_ijt; 

%mean utilities for vehicles
deltaf_jt = reshape(x(1:Cst.nprods*Cst.nmarkets), Cst.nprods, Cst.nmarkets);
%utilities with preference heterogeneity
deltaf_ijt = repmat(deltaf_jt,1,1,Cst.nconsumers) + Data.mu_ijt;

%mean utilities for lottery participation
eta_t     = x(Cst.nprods*Cst.nmarkets+1:end);

%initial delta and phi
delta_it = zeros(Cst.nperiods, Cst.nmarkets, Cst.nconsumers);
phi_it   = zeros(Cst.nperiods_af, Cst.nconsumers);

%matrix for shares/particiaption rates/ownership
share_ijt = zeros(Cst.nprods,Cst.nmarkets,Cst.nconsumers);
first_ijt = zeros(Cst.nprods,Cst.nmarkets,Cst.nconsumers);
parate_it  = zeros(Cst.nperiods_af,Cst.nconsumers);
owner_it = zeros(Cst.nperiods,Cst.nmarkets,Cst.nconsumers);

parfor ii=1:Cst.nconsumers

%Step 2: value function iterations

%Update ev_value and delta_it
 [ev_value(:,:,:,:,:,ii), evi, evifwd, statei] = f_ev(theta2m, deltaf_ijt(:,:,ii), rho_ijt(:,:,ii), rho(:,:,:,ii), ev_value(:,:,:,:,:,ii), Data, Cst);

if Cst.nmarkets==1 %if there is only one market, by default it is Beijing
  delta_it(:,:,ii)  = statei.BJ.delta_it;
else
  delta_it(:,:,ii) = [statei.BJ.delta_it, statei.TJ.delta_it];
end
evfwd(:,:,:,:,ii) = evifwd;

 statei.inc_groupi = squeeze(Data.inc_group_t(:,1,ii));

%Update ev0_value
 [ev0_value(:,:,:,:,ii), ev0i, ev0ifwd] = f_ev0(theta2m, statei, Data, Cst);
 ev0(:,:,:,ii) =ev0i;
 ev0fwd(:,:,ii) = ev0ifwd;
 
%Update ew and phi
 ewifwd = ewfwd(:,ii);
 ewi_value = ew_value(:,:,:,ii);
 [ewi_value, ewifwd, statei] = f_ew(theta2m,eta_t,statei,ewifwd, ewi_value, ev0i,Data, Cst);
 ew_value(:,:,:,ii) = ewi_value;
 ewfwd(:,ii) = ewifwd;
 phi_it(:,ii) = statei.BJ.phi_it;
 
 %Step 3: Calculate market share, and update deltaf
 statei.owni = Data.owner_draws(ii,:);
 statei.inci = squeeze(Data.inc_draws_t(:,:,ii));
 statei.incgi = squeeze(Data.inc_group_t(:,:,ii));
 [share_ijt(:,:,ii), parate_it(:,ii), first_ijt(:,:,ii),owner_it(:,:,ii)] = f_DynamicShareHat(theta2m, deltaf_ijt(:,:,ii), rho_ijt(:,:,ii), statei, evifwd, ev0ifwd, ewifwd, Data, Cst);

end

sharehat      =  mean(share_ijt,3);
paratehat     =  mean(parate_it,2);

deltaf_jt_new    = deltaf_jt +(log(Data.share(:,1:Cst.nmarkets))-log(sharehat));
eta_t_new       = eta_t  +(log(Data.parate_bj)-log(paratehat));

%update mean utilities
deltaf_jt = deltaf_jt_new;
eta_t  = eta_t_new;

%stack up mean utilities
x = [deltaf_jt(:);eta_t];
Iter.x  = x;
Iter.parate_it = parate_it;
Iter.share_ijt  = share_ijt;
Iter.first_ijt  = first_ijt;
Iter.owner_it   = owner_it;
Iter.deltaf_jt = deltaf_jt;
Iter.eta_t     = eta_t;
Iter.delta_it  = delta_it;
Iter.phi_it    = phi_it;
Iter.ew_value  = ew_value;
Iter.ev_value  = ev_value;
Iter.ev0_value = ev0_value;
Iter.ewfwd     = ewfwd;
Iter.ev0  = ev0;
Iter.ev0fwd = ev0fwd;
Iter.evfwd  = evfwd;
end
