function [evi_value, ev, evfwd, statei] = f_ev(theta2m,deltaf_ijt, rho_ijt, rho, evi_value, Data, Cst)

beta = Cst.beta;
nperiods = Cst.nperiods;
nmonths= Cst.nmonths;
nmarkets = Cst.nmarkets;
nincgroups = Cst.nincgroups;
ntrans   = Cst.ntrans;
dmonthly = Data.d_monthly;
ndelta_it   = Cst.ndelta_it;

%Kappa
kappa  = theta2m(1:4,1);

evfwd = zeros(2,4,nperiods,nmarkets);
ev    = zeros(2,4,nperiods,nmarkets);



for mm=1:nmarkets

%State Space
%delta
delta = log(Data.TimeProduct*exp(deltaf_ijt(:,mm)+squeeze(rho(:,mm,:))));
delta_it = log(Data.TimeProduct*exp(deltaf_ijt(:,mm)+squeeze(rho_ijt(:,mm))));
deltafwd = zeros(ndelta_it, ntrans, nmonths, nincgroups);
delta_hi = zeros(nperiods,nincgroups); delta_lo = zeros(nperiods,nincgroups); pr_delta_hi = zeros(nperiods, nincgroups);
deltafwd_hi = zeros(ndelta_it, ntrans, nmonths,nincgroups); deltafwd_lo = zeros(ndelta_it, ntrans, nmonths,nincgroups); pr_deltafwd_hi = zeros(ndelta_it, ntrans,nmonths,nincgroups);
delta_trans_mat = zeros(ndelta_it*nmonths, ndelta_it*nmonths, nincgroups);
mid_delta = zeros(ndelta_it, 4);



for g=1:nincgroups


%Delta grid
min_delta_it = Cst.min_delta_it(g);%min(delta(:,g)-1);
max_delta_it =Cst.max_delta_it(g);%max(delta(:,g)+1);%
mid_delta_it =linspace(min_delta_it, max_delta_it, ndelta_it)';
mid_delta(:,g) = mid_delta_it;
step = mid_delta_it(2)-mid_delta_it(1);

[~,~,deltafwd(:,:,:,g)] = f_CalcDeltaEvolution(delta(:,g), Data, Cst, mid_delta_it);

[delta_hi(:,g), delta_lo(:,g), pr_delta_hi(:,g)] = assignStates(delta(:,g),mid_delta_it,ndelta_it,step,min_delta_it,max_delta_it,0);  
[deltafwd_hi(:,:,:,g), deltafwd_lo(:,:,:,g), pr_deltafwd_hi(:,:,:,g)] = assignStates(deltafwd(:,:,:,g),mid_delta_it,ndelta_it,step,min_delta_it,max_delta_it,0);  



for id=1:ndelta_it
for tt=1:nmonths
ss= tt<nmonths;
sims_hi = squeeze(deltafwd_hi(id,:,tt,g));
sims_lo = squeeze(deltafwd_lo(id,:,tt,g));

uniq_hi = unique(sims_hi);
uniq_lo = unique(sims_lo);

for idx = uniq_hi
delta_trans_mat((tt-1)*ndelta_it+id,ss*tt*ndelta_it+idx,g) = sum(pr_deltafwd_hi(id, sims_hi==idx, tt,g));
end
for idx = uniq_lo
delta_trans_mat((tt-1)*ndelta_it+id,ss*tt*ndelta_it+idx,g) = delta_trans_mat((tt-1)*ndelta_it+id,ss*tt*ndelta_it+idx,g) +sum(1-pr_deltafwd_hi(id, sims_hi==(idx+1), tt,g));
end


end
end
end
delta_trans_mat = delta_trans_mat/100;

%income
inc_trans_mat = Data.inc_trans(:,:,mm);


% Value function iterations
evi_value0 = squeeze(evi_value(1,:,:,:,mm));
evi_value1 = squeeze(evi_value(2,:,:,:,mm));


diff = 1e20;
iterV = 1;
while(diff>Cst.tolV)

   %nonowner, not purchase
    v_value00 = beta*evi_value0;
   %nonowner, purchase
    v_value01 = kappa+repmat(mid_delta_it',4,1);
   %owner, not purchase
    %v_value10 = kappa+beta*evi_value1;
   v_value10 = beta*evi_value1; %(1-beta)*own+
   %owner, purchase
    %v_value11 = kappa+own+repmat(mid_delta_it',4,1);
   v_value11 = repmat(mid_delta_it',4,1);
   %Take expecation over epsilons 
    v_value0 = log(exp(v_value00)+exp(v_value01));
    v_value1 = log(exp(v_value10)+exp(v_value11));
    
   %Take expecation over states
   evi_value0g = zeros(4, ndelta_it*nmonths);
   evi_value1g = zeros(4, ndelta_it*nmonths);
   v_value0_rs = reshape(v_value0(:,:,:), nincgroups, ndelta_it*nmonths);
    v_value1_rs = reshape(v_value1(:,:,:), nincgroups, ndelta_it*nmonths);
   for g=1:nincgroups
    evi_value0g(g,:) = delta_trans_mat(:,:,g)*v_value0_rs(g,:)';
    evi_value1g(g,:) = delta_trans_mat(:,:,g)*v_value1_rs(g,:)';
   end

    evi_value0_new = reshape(inc_trans_mat*evi_value0g,4,ndelta_it,nmonths);
    evi_value1_new = reshape(inc_trans_mat*evi_value1g,4,ndelta_it,nmonths);

    diff1 = mean((evi_value0_new(:)-evi_value0(:)).^2);
    diff2 = mean((evi_value1_new(:)-evi_value1(:)).^2);
    diff = max(diff1,diff2);
    
    evi_value0 = evi_value0_new;
    evi_value1 = evi_value1_new;
    
    if iterV>Cst.iterMaxV
    fprintf("Value Function does not coverge, difference = %3.12f \n", diff);
    break;
    end
    

    iterV= iterV+1;
end
evi_value(1,:,:,:,mm) = evi_value0;
evi_value(2,:,:,:,mm) = evi_value1;

for i=1:2
for g=1:4
evi_value_hi = sum(squeeze(evi_value(i,g,delta_hi(:,g),:,mm)).*dmonthly,2);
evi_value_lo = sum(squeeze(evi_value(i,g,delta_lo(:,g),:,mm)).*dmonthly,2);
evfwd(i,g,:,mm) = pr_delta_hi(:,g).*evi_value_hi+(1-pr_delta_hi(:,g)).*evi_value_lo;
u0 = (i-1)*kappa(g)+beta*squeeze(evfwd(i,g,:,mm)); %(i-1)*(1-beta)*own+
u1 = delta(:,g) + kappa(g) ;
ev(i,g,:,mm) = log(exp(u0)+exp(u1));
end
end

%save state space variables
if mm==1
statei.BJ.delta_it = delta_it; statei.BJ.deltafwd = deltafwd;
statei.BJ.delta_hi = delta_hi; statei.BJ.delta_lo = delta_lo; statei.BJ.pr_delta_hi = pr_delta_hi;
statei.BJ.deltafwd_hi = deltafwd_hi; statei.BJ.deltafwd_lo = deltafwd_lo; statei.BJ.pr_deltafwd_hi = pr_deltafwd_hi;
statei.BJ.delta_trans_mat = delta_trans_mat;
statei.BJ.mid_delta = mid_delta;
else
statei.TJ.delta_it = delta_it; statei.TJ.deltafwd = deltafwd;

statei.TJ.delta_hi = delta_hi; statei.TJ.delta_lo = delta_lo; statei.TJ.pr_delta_hi = pr_delta_hi;
statei.TJ.deltafwd_hi = deltafwd_hi; statei.TJ.deltafwd_lo = deltafwd_lo; statei.TJ.pr_deltafwd_hi = pr_deltafwd_hi;
statei.TJ.delta_trans_mat = delta_trans_mat;
statei.TJ.mid_delta = mid_delta;
end
end

end
