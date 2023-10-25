function [ewi_value, ewifwd, statei] = f_ew(theta2m,eta_t,statei,ewifwd,ewi_value, ev0i,Data, Cst)

beta = Cst.beta;
nperiods_bf = Cst.nperiods_bf;
nperiods_af = Cst.nperiods_af;
ntrans = Cst.ntrans;
sigma = theta2m(1,5);
nphi = Cst.nphi;

%delta
delta_it= statei.BJ.delta_it; 


%Winning odds
pi_t = Data.wodd_BJ;

%income
inc_trans_mat = Data.inc_trans(:,:,1);
inc_groupi  = statei.inc_groupi;


%Iteration of value function
phi_trans_mat = zeros(nphi, nphi, nperiods_af);
trans_mat = zeros(4, nphi,nphi, nperiods_af);

%profile on
%Compute phi inclusive value

phi_it = zeros(nperiods_af,1);
ewifwd_new = zeros(nperiods_af,1);

phi_it_hi = zeros(nperiods_af,1);
phi_it_lo = zeros(nperiods_af,1);

for tt=1:nperiods_af

phi_t = eta_t + pi_t(tt).*squeeze(ev0i(4,1,nperiods_bf+1:end))';
phi_it(tt) = phi_t(tt,inc_groupi(nperiods_bf+tt));


for g=1:4
min_phi = floor(min(phi_t(:,g)));
max_phi = floor(max(phi_t(:,g)))+1; 
mid_phi = linspace(min_phi,max_phi,nphi)';
step = mid_phi(2)-mid_phi(1);
 
[phi_hi, phi_lo, pr_phi_hi] = assignStates(phi_t(tt,g),mid_phi,nphi,step,min_phi,max_phi,0);

[~,~,phifwd] = f_CalcPhiEvolution(phi_t(:,g),Data, Cst, mid_phi);
[phifwd_hi, phifwd_lo, pr_phifwd_hi] = assignStates(phifwd,mid_phi,nphi,step,min_phi,max_phi,0); 


for ip = 1:nphi

    simus_hi = squeeze(phifwd_hi(ip,:));
    uniq_hi = unique(simus_hi);
    for k= 1:length(uniq_hi)
    s = uniq_hi(k);
    phi_trans_mat(ip, s,tt) = sum(squeeze(pr_phifwd_hi(ip, simus_hi==s)));
    end
    
    simus_lo = squeeze(phifwd_lo(ip,:));
    uniq_lo = unique(simus_lo);
    for k= 1:length(uniq_lo)
    s = uniq_lo(k); 
    phi_trans_mat(ip, s, tt) = phi_trans_mat(ip, s, tt)+sum(squeeze(1-pr_phifwd_hi(ip, simus_hi==(s+1))));
    end
    
end

phi_trans_mat(:,:,tt) = phi_trans_mat(:,:,tt)/ntrans;

 
trans_mat(g,:,:,tt) = squeeze(phi_trans_mat(:,:,tt));
end

diffO = 20;
iterVO = 1;

while (diffO>Cst.tolV)
 ewi_value_gt = zeros(4,Cst.nphi);
 for g=1:4
 u1 = mid_phi+(1-pi_t(tt))*beta*reshape(ewi_value(g,:,tt), Cst.nphi,1);
 u0 = beta*reshape(ewi_value(g,:,tt),Cst.nphi,1);
 ewi_value_gt(g,:) = squeeze(trans_mat(g,:,:,tt))*(sigma*log(exp(u1/sigma)+exp(u0/sigma)));
 end
 %Integrate over incomes
 ewi_value(:,:,tt) = inc_trans_mat*ewi_value_gt;
 
 idx_hi = inc_groupi(tt)+(phi_hi-1)*4+(tt-1)*4*Cst.nphi;
 idx_lo = inc_groupi(tt)+(phi_lo-1)*4+(tt-1)*4*Cst.nphi;
 ewifwd_new(tt) = ewi_value(idx_hi).*pr_phi_hi + ewi_value(idx_lo).*(1-pr_phi_hi);
 diffO = max(abs(ewifwd_new(tt)-ewifwd(tt)));
 ewifwd(tt) = ewifwd_new(tt);
  
  if iterVO> Cst.iterMaxV
    break;
  end
  iterVO= iterVO+1;
end
end


statei.BJ.phi_it = phi_it;

end
