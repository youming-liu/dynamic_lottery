function [ ev0i_value, ev0, ev0fwd] = f_ev0(theta2m,statei, Data, Cst)
beta = Cst.beta;
nperiods    = Cst.nperiods;
T = Cst.T;
nmonths= Cst.nmonths;
month  = Data.month;
ndelta_it   = Cst.ndelta_it;
Time = [1:nperiods]';

%Kappa
kappa  = theta2m(1:4,1);

%income
inc_trans_mat = Data.inc_trans(:,:,1);
inc_groupi  = statei.inc_groupi;

delta_hi = statei.BJ.delta_hi;
delta_lo = statei.BJ.delta_lo;
pr_delta_hi = statei.BJ.pr_delta_hi;
delta_it_hi = delta_hi(Time + (inc_groupi-1)*nperiods);
delta_it_lo = delta_lo(Time + (inc_groupi-1)*nperiods);
pr_delta_it_hi = pr_delta_hi(Time + (inc_groupi-1)*nperiods); 

deltafwd_hi = statei.BJ.deltafwd_hi;
deltafwd_lo = statei.BJ.deltafwd_lo;
pr_deltafwd_hi = statei.BJ.pr_deltafwd_hi;


ev0i_value = zeros(T,4,ndelta_it,nmonths);
v0i_value = zeros(T,4,ndelta_it,nmonths);

%Value function iterations
for ss=T:-1:1
  %utility of not purchasing
      u0 = squeeze(beta*ev0i_value(ss,:,:,:));
      %utility of purchasing
      mid_delta_it = statei.BJ.mid_delta;
      u1 = kappa+mid_delta_it';
      v0i_value(ss,:,:,:) = log(exp(u0)+exp(u1));
   
   if ss==1
   continue; 
   end

   for tt=1:nmonths

  
  %integrate over deltas
   ev0i_value_gt = zeros(4,ndelta_it); %next period value functions integrated over deltas
     for gg=1:4
      
      if tt<nmonths
      v0i_value_gt = squeeze(v0i_value(ss,gg,:,tt+1));
      else
      v0i_value_gt = squeeze(v0i_value(ss,gg,:,1));    
      end
      deltafwd_hi_t = deltafwd_hi(:,:,tt,gg);
      deltafwd_lo_t = deltafwd_lo(:,:,tt,gg);
      pr_deltafwd_hi_t = pr_deltafwd_hi(:,:,tt,gg);
      ev0i_value_gt(gg,:) = mean(v0i_value_gt(deltafwd_hi_t).*pr_deltafwd_hi_t+v0i_value_gt(deltafwd_lo_t).*(1-pr_deltafwd_hi_t),2);
     end
  %integrate over income transitions
    ev0i_value(ss-1,:,:,tt) = inc_trans_mat*ev0i_value_gt;
   end

end

ev0 = zeros(T,4,nperiods);
ev0fwd = zeros(T,nperiods);
for s=1:T
   for g=1:4
   delta_g_hi= delta_hi(Time + (g-1)*nperiods);
   delta_g_lo= delta_lo(Time + (g-1)*nperiods);
   ev0i_value_hi = v0i_value(s+(g-1)*T+(delta_g_hi-1)*T*4+(month-1)*ndelta_it*T*4);
   ev0i_value_lo = v0i_value(s+(g-1)*T+(delta_g_lo-1)*T*4+(month-1)*ndelta_it*T*4);
   ev0(s,g,:) = pr_delta_it_hi.*ev0i_value_hi+(1-pr_delta_it_hi).*ev0i_value_lo;
   end
   
   ev0ifwd_hi = v0i_value(s+(inc_groupi-1)*T+(delta_it_hi-1)*T*4+(month-1)*ndelta_it*T*4);
   ev0ifwd_lo = v0i_value(s+(inc_groupi-1)*T+(delta_it_lo-1)*T*4+(month-1)*ndelta_it*T*4);
   ev0fwd(s,:) = pr_delta_it_hi.*ev0ifwd_hi+(1-pr_delta_it_hi).*ev0ifwd_lo;
end

end
