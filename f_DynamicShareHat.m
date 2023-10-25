function [s_ijt, l_it, f_ijt, o_it] = f_DynamicShareHat(theta2m, deltaf_ijt,rho_ijt, statei, evifwd,  ev0ifwd, ewifwd, Data, Cst)
beta = Cst.beta;
nmarkets = Cst.nmarkets;
nperiods = Cst.nperiods;
nprods   = Cst.nprods;
nperiods_bf = Cst.nperiods_bf;
nperiods_af = Cst.nperiods_af;
onwership   = Data.ownership_rate;
T = Cst.T;


Time = Data.Time;
wodd = Data.wodd_BJ;



phi_it  = statei.BJ.phi_it;

inc = statei.inci;
incg = statei.incgi;
%own = theta2m(1,4);
own = 1/(1+exp(theta2m(1,4)));
sigma = theta2m(1,5);
kappas = theta2m(:,1);
%kappas(2:end) = kappas(1)+kappas(2:end);
%kappa_it = reshape(kappas(incg),nperiods,2);




s_ijt = zeros(nprods, nmarkets);
f_ijt = zeros(nprods, nmarkets);
l_it  = zeros(nperiods_af,1);
o_it = zeros(nperiods,nmarkets);
remain = [ones(1, nmarkets); zeros(nperiods, nmarkets)];
winner = zeros(nperiods+1, T+1); %Last column means given up lotteries

for mm=1:nmarkets
    if mm==1
     delta_it = statei.BJ.delta_it;
    else
     delta_it = statei.TJ.delta_it;
    end
delta_ijt = delta_it(Time,:); 
%Probability of buying vehicle j conditional on purchasing
cond_p_ij = exp(deltaf_ijt(:,mm)+rho_ijt(:,mm))./exp(delta_ijt);
for tt=1:nperiods
  remain_t = remain(tt,mm);
    
    % non-owner
    evifwd_t = evifwd(1, incg(tt,mm), tt,mm);
    p0_it     = exp(delta_it(tt)+kappas(incg(tt,mm)))./(exp(delta_it(tt)+kappas(incg(tt,mm)))+ exp(beta*evifwd_t)); 
    
    %owner
    evifwd_t = evifwd(2, incg(tt,mm),tt,mm);
    p1_it     = own*exp(delta_it(tt))./(exp(delta_it(tt))+ exp(beta*evifwd_t)); %(1-beta)*own+

    %Probability of the need for new license
    %initial ownership
    if mm==1
        temp = 0.68*exp(-2.348115+ 0.1033721*inc(tt,mm));
    else
        temp = 0.67*exp(-1.493984+ 0.1111152*inc(tt,mm));
    end
    if tt==1
    L_it = 1;%1/(1+temp);
    else
    L_it = 1-o_it(tt-1,mm);
    end

    if mm~=1 || tt<= nperiods_bf %|| own(mm)==1
    p_it = L_it*p0_it+ (1-1/(1+temp))*p1_it;
    s_ijt(Time==tt,mm) = s_ijt(Time==tt,mm)+ p_it*cond_p_ij(Time==tt);
    f_ijt(Time==tt,mm) = f_ijt(Time==tt,mm)+ L_it*p0_it*cond_p_ij(Time==tt);
   % remain(tt+1,mm) = remain(tt,mm).*(1-p_it);
     if tt==1
       o_it(tt,mm) = L_it*p0_it;%+ 1-L_it;
     else
       o_it(tt,mm) = o_it(tt-1,mm) + L_it *p0_it;
     end
    else
    
    %Purchases from a non-onwer whom does not need a new license
    s_ijt(Time==tt,mm) = s_ijt(Time==tt,mm)+ (1-1/(1+temp))*p1_it*cond_p_ij(Time==tt);
    
    %Lottery Participation
    wodd_t = wodd(tt-nperiods_bf);
    v = phi_it(tt-nperiods_bf)-wodd_t*beta*ewifwd(tt-nperiods_bf);
    l_it(tt-nperiods_bf)  = L_it*(1-1/(exp(v/sigma)+1));%exp(phi_it(tt-nperiods_bf)/sigma)/(1+exp(phi_it(tt-nperiods_bf)/sigma));
    winner(tt,1) = l_it(tt-nperiods_bf)*wodd_t; 
    %remain(tt+1,mm) = remain(tt,mm)-winner(tt+1,1);
    
    %Lottery Winner Purchasing
    ev0ifwd_t = ev0ifwd(:,tt);
    p0_it   = exp(delta_it(tt)+kappas(incg(tt,mm)))./(exp(delta_it(tt)+kappas(incg(tt,mm)))+ exp(beta*ev0ifwd_t));
    
    for ss=1:T
    f_ijt(Time==tt,mm) = f_ijt(Time==tt,mm)+ winner(tt,ss).*p0_it(ss).*cond_p_ij(Time==tt);
    winner(tt+1,ss+1) = winner(tt, ss)*(1-p0_it(ss)); 
    end
    s_ijt(Time==tt,mm) =  s_ijt(Time==tt,mm) +f_ijt(Time==tt,mm);
    o_it(tt,mm) = o_it(tt-1,mm)+sum(f_ijt(Time==tt,mm)); 
  end

end
end

end
