function [moments, paras, Cst] = f_moments(theta2m, Data, Cst, Iter)

nprods = Cst.nprods;
nperiods = Cst.nperiods;
nperiods_af = Cst.nperiods_af;
nmarkets = Cst.nmarkets;
nyears   = Cst.nyears;
nmonths  = Cst.nmonths;   
%Moments

if Iter.flag_NAN

else
deltaf_jt = Iter.deltaf_jt;
eta_t     = Iter.eta_t;


% Moment  1
% 2SLS to get characteristics and price coefficients
% Find residual of the regression:
% If using the aggregated (by segment) data 
% deltaf_{sct} = X_{st}beta+ alpha P_{st} + xi_{sct}
% If using the model-level data
% deltaf_{jct} = X_{jt}beta+ alpha P_{jt} + xi_{jct}
% Our IV for price is tax discount dummy and consumer tax rate

Y = deltaf_jt(:);
X = Data.X_mom11;
W = [Data.X_mom11(:,1:end-1), Data.IV_mom11];
[b11, std] = tsls(Y,X,W);
paras.lin_paras = b11;
paras.std_lin_paras = std;

xi = Y - X*b11;

%Construct moments using the DID assumption
% Find residual of the regression:
% If using the aggregated (by segment) data 
% xi_{sct} = d_{s, yr} + d_t + trend_{c} + trend_{c}^2 + e_{sct}
% If using the model-level data
% xi_{jct} = d_{j, yr} + d_t + d_{cs} + trend_{c} + trend_{c}^2 + e_{jct}
% d_t is dummy for each (year-month)
% d_{s, yr} is dummy for each segment (s - year)
% trend_c is a log time trend for city c
% trend_c^2 is log time trend squared for city c
% Our IV for mom1 is city-year dummy, d_cy

Y = demean(xi, Data.model_vin);
X = Data.X_mom12;
b12 = OLS(Y,X);
e_mom1 = Y-X*b12;

if FlagSpecs.isNoLast2==1
%If drop last two months before lottery and a month after lottery
idx = Time == 35 | Time ==36 | Time==36;
idx = repmat(idx,2,1);
moments.mom1 = e_mom1(~idx).'*Data.IV_mom12(~idx,:)/(nprods*nmarkets-sum(idx));
else
moments.mom1 = e_mom1.'*Data.IV_mom12/(nprods*nmarkets);
end

if Cst.optw
moments.W_mom1 = Data.OW_mom1;
else
moments.W_mom1 = Data.IV_mom12'*Data.IV_mom12;
end
shat=e_mom1*ones(1,size(Data.IV_mom12,2)).*Data.IV_mom12;
moments.OW_mom1=shat'*shat/length(shat);


% Moment 2 - Lottery participation unobservables
%Regress eta_t on constant and policy change
Y = eta_t;
X = Data.X_mom2;
b21 = OLS(Y,X);
res = Y-X*b21;

%Fit AR(1) on residuals
Y = res(2:end);
X = res(1:end-1);
b22 = OLS(Y,X);

e_mom2 = Y-X*b22;
moments.mom2 = e_mom2.'*Data.IV_mom2/(nperiods_af-1);

if Cst.optw
moments.W_mom2 = Data.OW_mom2;
else
moments.W_mom2 = Data.IV_mom2'*Data.IV_mom2;
end
lhat=e_mom2*ones(1,size(Data.IV_mom2,2)).*Data.IV_mom2;
moments.OW_mom2=lhat'*lhat/length(lhat);



% Moment 3 - Micro-moments
incmom = Data.incmom;
incsegmom = Data.incsegmom;
first_share = Data.first_share;
replace_share = Data.replace_share;
owner_by_income = Data.owner_by_income;

share_ijt = Iter.share_ijt;
first_ijt = Iter.first_ijt;
replace_ijt = share_ijt-first_ijt;
parate_it = Iter.parate_it;
owner_it = Iter.owner_it;

%Matching purchase by income groups
incmom_hat = zeros(4, nyears, nmarkets);
incsegmom_hat = zeros(4,3, nmarkets);
first_sharehat = zeros(4,nyears,nmarkets);
replace_sharehat = zeros(nyears,nmarkets);
owner_by_incomehat = zeros(4,nyears,nmarkets);

for mm=1:nmarkets
share_it = Data.TimeProduct*squeeze(share_ijt(:,mm,:));
first_it = Data.TimeProduct*squeeze(first_ijt(:,mm,:));
owner    = squeeze(owner_it(:,mm,:));
replace_it = Data.TimeProduct*squeeze(replace_ijt(:,mm,:));
tshare   = sum(share_it,2);

for g=1:4

for y=1:nyears
ind = squeeze(Data.inc_group_t(12*y,mm,:))==g;
temp1 = sum(share_it(:,ind),2);
temp2 = sum(first_it(:,ind),2);
temp3 = mean(owner(:,ind),2);
temp4 = sum(squeeze(share_ijt(:,mm,ind)),2);   
incmom_hat(g,y,mm) = sum(temp1((y-1)*nmonths+1:y*nmonths))./sum(tshare((y-1)*nmonths+1:y*nmonths));
first_sharehat(g,y,mm) = sum(temp2((y-1)*nmonths+1:y*nmonths))./sum(temp1((y-1)*nmonths+1:y*nmonths));
owner_by_incomehat(g,y,mm) = temp3(y*nmonths-2);

if y>1
for s=1:3
indseg1 = Data.seg_group == s & Data.year>2008;
indseg2 = Data.seg_group == s & Data.year==2007+y;
temp5 = sum(squeeze(share_ijt(indseg1,mm,:)),2);
incsegmom_hat(g,s,mm) = incsegmom_hat(g,s,mm)+sum(temp4(indseg2))./sum(temp5);
end
end

end

end

temp3 = sum(replace_it,2);
for y=1:nyears
replace_sharehat(y,mm) = sum(temp3((y-1)*nmonths+1:y*nmonths))./sum(tshare((y-1)*nmonths+1:y*nmonths));
end
end
end

e_mom31 = incmom_hat(2:end,2:end,:)-incmom;
e_mom32 = first_sharehat(:,2:end,:) - first_share;
e_mom33 = replace_sharehat - replace_share;
e_mom34 = squeeze(owner_by_incomehat(:,3,:))-owner_by_income;
e_mom35 = incsegmom_hat(2:end,:,:)-incsegmom;

moments.mom3 = [e_mom31(:); e_mom32(:); e_mom33(:); e_mom34(:); e_mom35(:)];
if Cst.optw==1
moments.W_mom3 = blkdiag(Data.OW_mom31,Data.OW_mom32, Data.OW_mom33, Data.OW_mom34, Data.OW_mom35);
else
moments.W_mom3 = blkdiag(Data.W_mom31, Data.W_mom32, Data.W_mom33, Data.W_mom34, Data.W_mom35);
end
save 'Interim/temp.mat'
end
