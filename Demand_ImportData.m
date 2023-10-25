%============================================================
    %IMPORT DATA
%============================================================
rng(12345,'twister');
%Data
data_temp = dataset('file','sales_attributes_new.txt','ReadVarNames',false);
%Income draws
load(strcat('inc_draws_c05_',num2str(Cst.nconsumers),'_nmonths_',num2str(Cst.nmonths)));
% Population data
load('nhouseholds');
% CPI
load('cpi');
% Income data (for micro moments)
load('shareincome');
samplesize_inc = samplesize;
load('firstincome');
samplesize_first = sample_size;
load('replaceshare');
samplesize_replace = sample_size;
load('owner_by_income');
samplesize_owner = sample_size;
load('incsegmom');
samplesize_incseg = sample_size;

% Quota data
load('beijing_quota') % quota by month from Jan 2011 to Dec 2012
% Pentdemand load
load('beijing_pentdemand')
beijing_pentdemand=beijing_pentdemand/10000;
load('beijing_drop')
beijing_drop=beijing_drop/10000;

% Column Sorter.  Car characteristics at end.
%Observation data
col.brand       = 1;
col.model       = 2;
col.model_vin   = 3; % Unique identifying for each (model-year pair).
col.time        = 4;
col.year        = 5;
col.month       = 6;
%Quantity and price 
col.shareB      = 7;
col.shareT      = 8;
col.shareS      = 9;
col.shareN      = 10;
col.priceS      = 11; 
col.priceB       = 12;
col.priceT      = 13;
col.priceN      = 14;

% Car Characteristics
col.segment     = 15;
col.type        = 16;
col.hybrid      = 17;
col.electric    =18;
col.imported    =19;
col.vehlength   =20;
col.vehwidth    =21;
col.vehdisp     =22;
col.vehweight   =23;
col.doors       =24;
col.seats       =25;
col.engine_ccm  =26;
col.horsepower  =27;
col.urban_1pm   =28;
col.suburban_1pm=29;
col.manual      = 30;
col.vehsize     = 31;
col.vehlpm      = 32;
col.fc          = 33;
col.collapse_ind = 34;
col.collapse_month = 35;
col.constant    = 36;

% Characteristic Indicator
char_ind = [col.vehlpm,col.vehsize, col.vehdisp];

%% Data Matrix
data = zeros(size(data_temp,1),col.constant);
% Observation = (Model,Time,Year)
% data(:,col.brand)  - change to numerical value below
data(:,col.model)        =  double(data_temp(:,28));
% data(:,col.model_vin)  - calculate below using model and year
data(:,col.time)         =  double(data_temp(:,9)); 
data(:,col.year)         =  double(data_temp(:,11));  
data(:,col.month)        =  double(data_temp(:,12)); 
% Sales and price data
data(:,col.shareB)       =  double(data_temp(:,13)); %which is beijing? % Sales data - adjust for shares %all shares adjusted by the first purchase only after 2008
sales_beijing   =  double(data_temp(:,13)); % Save for later
data(:,col.shareS)       =  double(data_temp(:,15)); % Sales data - adjust for shares
data(:,col.shareN)       =  double(data_temp(:,14)); % Sales data - adjust for shares
data(:,col.shareT)       =  double(data_temp(:,10)); % Sales data - adjust for shares
%data(:,col.priceS)      % Load Shanghai prices shortly
data(:,col.priceB)        =  double(data_temp(:,24));
data(:,col.priceT)        =  double(data_temp(:,24));
%Car characteristics
%data(:,col.segment) This is string.  Convert to number below.
data(:,col.type)         =  double(data_temp(:,4));
data(:,col.hybrid)       =  double(data_temp(:,6));
data(:,col.electric)     =  double(data_temp(:,7));
data(:,col.imported)     =  double(data_temp(:,8));
data(:,col.vehlength)    =  double(data_temp(:,16));
data(:,col.vehwidth)     =  double(data_temp(:,17));
data(:,col.vehdisp)      =  double(data_temp(:,18)); % Coding error?
data(:,col.vehweight)    =  double(data_temp(:,19));
data(:,col.doors)        =  double(data_temp(:,20));
data(:,col.seats)        =  double(data_temp(:,21));
data(:,col.engine_ccm)   =  double(data_temp(:,22));
data(:,col.horsepower)   =  double(data_temp(:,23)); 
data(:,col.urban_1pm)    =  double(data_temp(:,25));
data(:,col.suburban_1pm) =  double(data_temp(:,26));
data(:,col.manual)       =  double(data_temp(:,27));
data(:,col.constant)     =  ones(size(data(:,1)));


%% Correct Data Matrix 
% Change brand to numerical value
% Coded from 1 to 90 different brands.
brand = cellstr(data_temp(:,2));
uniq_brand = unique(brand);
for ii = 1:size(uniq_brand,1)
    data(strcmpi(brand,uniq_brand{ii}),col.brand) = ii;
end

% Change segment to numerical value
% 1 - large, 2 - luxury, 3 - medium, 4 - mini, 5 - small, 6 - upper medium.
segment = cellstr(data_temp(:,5));
uniq_seg = unique(segment);
for ii = 1:size(uniq_seg,1)
    data(strcmpi(segment,uniq_seg{ii}),col.segment) = ii;
end


% Vehicle size
data(:,col.vehsize) = data(:,col.vehlength).* data(:,col.vehwidth);

% Vehicle efficiency (urban and suburban)
data(:,col.vehlpm)  = (data(:,col.urban_1pm)+data(:,col.suburban_1pm))/2;
clear data_temp

%% Model-Vintage 
% model_vin is a unique identifyer for each (model,year)-pair.  This is
% used as part of the calculation for the first moment condition.
uniq_model = unique(data(:,col.model));
uniq_year = unique(data(:,col.year));
nmodel_vin = size(uniq_model,1)*size(uniq_year,1);
id = 1;
for mm = 1:size(uniq_model,1)
   for yy = 1:size(uniq_year,1)
       % Find all observations with model mm and year yy
       ind = find(data(:,col.model) == uniq_model(mm) & data(:,col.year) == uniq_year(yy));
       if isempty(ind) == 1
           continue;
       end
       data(ind,col.model_vin) = id;
       id = id+1;
   end
end
clear uniq_model uniq_year nmodel_vin
%Number of Products
Cst.nprods = size(data,1);
%% Market Size
% Assume three quater of households are on the market each year
nhouseholds=nhouseholds*10000*.75;

%% Market Shares
% If sales == 0, drop it.
data(data(:,col.shareB) ==0,:) = [];
data(data(:,col.shareN) ==0,:) = [];
data(data(:,col.shareS) ==0,:) = [];
data(data(:,col.shareT) ==0,:) = [];
sales=zeros(Cst.nyears*12,Cst.nmarkets);
id=1;
for ii=1:Cst.nyears
  for mm=1:12
  ind=find(data(:,col.year)==(2007+ii) & data(:,col.month)==mm);
  if isempty(ind)
      continue
  else
  sales(id,1)=sum(data(ind,col.shareB));
  sales(id,2)=sum(data(ind,col.shareT));
  sales(id,3)=sum(data(ind,col.shareS));
  sales(id,4)=sum(data(ind,col.shareN));
  id=id+1;
  end
  end
end
 nhouseholds = reshape(nhouseholds,Cst.nyears,4);
% population in beijing,tianjin,shanghai,nanjin in 2008 to 2012
% Columns are beijing,tianjin,shanghai,nanjin
% Divide sales by market size to get shares
for tt=1:Cst.nyears
    for mm=1:12
    % Find entries with year tt and month mm
    ind=find(data(:,col.year)==Cst.year0+tt & data(:,col.month)==mm);
        for jj=1:4
            data(ind,col.shareB + jj - 1) = data(ind,col.shareB + jj - 1)/( nhouseholds(tt,jj));
        end
    end
end

%% Prices and income
% Adjust price so it's in 1000's of yaun
data(:,col.priceB) =(data(:,col.priceB)/1000);
data(:,col.priceT) =(data(:,col.priceT)/1000);
%Income draws in 10000 of yaun

inc_draws=(inc_draws/10000);

% Adjust price for inflation
for tt = 1:Cst.nyears
    % Find entries with year == tt
    ind = find(data(:,col.year) == Cst.year0+tt);
    data(ind,col.priceB) =  data(ind,col.priceB)*100/cpi(tt);
    data(ind,col.priceT) =  data(ind,col.priceT)*100/cpi(tt);
end

% Adjust income draws for inflation
for t=1:Cst.nyears
   inc_draws(:,(t-1)*Cst.nmonths+1:t*Cst.nmonths,:)=100*inc_draws(:,(t-1)*Cst.nmonths+1:t*Cst.nmonths,:)/cpi(t);
end

% Adjust price for sales tax
%consumption tax rate
CTRate = zeros(size(data,1),1);
disp = data(:,col.vehdisp);
vintage = data(:,col.year);

ind = disp<=1;
CTRate(ind,1) = 0.01;
ind = disp<=1.5 & disp>1;
CTRate(ind,1) = 0.03;
ind = disp<=2 & disp>1.5;
CTRate(ind,1) = 0.05;
ind = disp<=2.5 & disp>2;
CTRate(ind,1) = 0.09;
ind = disp<=3 & disp>2.5;
CTRate(ind,1) = 0.12;
ind = disp<=4 & disp>3;
CTRate(ind,1) = 0.25;
ind = disp>4;
CTRate(ind,1) = 0.40;

STRate = zeros(size(data,1),1)+0.1;%standard 10% 
% for 1.6L and below, the sales tax rate is 5% in 2009 and 7.5% in 2010

ind = disp<=1.6 & vintage == 2009;
STRate(ind,1) = 0.05;
ind = disp<=1.6 & vintage == 2010;
STRate(ind,1) = 0.075;
taxdis = 0.1-STRate;

data(:,col.priceB) = data(:,col.priceB)+data(:,col.priceB).*(1-CTRate).*STRate;
data(:,col.priceT) = data(:,col.priceT)+data(:,col.priceT).*(1-CTRate).*STRate;
clear ind tt m


%% Representative group
if FlagSpecs.isGroupSeg
uniq_year = unique(data(:,col.year));
%segment: 1 - large, 2 - luxury, 3 - medium, 4 - mini, 5 - small, 6 - upper medium.
%type: 0 -sedan, 1-suv, 2-minivan

uniq_group  = unique(data(:,[col.imported,col.type,col.segment]),'rows');
ngroup      = length(uniq_group);
CTRate2     = [];
taxdis2   = [];
seg_group =[];
tt        = 1;
for yy=1:Cst.nyears
    for mm=1:12
        for ss=1:ngroup
%     ind = find(data(:,col.year) == uniq_year(yy) & data(:,col.month)==mm);
    ind = find(data(:,col.year) == uniq_year(yy) & data(:,col.month)==mm ...
        & data(:,col.imported)==uniq_group(ss,1) &data(:,col.type)==uniq_group(ss,2) & data(:,col.segment)==uniq_group(ss,3));
    if ~isempty(ind)
    data2(tt,:) = mean(data(ind,:),1);
    data2(tt,col.shareB)=sum(data(ind,col.shareB),1);
    data2(tt,col.shareS)=sum(data(ind,col.shareS),1);
    data2(tt,col.shareN)=sum(data(ind,col.shareN),1);
    data2(tt,col.shareT)=sum(data(ind,col.shareT),1);
    data2(tt,col.model) = ss;
    data2(tt,col.model_vin) = (yy-1)*ngroup+ss;
    data2(tt,col.type) = mode(data(ind,col.type),1);
    data2(tt,col.hybrid) = mode(data(ind,col.hybrid),1);
    data2(tt,col.electric) = mode(data(ind,col.electric),1);
    data2(tt,col.imported) = uniq_group(ss,1);
    data2(tt,col.doors) = mode(data(ind,col.doors),1);
    data2(tt,col.seats) = mode(data(ind,col.seats),1);
    data2(tt,col.manual) = mean(data(ind,col.manual),1);
    data2(tt,col.year)  = 2007+yy;
    data2(tt,col.month) = mean(data(ind,col.month),1);
    data2(tt,col.priceB) = mean(data(ind,col.priceB),1);
    data2(tt,col.priceT) = mean(data(ind,col.priceT),1);
    data2(tt,col.time)  = (yy-1)*12+mm;
    CTRate2 = [CTRate2; mean(CTRate(ind))];
    taxdis2 = [taxdis2; mean(taxdis(ind))];
    
% A coarse segment classification combining type and segment
% 1 - compact sedan (mini, small), 2 - medium/ upper medium sedan, large and luxury sedan, 3 - SUV/ Minivan
 if uniq_group(ss,2) == 0 && (uniq_group(ss,3) <=3 || uniq_group(ss,3) ==6)
 seg_group = [seg_group; 2];
 elseif uniq_group(ss,2) == 0 && (uniq_group(ss,3) ==4 || uniq_group(ss,3) ==5)
 seg_group = [seg_group; 1];
 elseif uniq_group(ss,2) ==  1 || uniq_group(ss,2) ==2
 seg_group = [seg_group; 3]; 
 end

    tt=tt+1;
    end
        end
    end
end
data=data2;
taxdis = taxdis2;
CTRate = CTRate2;
%Number of Products
Cst.nprods = size(data,1);
end
 Data.seg_group = seg_group;
%% Lottery data
% Get winning odd
  Data.wodd_BJ          = beijing_quota./beijing_pentdemand/10000;

% Get predicted participation rates for Beijing 
  Data.parate_bj   = beijing_pentdemand*10000./kron(nhouseholds(4:5,1),ones(Cst.nmonths,1));
% Get dropped rate for Beijing
  Data.droprate_bj = beijing_drop./beijing_quota;
            

%% Price, Share, and Characteristics
% price = [price of all other cities, price of shanghai
priceB = data(:,col.priceB);
priceT = data(:,col.priceT);
characteristics = [data(:,col.constant),log(data(:,char_ind))];
Data.share = [data(:,col.shareB),data(:,col.shareT),data(:,col.shareS),data(:,col.shareN)];
Data.model_vin = repmat(data(:,col.model_vin),Cst.nmarkets,1);
Data.year = data(:,col.year);


%% Moment Variables
% Moment 1 - vehicle demand
Data.model           = repmat(data(:,col.model),Cst.nmarkets,1);
d_time          = repmat(dummyvar(data(:,col.time)),Cst.nmarkets,1);
d_time(:,1:Cst.nmonths:Cst.nmonths*Cst.nyears) = []; % drop first month in every year to avoid multi-collinearity

% Use log(2) for 2008 to distinguish between first year and when city = 0.
trend_c         = kron(eye(Cst.nmarkets,Cst.nmarkets),log(data(:,col.year)-2007+1));
trend_c2        = kron(eye(Cst.nmarkets,Cst.nmarkets),log(data(:,col.year)-2007+1).^2);

% Drop first column of trend_c and trend_c2 to avoid multicollinearity
trend_c(:,1) = [];
trend_c2(:,1) = [];

d_cy            = dummyvar(data(:,col.year)-2007);
d_cy(:,1) = [];  % drop first column to avoid multi-collinearity
d_cy            = kron(eye(Cst.nmarkets),d_cy);
d_cy(:,1:Cst.nyears - 1) = []; % drop first four columns so that the demeaned (on model vintage) d_cy has full column rank.

d_cs            = dummyvar(data(:,col.segment));
d_cs(:,1) = [];   % drop first column to avoid multi-collinearity
d_cs            = kron(eye(Cst.nmarkets,Cst.nmarkets),d_cs);
d_cs(:,1:5) = []; % drop segments from first city to get full rank


% De-mean variables WRT model_vintage to remove d_j
Data.price = [priceB,priceT];
Data.chars    = repmat(characteristics,1,1,Cst.nmarkets);
Data.X_mom11 = [repmat(characteristics,Cst.nmarkets,1), log(Data.price(:))];
if FlagSpecs.isGroupSeg
Data.X_mom12 = demean([d_time trend_c trend_c2],Data.model_vin);
else
Data.X_mom12 = demean([d_time d_cs trend_c trend_c2],Data.model_vin);  
end
Data.IV_mom11 = repmat([CTRate,taxdis],2,1);
Data.IV_mom12 = [d_cy];

% Moment 2 - lottery particiaption
d_policychange = [zeros(15,1); ones(Cst.nperiods_af-15,1)]; %Policy change for participation
Data.X_mom2 = [ones(Cst.nperiods_af,1), d_policychange];

ave_chars = grpstats(characteristics(:,2:end), data(:,col.time), 'mean');
ave_prices = grpstats(priceB, data(:,col.time), 'mean');

Data.IV_mom2 = [Data.wodd_BJ(1:end-1)*100, (Data.wodd_BJ(2:end)-Data.wodd_BJ(1:end-1))*100, ...
                beijing_quota(1:end-1)/1e4, (beijing_quota(2:end)-beijing_quota(1:end-1))/1e3, ave_chars(Cst.nperiods_bf+1:end-1,:)/10, ave_prices(Cst.nperiods_bf+1:end-1)/100];


% Moment 3 - micro moments
incmom_BJ = reshape(incmom(1:4*(Cst.nyears-1)),4, Cst.nyears-1);
incmom_TJ = reshape(incmom(4*(Cst.nyears-1)+1:end),4, Cst.nyears-1);
%Remove on income group
incmom_BJ(1,:)= [];
incmom_TJ(1,:)= [];
incmom(1:4:end)=[];
samplesize_inc(1:4:end)=[];
Data.incmom = cat(3, incmom_BJ, incmom_TJ); 

incsegmom_BJ = reshape(incsegmom(1:4*3),4,3);
incsegmom_TJ = reshape(incsegmom(4*3+1:end),4,3);
incsegmom_BJ(1,:) = [];
incsegmom_TJ(1,:) = [];
incsegmom(1:4:end) = [];
samplesize_incseg(1:4:end) = [];
Data.incsegmom = cat(3,incsegmom_BJ, incsegmom_TJ);

%first_share moments start from 2009
first_shareBJ = reshape(first_by_income(1:4*(Cst.nyears-1),1),4,Cst.nyears-1);
first_shareTJ = reshape(first_by_income(4*(Cst.nyears-1)+1:end,1),4,Cst.nyears-1);
Data.first_share = cat(3, first_shareBJ, first_shareTJ);
Data.replace_share = reshape(share_replace,Cst.nyears,2); 
Data.owner_by_income = reshape(owner_by_income,4,2); %by the end of 2012, ownership by hh income

%weighting matrix for micro moments
temp = blkdiag(1-eye((Cst.nyears-1)*3),1-eye((Cst.nyears-1)*3));
Data.W_mom31 = diag(1./samplesize_inc);%diag(incmom.*(1-incmom))./samplesize_inc; 
Data.W_mom32 = diag(1./samplesize_first);%diag(first_by_income.*(1+1e-10-first_by_income)./samplesize_first);% 
Data.W_mom33 = diag(1./samplesize_replace);%diag(share_replace.*(1-share_replace)./samplesize_replace);%
Data.W_mom34 = diag(1./samplesize_owner);%diag(owner_by_income.*(1-owner_by_income)./samplesize_owner);%
Data.W_mom35 = diag(1./samplesize_incseg);
Data.OW_mom31 = diag(incmom.*(1-incmom))./samplesize_inc;
Data.OW_mom32 = diag(first_by_income.*(1+1e-10-first_by_income)./samplesize_first);
Data.OW_mom33 = diag(share_replace.*(1-share_replace)./samplesize_replace);
Data.OW_mom34 = diag(owner_by_income.*(1-owner_by_income)./samplesize_owner);
Data.OW_mom35 = diag(incsegmom.*(1-incsegmom))./samplesize_incseg;
Data.obs3  = [samplesize_inc; samplesize_replace;samplesize_owner;samplesize_incseg];
clear incmom yy ii mm d_cy trend_c trend_c2 d_time

%% Evolution Simulation
%Monthly Dummy
Data.d_monthly        = dummyvar(repmat((1:Cst.nmonths)',Cst.nyears,1));
Data.month            = repmat((1:Cst.nmonths)',Cst.nyears,1);
%% Income and House Draws
% Income draws is generating in gen_income_census05
% Halton Draws
house_draws=zeros(Cst.nconsumers,Cst.nrand,Cst.nperiods,Cst.nmarkets);  %60 months, 4 cities, 5 random characteristics: constant, size, displacement, weight, price
halt=haltonseq(Cst.nconsumers*Cst.nmarkets+10,Cst.nrand)+rand(Cst.nconsumers*Cst.nmarkets+10,Cst.nrand);
temp=find(halt>1);
halt(temp)=halt(temp)-1;
draws=norminv(halt,0,1);

draws(1:10,:)=[];

for t=1:Cst.nperiods
    for j=1:Cst.nmarkets
        for k=1:Cst.nrand
            house_draws(:,k, t,j)=draws((j-1)*Cst.nconsumers+1:j*Cst.nconsumers, k);
        end
    end
%     draws(1:Cst.nmarkets*Cst.nconsumers,:)=[];
end
clear draws halt temp

%Initial Ownership draws 
own_B = 207.94/565.9;
own_T = 66.61/333.42;
own = [own_B;own_T];

ownership=[
207.94/565.9  66.61/333.42;
244.27/605.9864  79.89/347.8753;
296.56/636.2857  100.01/356.9201;
371.51/668.1000  125.70/366.2000;
387.29/687.8602  155.46/383.3497;
]; 


Data.ownership_rate = kron(ownership, ones(Cst.nmonths,1));

halt = haltonseq(Cst.nconsumers*Cst.nmarkets+10,1)+rand(Cst.nconsumers*Cst.nmarkets+10,1);
temp=find(halt>1);
halt(temp)=halt(temp)-1;
halt(1:10) = [];
owner_draws = zeros(Cst.nconsumers,Cst.nmarkets);

for j=1:Cst.nmarkets
draws= halt((j-1)*Cst.nconsumers+1:j*Cst.nconsumers)<=own(j);
owner_draws(:,j) = draws;
end


% Reshape inc_draws and house_draws -> required for calculating utility
inc_draws_temp   = zeros(Cst.nprods,Cst.nmarkets,Cst.nconsumers);
inc_draws_t      = zeros(Cst.nperiods, Cst.nmarkets, Cst.nconsumers);
std_inc_draws_temp   = zeros(Cst.nprods,Cst.nmarkets,Cst.nconsumers);
inc_group_temp = zeros(Cst.nprods,Cst.nmarkets,Cst.nconsumers);
inc_group_t    = zeros(Cst.nperiods, Cst.nmarkets, Cst.nconsumers);
house_draws_temp = zeros(Cst.nprods,Cst.nrand,Cst.nmarkets,Cst.nconsumers);

for yy = 1:Cst.nyears
   for month = 1:Cst.nmonths
       % Find indices with (yy,month)
       ind = find(data(:,col.year) == 2007+yy & data(:,col.month) == month);
       % Save indices with particular (year,month) for later
       ind_ym{(yy-1)*Cst.nmonths+month} = ind;
      
       for mm = 1:Cst.nmarkets
           ind1 = inc_draws(:,(yy-1)*Cst.nmonths+month,mm)<=4.8;
           inc_group_temp(ind,mm,ind1') = 1;
           inc_group_t((yy-1)*Cst.nmonths+month, mm, ind1') = 1;
           ind2 = (inc_draws(:,(yy-1)*Cst.nmonths+month,mm)<=9.6)&~ind1;
           inc_group_temp(ind,mm,ind2') = 2;
           inc_group_t((yy-1)*Cst.nmonths+month, mm, ind2') = 2;
           ind3 = (inc_draws(:,(yy-1)*Cst.nmonths+month,mm)<=14.4)&~ind2&~ind1;
           inc_group_temp(ind,mm,ind3') = 3;
           inc_group_t((yy-1)*Cst.nmonths+month, mm, ind3') = 3;
           ind4 = inc_draws(:,(yy-1)*Cst.nmonths+month,mm)>14.4;
           inc_group_temp(ind,mm,ind4') = 4;
           inc_group_t((yy-1)*Cst.nmonths+month, mm, ind4') = 4;
           
           %income draws are standardized
           temp = inc_draws(:,(yy-1)*Cst.nmonths+month,mm);
           inc_draws_temp(ind,mm,:) = repmat(temp,1,size(ind,1)).';
           inc_draws_t((yy-1)*Cst.nmonths+month,mm,:) = reshape(temp,1,1,Cst.nconsumers);
           temp = (inc_draws(:,(yy-1)*Cst.nmonths+month,mm)-mean(inc_draws(:,(yy-1)*Cst.nmonths+month,mm)))./std(inc_draws(:,(yy-1)*Cst.nmonths+month,mm));
           std_inc_draws_temp(ind,mm,:) = repmat(temp,1,size(ind,1)).';
           
           for ii = 1:Cst.nrand
               house_draws_temp(ind,ii,mm,:) = repmat(house_draws(:,ii,(yy-1)*Cst.nmonths+month,mm),1,size(ind,1)).';
           end
       end
   end
end

%Income group transition matrix
%Monthly transition matrix
inc_transBJ = [ 0.97890 0.02002 0.00022 0.00094
                0.00308 0.97246 0.02109 0.00354
                0.00358 0.01306 0.95854 0.02864
                0.00185 0.00370 0.00101 0.98196];
            
inc_transTJ = [0.96035 0.03580 0.00376 0.00009
               0.01957 0.96739 0.01272 0.00032
               0.00661 0.00072 0.94457 0.04810
               0.00031 0.02372 0.02046 0.95551];

           
Data.inc_trans = cat(3, inc_transBJ, inc_transTJ);
Data.inc_group = inc_group_temp;
Data.inc_group_t = inc_group_t;
Data.std_inc_draws = std_inc_draws_temp;
Data.inc_draws = inc_draws_temp;
Data.inc_draws_t = inc_draws_t;
Data.house_draws = house_draws_temp;
Data.owner_draws = owner_draws;
% Trans Draws used for computing transition matrix.
Data.tran_draws  = norminv(linspace(.5/Cst.ntrans,1-.5/Cst.ntrans,Cst.ntrans),0,1);
clear inc_draws_temp house_draws_temp mm yy month halt temp

%% Time, Time Product
% Time tracks time period that product j is in
% TimeProduct is used to aggregate inclusive value
Time = zeros(Cst.nprods,1);
Data.TimeProduct = zeros(Cst.nyears*Cst.nmonths,Cst.nprods);
for tt = 1:Cst.nyears*Cst.nmonths
   Time(ind_ym{tt},1) = tt*ones(size(ind_ym{tt}));
   Data.TimeProduct(tt,ind_ym{tt}) = ones(1,size(ind_ym{tt},1));
end

Data.Time = Time;

