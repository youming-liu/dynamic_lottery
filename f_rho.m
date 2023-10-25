function [f,f_ijt] = f_rho(theta2m,Data,Cst)
rho_ijt = zeros(Cst.nprods,Cst.nmarkets,Cst.nconsumers);
rho = zeros(Cst.nprods,Cst.nmarkets,Cst.nincgroups, Cst.nconsumers);
for mm=1:Cst.nmarkets
x = log(Data.price(:,mm));
v_i = squeeze(Data.house_draws(:,1,mm,:));
%d_i = Data.inc_draws(:,;,mm);
%d_i = Data.std_inc_draws(:,:,mm);
d_i = Data.inc_group(:,mm,:);

%rho(:,:,mm)=x.*(d_i*theta2m(1,2)+v_i*theta2m(1,3));  
idx = reshape(d_i,[],1);
alphas = -theta2m(:,2);
alpha0 = -10;
alphas = alpha0-alphas;
%alphas(1) = -10;
%alphas(2:end) = alphas(1)-alphas(2:end);
alpha_inc = reshape(alphas(idx),Cst.nprods,Cst.nconsumers);
rho_ijt(:,mm,:)= x.*(alpha_inc+v_i*theta2m(1,3));
for g=1:Cst.nincgroups
alpha_g = alphas(g);
rho(:,mm,g,:)= x.*(alpha_g+v_i*theta2m(1,3));
end
end
f=rho;
f_ijt=rho_ijt;
