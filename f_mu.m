function f = f_mu(theta2m,Data,Cst)


k = size(Data.chars,2);
mu = zeros(Cst.nprods,Cst.nmarkets,Cst.nconsumers);

for mm=1:Cst.nmarkets
for kk=1:k
    x = squeeze(Data.chars(:,kk,mm));
    v_i = squeeze(Data.house_draws(:,kk+1,mm,:));
   mu(:,mm,:)=squeeze(mu(:,mm,:))+x.*v_i*theta2m(kk+1,3);
end
end
f=mu;
