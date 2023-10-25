function [se,rho,xfwd] = f_CalcPhiEvolution(x_it, Data, Cst, mid_phi)
nperiods_af = Cst.nperiods_af;
%nperiods_bf = Cst.nperiods_bf;

%mid_phi = Cst.mid_phi;
%nphi      = Cst.nphi;
%mid_delta_it = Cst.mid_delta_it;
ndelta  = Cst.ndelta_it;
%ntrans   = Cst.ntrans;
tran_draws = Data.tran_draws;



 
X = x_it(2:end);
lX = x_it(1:end-1);
 
YY = X;
%XX = [lX, delta_it(nperiods_bf+1:end-1), ones(nperiods_af-1,1)];
 XX = [lX, ones(length(lX),1)];
 rho   = pinv(XX'*XX)*XX'*YY;
 res     = YY-XX*rho;
 se    = sqrt(res'*res/(nperiods_af-size(XX,2)));
 %std   = sqrt(diag(se.^2*inv(XX'*XX)));
 %Next period phis
 %xfwd = mid_phi*rho(1)+rho(2)+se.*tran_draws;
xfwd = mid_phi*rho(1)+rho(2)+se.*tran_draws;



end
