function [se,rho,xfwd] = f_CalcDeltaEvolution(x_it, Data, Cst, mid_delta_it)


  %Estimate Delta's Evolution Process
     X   = x_it(2:end);
     lX  = x_it(1:end-1);
     d_monthly = Data.d_monthly;
     nmonths = Cst.nmonths;
     
     XX  = [lX,d_monthly(1:end-1,:)];

     YY = X;
     
     ind = [35:37];
     
     if Cst.FlagSpecs.isNoLast2 ==1 %if removing last two months of pre-lottery and first-month of lottery
     XX(ind,:) = [];
     YY(ind,:) = [];
     end
     rho  = (XX'*XX)\XX'*YY;
     res  = YY-XX*rho;
     se   = sqrt(res'*res/(Cst.nperiods-size(XX,2)));  
     
     %Next period deltas
     mid_x_it =mid_delta_it;
     nx  = Cst.ndelta_it;
     ntrans   = Cst.ntrans;
     tran_draws = Data.tran_draws;
     
     xfwd = zeros(nx,ntrans,nmonths);
     for mm=1:nmonths
     xfwd(:,:,mm) = repmat(mid_x_it,1,ntrans)*rho(1)+rho(mm+1)+se*repmat(tran_draws,nx,1);
     end
end
