function HYB = doHYB(HYB)

 [T,n]    = size(HYB.vars);

 VARSLAGS = lagmatrix(HYB.vars,1:HYB.p);
 TSHOCKS  = HYB.TSHOCKS(HYB.p+1:end,:);
 VARSLAGS = VARSLAGS(HYB.p+1:end,:);
 VARS     = HYB.vars(HYB.p+1:end,:);
 
 bet=[VARSLAGS TSHOCKS  HYB.DET(HYB.p+1:end,:)]\VARS; 
 
 irs(HYB.p+1,:) = bet(HYB.p*n+1,:);
 
 for j=2:HYB.irhor
 lvars = (irs(HYB.p+j-1:-1:j,:))';
 irs(HYB.p+j,:) = lvars(:)'*bet(1:HYB.p*n,:);     
  if j<=1+HYB.q   
     irs(HYB.p+j,:) = irs(HYB.p+j,:)+bet(HYB.p*n+j,:);
  end
 end
 HYB.irs = irs(HYB.p+1:end,:)*HYB.shocksize;