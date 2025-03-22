function HYBbs = doHYBbs(HYB,nboot,clevel)

[T,n]    = size(HYB.vars);

 VARSLAGS = lagmatrix(HYB.vars,1:HYB.p);
 TSHOCKS  = HYB.TSHOCKS(HYB.p+1:end,:);
 VARSLAGS = VARSLAGS(HYB.p+1:end,:);
 VARS     = HYB.vars(HYB.p+1:end,:);
 
 bet=[VARSLAGS TSHOCKS  HYB.DET(HYB.p+1:end,:)]\VARS; 
 
 res = VARS-[VARSLAGS TSHOCKS  HYB.DET(HYB.p+1:end,:)]*bet;
 
 
HYBBS = HYB;

 
  jj=1;
  while jj<nboot+1
      %rb   =  res(randsample(T-HYB.p,T-HYB.p,true),:);
      
      rr = 1-2*(rand(T-HYB.p,1)>0.5);
      rb = [rr.*res(:,1) rr.*res(:,2) rr.*res(:,3)];
      
      resb =  rb';
      
      varsb = zeros(size(HYB.vars));
      varsb(1:HYB.p,:)=HYB.vars(1:HYB.p,:);
     for j=HYB.p+1:T
     lvars = (varsb(j-1:-1:j-HYB.p,:))';
     varsb(j,:) = lvars(:)'*bet(1:HYB.p*n,:)+TSHOCKS(j-HYB.p,:)*bet(HYB.p*n+1:HYB.p*n+(HYB.q+1),:)...
         +HYB.DET(j,:)*bet(HYB.p*n+(HYB.q+1)+1:end,:)+resb(:,j-HYB.p)';     
     end

     HYBBS.vars = varsb;
     HYBBS = doHYB(HYBBS);
     
%     HYBBS.irs = transplot(HYBBS.irs,HYBBS.shocksize,TRY,GY);
     HYBbs.irs(:,jj) = HYBBS.irs(:);
     jj=jj+1;
  end
 
   HYBbs.irsH=reshape(quantile(HYBbs.irs',(1-clevel/100)/2),HYB.irhor,n);
   HYBbs.irsL=reshape(quantile(HYBbs.irs',1-(1-clevel/100)/2),HYB.irhor,n);
  
  
