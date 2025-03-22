function VARbs = doVARbs(VAR,nboot,clevel)

%[T,n]    = size(VAR.vars);

 jj=1;
  while jj<nboot+1      
     rr = 1-2*(rand(VAR.T,1)>0.5);
     eb     =  VAR.e.*(rr*ones(1,VAR.n));
     resb   =  VAR.D*eb';
     
     varsb = zeros(VAR.T,VAR.n);
     varsb(1:VAR.p,:)=VAR.vars(1:VAR.p,:);

     for j=VAR.p+1:VAR.p+VAR.T
        lvars = (varsb(j-1:-1:j-VAR.p,:))';
        varsb(j,:) = lvars(:)'*VAR.bet(1:VAR.p*VAR.n,:)+resb(:,j-VAR.p)';     %VAR.DET(j,:)*VAR.bet(VAR.p*VAR.n+1:end,:)+
     end
     
     VARBS = VAR;
     VARBS.vars = varsb;
     VARBS.d0 = VAR.d0;
     VARBS = doVAR(VARBS);

     if VARBS.rc==0;
     VARbs.irs(:,jj) = VARBS.irs(:);
     VARbs.irsg(:,jj) = VARBS.irsg(:);
     VARbs.thetaG(jj) = VARBS.thetaG;
     VARbs.zetaT(jj) = VARBS.zetaT;
     VARbs.zetaG(jj) = VARBS.zetaG;
     VARbs.sigmaT(jj) = VARBS.sigmaT;
     VARbs.sigmaG(jj) = VARBS.sigmaG;
     VARbs.sigmaY(jj) = VARBS.sigmaY;
     jj=jj+1;
     end
  end
  
  VARbs.irsH=reshape(quantile(VARbs.irs',(1-clevel/100)/2),VAR.irhor,VAR.n);
  VARbs.irsL=reshape(quantile(VARbs.irs',1-(1-clevel/100)/2),VAR.irhor,VAR.n);
  VARbs.irsgH=reshape(quantile(VARbs.irsg',(1-clevel/100)/2),VAR.irhor,VAR.n);
  VARbs.irsgL=reshape(quantile(VARbs.irsg',1-(1-clevel/100)/2),VAR.irhor,VAR.n);
  
  VARbs.thetaGci = [quantile(VARbs.thetaG',(1-clevel/100)/2) quantile(VARbs.thetaG',1-(1-clevel/100)/2)];
  VARbs.zetaTci  = [quantile(VARbs.zetaT',(1-clevel/100)/2) quantile(VARbs.zetaT',1-(1-clevel/100)/2)];
  VARbs.zetaGci  = [quantile(VARbs.zetaG',(1-clevel/100)/2) quantile(VARbs.zetaG',1-(1-clevel/100)/2)];
  VARbs.sigmaTci = [quantile(VARbs.sigmaT',(1-clevel/100)/2) quantile(VARbs.sigmaT',1-(1-clevel/100)/2)];
  VARbs.sigmaGci = [quantile(VARbs.sigmaG',(1-clevel/100)/2) quantile(VARbs.sigmaG',1-(1-clevel/100)/2)];
  VARbs.sigmaYci = [quantile(VARbs.sigmaY',(1-clevel/100)/2) quantile(VARbs.sigmaY',1-(1-clevel/100)/2)];

  
function [z,D] = Dsolve(x,thetaY,gammaT,gammaY,Sigma)
[rows,cols]=size(x);

n=1;
while n<=cols; 
thetaG = x(1,n); zetaT  = x(2,n); zetaG  = x(3,n);
st = x(4,n); sg = x(5,n); sy = x(6,n);

D  = 1/(1-zetaT*thetaY)*[1 thetaY*zetaG+thetaG thetaY;
                         0 (1-zetaT*thetaY) 0 ;
                         zetaT zetaT*thetaG+zetaG 1]*diag([st sg sy]);
DD = D*D';

sel=tril(ones(3,3));
dif=DD(sel==1)-Sigma(sel==1);
z(:,n) = dif(:);
n=n+1;      
end