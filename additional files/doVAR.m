function VAR = doVAR(VAR)
 VARSLAGS = lagmatrix(VAR.vars,1:VAR.p);
 VARSLAGS = VARSLAGS(VAR.p+1:end,:);
 VARS     = VAR.vars(VAR.p+1:end,:);

 [VAR.T,VAR.n]    = size(VARS);
 
% Run VAR
%%%%%%%%%
 VAR.bet=[VARSLAGS VAR.DET(VAR.p+1:end,:)]\VARS; 
 res = VARS-[VARSLAGS VAR.DET(VAR.p+1:end,:)]*VAR.bet;
 VAR.Sigma = (res'*res)/(VAR.T-VAR.n*VAR.p-1); 
 
%Identification
%%%%%%%%%%%%%%%%
 
         if VAR.d0 ==0
             d0 = [0; 0;0;diag(VAR.Sigma)];
         else
             d0 = VAR.d0;
         end
         [sol,VAR.rc]=csolve(@Dsolve,d0(:),[],1.e-8,1000,VAR.BP.thetaY...
            ,VAR.BP.gammaT,VAR.BP.gammaY,VAR.Sigma); 
        
         VAR.thetaG = sol(1); VAR.zetaT = sol(2); VAR.zetaG = sol(3);
         VAR.sigmaT = sol(4); VAR.sigmaG = sol(5); VAR.sigmaY = sol(6); VAR.thetaY=VAR.BP.thetaY;
         [z,VAR.A,VAR.B,VAR.D] = Dsolve(sol,VAR.BP.thetaY,VAR.BP.gammaT,VAR.BP.gammaY,VAR.Sigma);
         VAR.d0 = sol;
 
 
 
% Impulse Responses (Tax Shock)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 irs(VAR.p+1,:) = -100*VAR.D(:,1)/VAR.D(1,1)*VAR.tshocksize;

 for j=2:VAR.irhor
 lvars = (irs(VAR.p+j-1:-1:j,:))';
 irs(VAR.p+j,:) = lvars(:)'*VAR.bet(1:VAR.p*VAR.n,:);     
 end
 VAR.irs = irs(VAR.p+1:end,:);

 irsg(VAR.p+1,:) = 100*VAR.D(:,2)/VAR.D(2,2)*VAR.gshocksize;
 
 for j=2:VAR.irhor
 lvarsg = (irsg(VAR.p+j-1:-1:j,:))';
 irsg(VAR.p+j,:) = lvarsg(:)'*VAR.bet(1:VAR.p*VAR.n,:);     
 end
 VAR.irsg = irsg(VAR.p+1:end,:);
 
 VAR.e             = (VAR.D(:,:,1)\res')';  

function [z,A,B,D] = Dsolve(x,thetaY,gammaT,gammaY,Sigma)
[rows,cols]=size(x);

n=1;
while n<=cols; 
thetaG = x(1,n); zetaT  = x(2,n); zetaG  = x(3,n);
st = x(4,n); sg = x(5,n); sy = x(6,n);

A = [1 0 -thetaY; 0 1 -gammaY; -zetaT -zetaG 1];
B = [st thetaG*sg 0; gammaT*st sg 0; 0 0 sy];

D = A\B;

DD = D*D';

sel=tril(ones(3,3));
dif=DD(sel==1)-Sigma(sel==1);
z(:,n) = dif(:);
n=n+1;      
end




