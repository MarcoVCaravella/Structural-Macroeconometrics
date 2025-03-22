function MA = doMA(MA)
 
 [T,n]    = size(MA.vars);

 VARS    = diff(MA.vars);
 X    = [MA.TSHOCKS(2:end,:) MA.DET(2:end,:)];

 bet=X\VARS ;
 [T,k]=size(X);
 
 MA.irs(1,:) = bet(1,:);
 for i = 2:MA.q+1
 MA.irs(i,:) = MA.irs(i-1,:)+bet(i,:);
 end
 MA.irs=MA.irs*MA.shocksize;
 
 res = VARS-[MA.TSHOCKS(2:end,:) MA.DET(2:end,:)]*bet;
 Sigma = (res'*res)/(T-k);
 
 
 for jj = 1:3
 sigy = Sigma(jj,jj); 
 varbet = inv(X'*X)*sigy;
 SE(1) = sqrt(varbet(1,1)); 
 for j = 2:MA.q+1
 SE(j) = sqrt(SE(j-1)^2+varbet(j,j)+2*sum(varbet(j,1:j-1)));
 end
 MA.SE(:,jj) = SE';
 end
 
MA.irsH = MA.irs+2*MA.SE*100;
MA.irsL = MA.irs-2*MA.SE*100;