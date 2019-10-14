function Loss=pstrLoss(B,T,Target,W)
%pstrLoss calculates the objective function associated to the pstr m-file
%
%Author: Katrijn Van Deun
DEV=T*B-Target;
wDEV=W.*DEV;
Loss=sum(sum(wDEV.^2));