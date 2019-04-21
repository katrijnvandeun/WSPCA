function [loadings indexnonzero indexzero]=Update_loadingsCardinalityV2(W,DATA,c,s,T,P,CARD)
%UPDATE_LOADINGSCARDINALITY Calculates the update of the loadings under a
%cardinality constraint
%   Objective: minimize L(T)=||W°((X-1c')S-TP^T)|| s.t. P has given nr of
%   zeros
%   This constrained WLS problem is reformulated as a constrained OLS
%   problem using MM. The constrained PCA problem is solved as explained in
%   Shen & Huang (2008).
%
%Author: Katrijn Van Deun, SEPTEMBER 2018

%precision
%eps=1e-6;

[J,R]=size(P);
[I,R]=size(T);
if sum(CARD)==0  %Weighted unconstrained PCA
    wmax=max(max(W));
    XHAT=T*P';
    wWsq=W.*W*(wmax^(-2));
    RES=residual(DATA,c,s,T,P);
    YSTAR=ystar(XHAT,RES,wWsq);
    P=YSTAR'*T;
    loadings=P;
else
    TDATA=(DATA-(ones(I,1)*c')).*(ones(I,1)*s');
    XHAT=T*P';
    Wsq=W.^2;
    RES=TDATA-XHAT;
    YSTAR=XHAT+(Wsq.*RES);
    TYSTAR=YSTAR'*T;
    [~, i]=sort(abs(TYSTAR(:)),'descend');
    P=zeros(J,R);
    P(i(1:CARD))=TYSTAR(i(1:CARD));
    loadings=P;
end