function scores =Update_scores(W,DATA,c,s,T,P,LASSOT,B,orth)
%UPDATE_scores Calculates the update of the scores in WSPCA
%   Objective: minimize L(T)=||W°((X-1c')S-TP^T)||
%   The WLS problem is solved by solving the problem for a majorizing
%   function; the problem than becomes an oblique Procrustes rotation
%   problem which can be solved by restricted least squares
%
%Author: Katrijn Van Deun, APRIL 2013, Revised October 2014
eps=1e-9;
[I,R]=size(T);
RES=residual(DATA,c,s,T,P);
XHAT=T*P';
YSTAR=ystar(XHAT,RES,W);
Told=T;
if isempty(B)
    B=ones(size(T));
end;
if orth==1 %orthogonal (non-sparse) component scores
    [U S V]=svds(P'*YSTAR',size(P,2));
    T=V*U';
else
    for r=1:R
        %restricted least squares
        Tminr=T;
        Tminr(:,r)=[];
        Pminr=P;
        Pminr(:,r)=[];
        XHATr=Tminr*Pminr';
        E=YSTAR-XHATr;
        t=E*P(:,r);
        %soft thresh-hold
        T(:,r)=sign(t).*(abs(t)-((LASSOT/2)*B(:,r)));
        T(abs(t)<((LASSOT/2)*B(:,r)),r)=0;
        if sum(abs(T(:,r)))~=0%avoid division by zero
            T(:,r)=T(:,r)/norm(T(:,r));
        else 
            T(:,r)=zeros(I,1);
        end;
     end;
end;
scores=T;
