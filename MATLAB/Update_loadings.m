function [loadings]=Update_loadings(W,DATA,c,s,T,P,LASSO,B,CARD)
%UPDATE_LOADINGS Calculates the update of the offsets in WSPCA
%   Objective: minimize L(T)=||W.*((X-1c')S-TP^T)||+\lambda|B.*P|_1
%   The WLS problem is solved by coordinate descent (Friedman et al., 2007).
%   B is a matrix of weights: The OLS estimates of P (adaptive lasso, Zou, 2006) or
%   random (randomized lasso; Buhlmann & Meinshausen, 2010)
%
%Author: Katrijn Van Deun, APRIL 2013

%precision
%eps=1e-6;

[J,R]=size(P);
[I,R]=size(T);
Wsq=W.^2;
if LASSO<=0 & isempty(CARD) %Weighted unconstrained PCA
    wmax=max(max(W));
    XHAT=T*P';
    wWsq=Wsq.*(wmax^(-2));
    RES=residual(DATA,c,s,T,P);
    YSTAR=ystar(XHAT,RES,wWsq);
    P=YSTAR'*T;
    loadings=P;
else
    TDATA=(DATA-(ones(I,1)*c')).*(ones(I,1)*s');
    if isempty(CARD) %penalized approach
        for r=1:R
            Tminr=T;
            Tminr(:,r)=[];
            Pminr=P;
            Pminr(:,r)=[];
            tsq=(T(:,r)).^2;
            Twsq=tsq'*Wsq;
            RES=TDATA-Tminr*Pminr';
            TX=(T(:,r)'*(Wsq.*RES))'./(Twsq');
            TX(Twsq==0)=0;
            Bw=B(:,r)./(Twsq');
            Bw(Twsq==0)=0;
            P(:,r)=sign(TX).*(abs(TX)-((LASSO/2)*Bw));
            P(abs(TX)<((LASSO/2)*(Bw)),r)=0;
        end;
    else
        XHAT=T*P';
        RES=TDATA-XHAT;
        YSTAR=XHAT+(Wsq.*RES);
        TYSTAR=YSTAR'*T;
        P=zeros(J,R);
        if length(CARD)>1 %constraint per component
            for r=1:R
                %TWO OPTIONS: strict cardinality constraint or internal
                %updating of LASSO (Shen & Huang, 2007; Lee et al., 2010)
                %The strict cardinality constraint relies on the
                %majorization scheme and gave bad results in the simulation
                %OPTION 1
                %[~, i]=sort(abs(TYSTAR(:,r)),'descend');
                %P(i(1:CARD(r)),r)=TYSTAR(i(1:CARD(r)),r);
                %OPTION 2
                Tminr=T;
                Tminr(:,r)=[];
                Pminr=P;
                Pminr(:,r)=[];
                tsq=(T(:,r)).^2;
                Twsq=tsq'*Wsq;
                RES=TDATA-Tminr*Pminr';
                TX=(T(:,r)'*(Wsq.*RES))'./(Twsq');
                TX(Twsq==0)=0;
                Bw=B(:,r)./(Twsq');
                Bw(Twsq==0)=0;
                TXbw=TX./Bw;
                [~, i]=sort(abs(TXbw),'descend');
                LASSO=2*abs(TXbw(i(CARD(r))));
                P(:,r)=sign(TX).*(abs(TX)-((LASSO/2)*Bw));
                P(abs(TX)<((LASSO/2)*(Bw)),r)=0;
            end;
        else %constraint over components
            [~, i]=sort(abs(TYSTAR(:)),'descend');
            P(i(1:CARD))=TYSTAR(i(1:CARD));
        end;
        loadings=P;
    end;
end;