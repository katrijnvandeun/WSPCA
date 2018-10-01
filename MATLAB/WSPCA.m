function [T,P,c,s,B,loss]=WSPCAfast(DATA,W,R,LP,offset,scaling,type,MAXITER,CONVERGENCE,HISTORY,LT,Tinit,Pinit,orth,INIT,CARD)
%WSPCAfast performs weighted sparse PCA with estimation of offset
%and scale.
%Objective functions:   ||W.*((X-1c')S-TP'||^2+LP|P|_1+LT|T|_1
%                   OR
%                       ||W.*(X-TP'||^2 s.t. T'T=I and card(P)=p*
%INPUT: DATA:  IxJ matrix of data
%       W:     IxJ matrix of known weights: nonnegative & smaller than one
%       LP:    nonnegative value to tune the L1 penalty for the loadings
%       LT:    nonnegative value to tune the L1 penalty for the scores
%OUTPUT: T: IxR matrix of component scores
%        P: JxR matrix of loadings
%        c: Jx1 vector of offsets
%        s: Jx1 vector of scaling parameters
%Author: K. Van Deun, Dept. Methodology & Statistics, Tilburg University
%and Department of Quantitative Psychology, KU Leuven
%First version: APRIL 2013, REVISED SEPTEMBER 2018

[I,J]=size(DATA);

%Initialize parameters subject to active constraints
switch INIT
    case 'rational'
        DATAinit=DATA;
    case 'random'
        DATAinit=randn(I,J);
    case 'semi-rational'
        DATAinit=DATA+randn(I,J);
end;
switch offset
    case 'off'
        c=zeros(1,J);
    case 'on'
        if (length(c)==0)
            c=mean(DATAinit);
        end;
end;
c=c';
Pzero=find(Pinit==0);%find zero constraints
P=Pinit;
T=Tinit;
if (length(Tinit)==0 | length(Pinit)==0)
    [U S V]=svds(DATAinit,R);
    if length(Tinit)==0
        T=U;
    end;
    if length(Pinit)==0
        P=V*S;
        P(Pzero)=0;
    end;
end;

switch scaling
    case 'on'
        s=sqrt(var(DATAinit));
        s=sqrt(J)*s/norm(s);
    case 'off'
        s=ones(J,1);
end;

%Calculation of elements remaining constant in the iterative procedure
wmax=max(max(W));
wWsq=W.*W*(wmax^(-2));
v=ones(I,1);

%Penalty weights: adaptive lasso, randomized lasso, ordinary lasso
if isempty(CARD)
    switch type
        case 'ordinary'
            B=ones(J,R);
            B_T=ones(I,R);
        case 'adaptive'
            for i=1:2
                P=Update_loadings(W,DATA,c,s,T,P,0,[],[]);
            end;
            B=abs(P).^(-1);
            for i=1:2
                T=Update_scores(wWsq,DATA,c,s,T,P,0,[],orth);
            end;
            B_T=abs(T).^(-1);
        case 'randomized'
            weakness=0.2;%value "recommended" by Meinshausen & Buhlmann, 2010
            B=round(rand(J,R));
            B(B==0)=weakness;
            B_T=round(rand(I,R));
            B_T(B_T==0)=weakness;
    end;
else
    if orth~=1
        fprintf('WARNING: cardinality constraint only works with orthogonal components \n')
        orth=1;
    end;
    LP=0;%no penalty in calculation loss
    B=ones(J,R);%required for calculation loss
    B_T=ones(I,R);
    P=Update_loadings(W,DATA,c,s,T,P,0,B,CARD);
end;

%calculation of initial loss
lossold=WSPCALOSS(DATA,W,LP,T,P,B,c,s,LT,B_T);
conv=0;%convergence
iter=1;%iteration counter
P(abs(P)<eps)=0;%????not needed

while conv~=1
    LOSS=lossold;
    
    %1. Conditional estimation of c
    if strncmp(offset,'on',2)
        c=Update_offset(wWsq,DATA,c,s,T,P);
        %else
        %    c=zeros(size(c));%!!!
    end;
    loss=WSPCALOSS(DATA,W,LP,T,P,B,c,s,LT,B_T);
    if HISTORY==1
        fprintf('offset: \t %4.0f Loss: %8.4f Diff: %20.12f \n',iter,loss,lossold-loss);
    end;
    if lossold-loss<0
        lossold-loss
    end;
    lossold=loss;
    
    %2. Conditional estimation of s
    if iter<3 & strncmp(scaling,'on',2)
        [s P]=Update_scale(wWsq,DATA,c,s,T,P);
    end;
    loss=WSPCALOSS(DATA,W,LP,T,P,B,c,s,LT,B_T);
    if HISTORY==1
        fprintf('scale: \t \t %4.0f Loss: %8.4f Diff: %20.12f \n',iter,loss,lossold-loss);
    end;
    lossold=loss;
    
    %3. Conditional estimation of T
    T=Update_scores(wWsq,DATA,c,s,T,P,LT,B_T,orth);
    
    loss=WSPCALOSS(DATA,W,LP,T,P,B,c,s,LT,B_T);
    if HISTORY==1
        fprintf('scores: \t %4.0f Loss: %8.4f Diff: %20.12f \n',iter,loss,lossold-loss);
    end;
    lossold=loss;
    
    %4. Conditional estimation of P
    P=Update_loadings(W,DATA,c,s,T,P,LP,B,CARD);
    P(Pzero)=0;
    
    loss=WSPCALOSS(DATA,W,LP,T,P,B,c,s,LT,B_T);
    if HISTORY==1
        fprintf('loadings: \t %4.0f Loss: %8.4f Diff: %20.12f \n',iter,loss,lossold-loss);
    end;
    lossold=loss;
    
    %5. Check stop criteria
    %5.1. max iterations?
    if iter==MAXITER
        conv=1;
    end;
    iter=iter+1;
    %5.2. convergence?
    if LOSS-lossold<CONVERGENCE
        conv=1;
    end;
    
end;