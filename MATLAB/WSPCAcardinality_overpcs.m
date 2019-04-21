function [T,P,c,s,B,loss]=WSPCAcardinality_overpcs(DATA,W,R,LP,offset,scaling,type,MAXITER,CONVERGENCE,HISTORY,LT,T,P,orth,INIT,CARD)
%WSPCAfast performs weighted sparse PCA with estimation of offset
%and scale.
%Objective function: ||W.*((X-1c')S-TP'||^2+LP|P|_1 s.t. Card(P)=CARD
%INPUT: DATA:  IxJ matrix of data
%       W:     IxJ matrix of known weights: nonnegative & smaller than one
%       LP:    nonnegative value to tune the L1 penalty for the loadings
%       LT:    nonnegative value to tune the L1 penalty for the scores
%       CARD:  number of non-zero coefficients OVER all components
%OUTPUT: T: IxR matrix of component scores
%        P: JxR matrix of loadings
%        c: Jx1 vector of offsets
%        s: Jx1 vector of scaling parameters
%Author: K. Van Deun, Dept. Methodology & Statistics, Tilburg University
%First version: SEPTEMBER 2018

[I,J]=size(DATA);
weakness=0.2;%value "recommended" by Meinshausen & Buhlmann, 2010
LP=0;
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
if (length(T)==0 | length(P)==0)
    [U S V]=svds(DATAinit,R);
    if length(T)==0
        T=U;
    end;
    if length(P)==0
        P=V*S;
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
%ordinary lasso
switch type
    case 'ordinary'
        B=ones(J,R);
        B_T=ones(I,R);
    case 'adaptive'
        for i=1:2
            P=Update_loadings(W,DATA,c,s,T,P,0,[]);
        end;
        B=abs(P).^(-1);
        for i=1:2
            T=Update_scores(wWsq,DATA,c,s,T,P,0,[],orth);
        end;
        B_T=abs(T).^(-1);
    case 'randomized'
        B=round(rand(J,R));
        B(B==0)=weakness;
        B_T=round(rand(I,R));
        B_T(B_T==0)=weakness;
end;

%create P subject to active constraints
P=Update_loadingsCardinalityV2(W,DATA,c,s,T,P,CARD);

%calculation of initial loss
lossold=WSPCALOSS(DATA,W,LP,T,P,B,c,s,LT,B_T);
conv=0;%convergence
iter=1;%iteration counter
P(abs(P)<eps)=0;

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
    P=Update_loadingsCardinalityV2(W,DATA,c,s,T,P,CARD);
    
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