function [T,P,c,s,B,loss]=WSPCAnoupdateP(DATA,W,R,LASSO,offset,scaling,lasso,MAXITER,CONVERGENCE,HISTORY,LASSOT,T,P,c,orth,INIT)
%WSPCAfast performs weighted sparse PCA with estimation of offset
%and scale using warm restarts and calculations on non-zero loading vectors
%Objective function: ||W.*((X-1c')S-TP'||^2+lasso|P|_1
%INPUT: DATA:  IxJ matrix of data
%       W:     IxJ matrix of known weights: nonnegative & smaller than one
%       lasso: nonnegative value to tune the adaptive L1 penalty
%OUTPUT: T: IxR matrix of component scores
%        P: JxR matrix of loadings
%        c: Jx1 vector of offsets
%        s: Jx1 vector of scaling parameters
%Author: K. Van Deun, Dept. Psychology, KU Leuven
%First version: APRIL 2013

[I,J]=size(DATA);
eps=1e-9;%numerical precision threshold
minspread=1/(10*J);%minimal spread threshold
weakness=0.2;%value "recommended" by Meinshausen & Buhlmann

%Initialize parameters subject to active constraints
switch INIT
    case 'rational'
        DATAinit=DATA;
    case 'random'
        DATAinit=randn(I,J);
    case 'semi-rational'
        DATAinit=DATA+randn(I,J);
end;
%Initialize parameters subject to active constraints
switch offset
    case 'off'
        c=zeros(1,J);
    case 'on'
        if (length(c)==0)
            c=mean(DATAinit);
        end;
end;
c=c';
if (length(T)==0 | length(Pinit)==0)
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
switch lasso
    case 'ordinary'
        B=ones(J,R);
        B_T=ones(I,R);
    case 'adaptive'
%        [U S V]=svds(DATA,R);
        for i=1:2
            P=Update_loadings(W,DATA,c,s,T,P,0,[]);%estim
        end;
        B=abs(P).^(-1);
        for i=1:2
            T=Update_scores(wWsq,DATA,c,s,T,P,0,[],orth);%estim
        end;
       B_T=abs(T).^(-1);
    case 'randomized'
        B=round(rand(J,R));
        B(B==0)=weakness;
        B_T=round(rand(I,R));
        B_T(B_T==0)=weakness;
end;
%calculation of initial loss
lossold=WSPCALOSS(DATA,W,LASSO,T,P,B,c,s,LASSOT,B_T);
conv=0;%convergence
iter=1;%iteration counter
%Pinit(abs(Pinit)<eps)=0;
%P=Pinit;

while conv~=1
    LOSS=lossold;
    
        %1. Conditional estimation of c
        if strncmp(offset,'on',2)
            c=Update_offset(wWsq,DATA,c,s,T,P);
        %else
        %    c=zeros(size(c));%!!!
        end;
        loss=WSPCALOSS(DATA,W,LASSO,T,P,B,c,s,LASSOT,B_T);
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
        loss=WSPCALOSS(DATA,W,LASSO,T,P,B,c,s,LASSOT,B_T);
        if HISTORY==1
            fprintf('scale: \t \t %4.0f Loss: %8.4f Diff: %20.12f \n',iter,loss,lossold-loss);
        end;
        lossold=loss;
    
    %3. Conditional estimation of T
    T=Update_scores(wWsq,DATA,c,s,T,P,LASSOT,B_T,orth);
    %fix length T
    %    normT=sqrt(min(sum(T.^2)));
    %     T=T/normT;
    %     P=P*normT;
    
    loss=WSPCALOSS(DATA,W,LASSO,T,P,B,c,s,LASSOT,B_T);
    if HISTORY==1
        fprintf('scores: \t %4.0f Loss: %8.4f Diff: %20.12f \n',iter,loss,lossold-loss);
    end;
    lossold=loss;
    
    %4. Check stop criteria
    %4.1. max iterations?
    if iter==MAXITER
        conv=1;
    end;
    iter=iter+1;
    %4.2. convergence?
    if LOSS-lossold<CONVERGENCE
        conv=1;
    end;
    
end;