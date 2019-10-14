function [B Loss]=pstr(T,Target,W,maxiter,convergence)
%PSTR finds a (orthogonal) rotation to a partially specified target
%The following objective function is used:
%       min_B ||Wo(TB-Target)||²
%with W a binary matrix of weights, T the matrix to rotate and Target the
%target to reach.
%The procedure to minimize the weighted least squares problem relies on
%iterative majorization.
%INPUT   T: matrix to rotate
%        Target: target to reach
%        W: weight matrix (0/1 for resp. unspecified and specified elements of the target)
%        maxiter: maximum number of iterations
%        convergence: minimal loss between current and previous iteration
%OUTPUT  B: an orthogonal rotation matrix
%
%Author: Katrijn Van Deun

[n m]=size(T);

L=[];
BMAT=[];
REFL=reflexmat(m);
for i=1:size(REFL,1)
    k=find(REFL(i,:)==-1);
    Binit=eye(m);
    Binit(:,k)=-1*Binit(:,k);
    
    B1=T'*T;
    alpha=max(eig(B1));
    iter=1;
    stop=0;
    if i==1
        Bcurrent=Binit;%Binit1;
        Lossc=pstrLoss(Binit,T,Target,W);
    else
        Bcurrent=Binit;%Binit2;
        Lossc=pstrLoss(Binit,T,Target,W);
    end;

    while stop==0
        Tw=W.*Target+T*Bcurrent-W.*(T*Bcurrent);
        A=-2*Tw'*T;
        F=A+2*Bcurrent'*B1'-2*alpha*Bcurrent';
        [U S V]=svd(-F);
        B=V*U';
        if iter==maxiter
            stop=1;
        end;
        Loss=pstrLoss(B,T,Target,W);
        Diff=Lossc-Loss;
        if abs(Diff)<convergence
            stop=1;
        end;
        %fprintf('Iteration Nr: %3.0f \t Loss: %5.4f \t  Diff: %5.4f \n',iter,Loss,Diff);
        iter=iter+1;
        Lossc=Loss;
        Bcurrent=B;
    end;
    %fprintf('Iteration Nr: %3.0f \t Loss: %5.4f \t  Diff: %5.4f \n',iter,Loss,Diff);
    L=[L Lossc];
    BMAT=[BMAT Bcurrent];
end;
k=find(L==min(L));
Loss=L(k(1));
B=BMAT(:,m*(k(1)-1)+1:m*k(1));
%fprintf('Iteration Nr: %3.0f \t Loss: %5.4f \t  Diff: %5.4f \n',iter,Loss,Diff);