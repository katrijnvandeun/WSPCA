function [PrMAT]=subsamplingWSPCAcard(DATA,W,R,LASSO,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,initP,initT,orth,MULTISTART,nsamples,splitratio,q)
%SUBSAMPLING Stability selection for weighted sparse 
%Calculates the relative frequency of obtaining a non-zero component loading for WSPCA.
%
%INPUT: Input parameters for sparse covariates regression
%       L_int: interval of lasso penalty tuning values (descending)
%       nsamples: number of resamples
%       splitratio: ratio of sample size resampled data to sample size
%       original data
%
%OUTPUT: PrMAT: matrix with relative frequency of non-zero component weight
%
%K. Van Deun, OCT2015; adapted to WSPCA on 28th SEPTEMBER 2018
[I J]=size(DATA);
sample_size=round(splitratio*I);

PrMAT=zeros(J,R);
%reference for determ. comp.scores
[Ti,Pi,ci,si,Bi,Lossi]=WSPCA(DATA,W.^(0.5),R,LASSO,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,MULTISTART,[]);
refmatrix=Ti;
for n=1:nsamples
    v=randperm(I,round(splitratio*I));
    X_s=DATA(v,:);
    W_s=W(v,:);
    [Ti_s,Pi,ci,si,Bi,Lossi]=WSPCA(X_s,W_s.^(0.5),R,LASSO,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,MULTISTART,q);
    if R>1
        newmatrix=Ti_s;
        [perm reflex tucker tuckervector]=tuckercongruence_pr2(refmatrix(v,:),newmatrix);
        Pi=Pi(:,perm);
        PrMAT=PrMAT+(Pi~=0);
    else
        PrMAT=PrMAT+(Pi~=0);
    end;
end;
PrMAT=PrMAT/nsamples;
fname=['PrMAT_N250card'];
save([fname],'PrMAT')