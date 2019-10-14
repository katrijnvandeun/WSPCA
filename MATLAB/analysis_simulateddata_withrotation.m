%fixed model & algorithm parameters
NRSTARTS=10;       %number of starts in the multistart procedure
MAXITER=200;     %maximum number of iterations
SCALING='off';%'on' or 'off'
OFFSET='off';
LASSOTYPE='ordinary';
CONVERGENCE=1e-4;
HISTORY=0;
LASSOt=0;%no penalty on component scores here;
orth=1;%1 or 0 for orthogonal or oblique

RESULTunif=[];
%load RESULTreviewer2 %rotation of (W)PCA to true P suggested by reviewer
for teller=1:720;
    fname=sprintf('../DATA/DATAnomissc%d.dat',teller)
    DATA=importdata(fname);
    k=find(abs(DATA)==Inf);
    DATA(k)=0;
    fname=sprintf('../DATA/WEIGHTSnomissc%d.dat',teller);
    W=importdata(fname);
    W(k)=0;
    fname=sprintf('../DATA/TRUEnomissc%d.dat',teller)
    TRUE=importdata(fname);
    fname=sprintf('../DATA/TRUEPnomissc%d.dat',teller)
    PTRUE=importdata(fname);
    fname=sprintf('../DATA/TRUETnomissc%d.dat',teller)
    TTRUE=importdata(fname);
    fname=sprintf('../DATA/STRUCTnomissc%d.mat',teller)
    load(fname)
    R=setting.R;
    propsparse=setting.PrSparse;
    propmiss=setting.PrMiss;
    noiseadd=setting.AddNoise;
    noisemult=setting.MultNoise;
    
    W=W/(max(max(W)));
    
    
    %***********************************
    %            ANALYSES
    %***********************************
    DATAmiss=DATA;
    DATAmiss(W==0)=0;
    Wmiss=W;
    Wmiss(W~=0)=1;
    
    sum(sum(W==0))
    %********
    %1. Non-sparse, non-weighted PCA WITH ROTATION
    W_unw=ones(size(DATA));
    W_unw(W==0)=0;
    LASSOP=0;
    [Ti,Pi,ci,si,Bi,Lossi]=WSPCAfast(DATAmiss,W_unw,R,LASSOP,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'rational');
    LOSS=Lossi
    resultPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si);
    for nrstart=1:NRSTARTS;
        [Ti,Pi,ci,si,Bi,Lossi]=WSPCAfast(DATAmiss,W_unw,R,LASSOP,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'random');
        if Lossi<LOSS
            resultPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si,'Loss',LOSS);
        end;
    end;
    %rotation to true P
    [B Loss]=pstr(resultPCA.P,PTRUE,ones(size(PTRUE)),200,1e-4);
    Protated = resultPCA.P*B;
    Trotated = resultPCA.T*B;
    DATAhat=Trotated*Protated';
    devtrue=Wmiss.*(TRUE-DATAhat);
    presstrue_pca=sum(sum(devtrue.^2))/sum(sum(((Wmiss).*TRUE).^2))
    [~,~,tucker_pca_rotated,~]=tuckercongruence_pr2(TTRUE,Trotated)
    [~,~,tucker_pcaP_rotated,~]=tuckercongruence_pr2(PTRUE,Protated)
    
    %********
    %2. Non-sparse, weighted PCA  WITH ROTATION
    LASSOP=0;
    [Ti,Pi,ci,si,Bi,Lossi]=WSPCAfast(DATAmiss,W.^(0.5),R,LASSOP,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'rational');
    LOSS=Lossi;
    resultwPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si);
    for nrstart=1:NRSTARTS;
        [Ti,Pi,ci,si,Bi,Lossi]=WSPCAfast(DATAmiss,W.^(0.5),R,LASSOP,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'random');
        if Lossi<LOSS
            LOSS=Lossi;
            resultwPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si,'Loss',LOSS);
        end;
    end;
    [B Loss]=pstr(resultwPCA.P,PTRUE,ones(size(PTRUE)),200,1e-4);
    wProtated = resultwPCA.P*B;
    wTrotated = resultwPCA.T*B;
    DATAhat=wTrotated*wProtated';
    devtrue=Wmiss.*(TRUE-DATAhat);
    presstrue_wpca=sum(sum(devtrue.^2))/sum(sum(((Wmiss).*TRUE).^2))
    [~,~,tucker_wpca_rotated,~]=tuckercongruence_pr2(TTRUE,wTrotated)
    [~,~,tucker_wpcaP_rotated,~]=tuckercongruence_pr2(PTRUE,wProtated)
    

    
    %Store results per 10 datasets
    [I J]=size(DATA);
    RESULTunif=[RESULTunif;noiseadd noisemult propmiss propsparse size(DATA,1) size(DATA,2) R ...
        tucker_pca_rotated tucker_wpca_rotated tucker_pcaP_rotated tucker_wpcaP_rotated ...
        (sum(sum(W==0)))/(I*J) (sum(sum(W>0.2)))/(I*J) (sum(sum(W>0.5)))/(I*J) (sum(sum(W>0.9)))/(I*J)];
    if rem(teller,10)==0
        save RESULTreviewer2 RESULTunif
    end;
    teller=teller+1
end;
dlmwrite(['../R/reviewer2.txt'],RESULTunif,'delimiter','\t')