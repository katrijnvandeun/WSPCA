%fixed model & algorithm parameters
NRSTARTS=10;       %number of starts in the multistart procedure
MAXITER=200;     %maximum number of iterations
SCALING='off';%'on' or 'off'
OFFSET='off';
LASSOTYPE='ordinary';
CONVERGENCE=1e-4;
HISTORY=0;
LASSOt=0;%sqrt(n1+n2)*quantile(abs(scores(:)),quants(lt));
orth=1;%1 or 0 for orthogonal or oblique

RESULTunif=[];
%load RESULTrandom
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
    
    %********
    %1. Non-sparse, non-weighted PCA
    W_unw=ones(size(DATA));
    W_unw(W==0)=0;
    LASSOP=0;
    [Ti,Pi,ci,si,Bi,Lossi]=WSPCA(DATAmiss,W_unw,R,LASSOP,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'rational',[]);
    LOSS=Lossi
    resultPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si);
    for nrstart=1:NRSTARTS;
        [Ti,Pi,ci,si,Bi,Lossi]=WSPCA(DATAmiss,W_unw,R,LASSOP,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'random',[]);
        if Lossi<LOSS
            resultPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si,'Loss',LOSS);
        end;
    end;
    DATAhat=resultPCA.T*resultPCA.P';
    devtrue=Wmiss.*(TRUE-DATAhat);
    presstrue_pca=sum(sum(devtrue.^2))/sum(sum(((Wmiss).*TRUE).^2))
    [~,~,tucker_pca,~]=tuckercongruence_pr2(TTRUE,resultPCA.T)
    [~,~,tucker_pcaP,~]=tuckercongruence_pr2(PTRUE,resultPCA.P)
    
    %********
    %2. Non-sparse, weighted PCA
    LASSOP=0;
    [Ti,Pi,ci,si,Bi,Lossi]=WSPCA(DATAmiss,W.^(0.5),R,LASSOP,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'rational',[]);
    LOSS=Lossi;
    resultwPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si);
    for nrstart=1:NRSTARTS;
        [Ti,Pi,ci,si,Bi,Lossi]=WSPCA(DATAmiss,W.^(0.5),R,LASSOP,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'random',[]);
        if Lossi<LOSS
            LOSS=Lossi;
            resultwPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si,'Loss',LOSS);
        end;
    end;
    DATAhat=resultwPCA.T*resultwPCA.P';
    devtrue=Wmiss.*(TRUE-DATAhat);
    presstrue_wpca=sum(sum(devtrue.^2))/sum(sum(((Wmiss).*TRUE).^2))
    [~,~,tucker_wpca,~]=tuckercongruence_pr2(TTRUE,resultwPCA.T)
    [~,~,tucker_wpcaP,~]=tuckercongruence_pr2(PTRUE,resultwPCA.P)
    
    %********
    %3. Sparse, weighted PCA
    %**Cardinality constraint imposed
    if propsparse==0.5
        CARD=sum(PTRUE==0,1);
    else
        J = size(PTRUE,1);
        CARD = 0.5*[J J];
    end;
    %Note: orthogonal T and sparse P => estimation of each p_jr is independent
    %(univariate soft thresholding)
    [Ti,Pi,ci,si,Bi,Lossi]=WSPCA(DATAmiss,W.^(0.5),R,0,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'rational',CARD);
    LOSS=Lossi;
    resultwSPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si);
    for nrstart=1:NRSTARTS;
        [Ti,Pi,ci,si,Bi,Lossi]=WSPCA(DATAmiss,W.^(0.5),R,0,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'random',CARD);
        if Lossi<LOSS
            LOSS=Lossi;
            resultwSPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si,'Loss',LOSS);
        end;
    end;
    
    DATAhat=resultwSPCA.T*resultwSPCA.P';
    devtrue=Wmiss.*(TRUE-DATAhat);
    presstrue_wspca=sum(sum(devtrue.^2))/sum(sum(((Wmiss).*TRUE).^2))
    [~,~,tucker_wspca,~]=tuckercongruence_pr2(TTRUE,resultwSPCA.T)
    [~,~,tucker_wspcaP,~]=tuckercongruence_pr2(PTRUE,resultwSPCA.P)
    
    %Store results per 10 datasets
    [I J]=size(DATA);
    RESULTunif=[RESULTunif;noiseadd noisemult propmiss propsparse size(DATA,1) size(DATA,2) R presstrue_pca ...
        presstrue_wpca presstrue_wspca ...
        tucker_pca tucker_wpca tucker_wspca tucker_pcaP tucker_wpcaP tucker_wspcaP ...
        (sum(sum(W==0)))/(I*J) (sum(sum(W>0.2)))/(I*J) (sum(sum(W>0.5)))/(I*J) (sum(sum(W>0.9)))/(I*J)];
    if rem(teller,10)==0
        save RESULTnomisssqrtPunifb RESULTunif
    end;
    teller=teller+1
end;
dlmwrite(['../R/RESULTnomiss_sqrtW_Puniofb.txt'],RESULTunif,'delimiter','\t')