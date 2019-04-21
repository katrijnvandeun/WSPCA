%last update: 19th February 2019

%***********************************
%    LOADING OF DATA
%***********************************
load ../DATA/DATA2008
DATA=DATA_diffD3';
load ../DATA/DATA2007
DATA2007=DATA_diffD3';
[I J]=size(DATA);
[I2007 J2007]=size(DATA2007);
%Load Criterion variable
load ../DATA/TIVtiter2008
y=log2(Y(:,4));%4th column is max of three titer types
load ../DATA/TIVtiter2007
y2007=log2(Y(:,4));%4th column is max of three titer types

%load Rocke Lorenzato/Durbin variance parameters and D0 and D3 data
load ../R/GSE29617/rawdata/GSE29617_RLparsD0.txt;
load ../R/GSE29617/rawdata/GSE29617_RLparsD3.txt;
load ../DATA/DATA_D0
load ../DATA/DATA_D3
%Calculate WEIGHTS according to Rocke Lorenzato model: this is difference
%of expression scores (D3-D0) so we use 4.20 of the Rocke & Durbin (2001) paper
%Because of background correction, we consider offsets to be equal to zero
epsilonD0=GSE29617_RLparsD0(2);
etaD0=GSE29617_RLparsD0(3);
epsilonD3=GSE29617_RLparsD3(2);
etaD3=GSE29617_RLparsD3(3);
varlnexpr=etaD0^2+etaD3^2+(epsilonD0^2*(DATA_D0).^(-1))+(epsilonD3^2*(DATA_D3).^(-1));
WEIGHTS=(1./varlnexpr);
WEIGHTS=WEIGHTS';

%Select shared probesets between 2008 and 2007 data
DATA=DATA(:,1:J2007);
WEIGHTS=WEIGHTS(:,1:J2007);
%%plot of weights agains difference scores
%plot(DATA,WEIGHTS,'.')

%The WSPCA assumes weights to be between zero and one
W=WEIGHTS/max(max(WEIGHTS));

%***********************************
%    PLOTS: WEIGHTS, SCREE
%***********************************

%Histogram of weights, to include in the manuscript
%Create Figure 2 in the paper: 2 lines are commented out because these
%require to download the fig.m and exportfig.m functions
%fig('units','centimeters','width',15,'height',15,'font','Arial')
hist(W(:),50)
axis([0 1 0 1.4*10^5])
xlabel('Weight')
ylabel('Frequency')
h=gcf();
exportfig(h,'../PlotHistogram.pdf','Resolution',300)

%***********************************
%    FIXING INPUT ARGUMENTS FOR WSPCA
%***********************************

%fixed model & algorithm parameters
NRSTARTS=10;       %number of starts in the multistart procedure
START='semi-rational'; %type of starting configuration: 'random', 'rational', or 'semi-rational'
MAXITER=200;     %maximum number of iterations
INIT=[]; %constraints
weakness=0.2;%for randomized lasso
SCALING='off';%'on' or 'off'
OFFSET='off';
LASSOTYPE='ordinary';
CONVERGENCE=1e-4;
HISTORY=1;
LASSOt=0;
nLambda=60; %number of values for tuning parameter!
orth=1;%1 or 0 for orthogonal or oblique

%SELECTING NUMBER OF COMPONENTS: SCREE PLOT FOR WPCA
R = 10;
[Ti,Pi,ci,si,Bi,Lossi]=WSPCAfast(DATA,W.^(0.5),R,0,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'rational');
LOSS=Lossi;
resultwPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si);
for nrstart=1:NRSTARTS;
    [Ti,Pi,ci,si,Bi,Lossi]=WSPCAfast(DATA,W.^(0.5),R,0,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'random');
    if Lossi<LOSS
        LOSS=Lossi;
        resultwPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si,'Loss',LOSS);
    end;
end;
VAF=[];
SUMVAF=[];
tVAF=[];
for r=1:R
    DATAhat_r=resultwPCA.T(:,r)*resultwPCA.P(:,r)';
    DATAhat_allr=resultwPCA.T(:,1:r)*resultwPCA.P(:,1:r)';
    VAF_r=1-sum(sum((DATA-DATAhat_r).^2))/(sum(sum(DATA.^2)));
    sumVAF_r=1-sum(sum((DATA-DATAhat_allr).^2))/(sum(sum(DATA.^2)));
    VAF=[VAF VAF_r];%note: orthogonal components so total VAF=sum(VAF)
    SUMVAF=[SUMVAF sumVAF_r];
end;
subplot(1,2,1)
bar(VAF)
xlabel('Component')
ylabel('Variance Accounted For')
subplot(1,2,2)
bar(SUMVAF)
xlabel('Component')
ylabel('Variance Accounted For')
h=gcf()
exportfig(h,'../PlotPVE.pdf','Resolution',300)

%results in R=3!
R=3;

%fit/pred.perf. for non-sparse (w)PCA
%PCA
Wunw = W~=0;
[Ti,Pi,ci,si,Bi,Lossi]=WSPCAfast(DATA,Wunw,R,0,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'rational');
LOSS=Lossi;
resultPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si);
XHAT=resultPCA.T*resultPCA.P';
RES=DATA-XHAT;
errspca=sum(sum(RES.^2))/sum(sum(DATA.^2));
corrcoef(resultPCA.T)
[b1,bint,r,rint,stats] = regress(y,[ones(I,1) resultPCA.T]);
rsq2008=stats(1)
yhat = [ones(I,1) resultPCA.T]*b1;
MSEtrain = sum((y-yhat).^2)/length(y)
%predictive quality for the continuous outcome
T2007=DATA2007*pinv(resultPCA.P')
y2007pred=[ones(I2007,1) T2007]*b1;
rsq2007=corr(y2007,y2007pred)^2
MSEtest = sum((y2007-y2007pred).^2)/length(y2007)
%wPCA
[Ti,Pi,ci,si,Bi,Lossi]=WSPCAfast(DATA,W.^(0.5),R,0,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'rational');
LOSS=Lossi;
resultwPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si);
for nrstart=1:NRSTARTS;
    [Ti,Pi,ci,si,Bi,Lossi]=WSPCAfast(DATA,W.^(0.5),R,0,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'random');
    if Lossi<LOSS
        LOSS=Lossi;
        resultwPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si,'Loss',LOSS);
    end;
end;
XHAT=resultwPCA.T*resultwPCA.P';
RES=DATA-XHAT;
errspca=sum(sum(RES.^2))/sum(sum(DATA.^2));
corrcoef(resultwPCA.T)
[b1,bint,r,rint,stats] = regress(y,[ones(I,1) resultwPCA.T]);
rsq2008=stats(1)
yhat = [ones(I,1) resultwPCA.T]*b1;
MSEtrain = sum((y-yhat).^2)/length(y)
%predictive quality for the continuous outcome
T2007=DATA2007*pinv(resultwPCA.P')
y2007pred=[ones(I2007,1) T2007]*b1;
rsq2007=corr(y2007,y2007pred)^2
MSEtest = sum((y2007-y2007pred).^2)/length(y2007)

%SAVE LOADINGS
[k ind] = sort(abs(resultwPCA.P(:)),'descend');
resultwPCA.P(abs(resultwPCA.P)<k(628)) = 0;
dlmwrite(['../DATA/wpcaloadings.txt'],resultwPCA.P,'delimiter','\t')


%Determine q, number of non-zero coefficients 
%Based on Meinshausen & Buhlmann (2010)
nsamples=250;%number of resamples
splitratio=0.9;%proportional sample size for resamples
EV=1;%error rate: expected number of false positives
pi_thr=0.9;%probability trhreshold
J=size(DATA,2);
qlambda=R*sqrt(EV*(2*pi_thr-1)*J);
CARD=qlambda;

%***********************************
%    PLAIN WSPCA ANALYSIS: CARD CONSTRAINT
%***********************************

%NO stability selection
[Ti,Pi,ci,si,Bi,Lossi]=WSPCAcardinality_overpcs(DATA,W.^(0.5),R,0,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'rational',round(qlambda));
LOSS=Lossi;
resultwSPCA2=struct('T',Ti,'P',Pi,'c',ci,'s',si);
%more straightforward in comparison with pma is rational start only
% for nrstart=1:NRSTARTS;
%     [Ti,Pi,ci,si,Bi,Lossi]=WSPCAcardinality_overpcs(DATA,W.^(0.5),R,0,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'semi-rational',round(qlambda));
%     if Lossi<LOSS
%         LOSS=Lossi;
%         resultwSPCA2=struct('T',Ti,'P',Pi,'c',ci,'s',si,'Loss',LOSS);
%     end;
% end;
XHAT=resultwSPCA2.T*resultwSPCA2.P';
RES=DATA-XHAT;
errspca=sum(sum(RES.^2))/sum(sum(DATA.^2));
corrcoef(resultwSPCA2.T)
[b1,bint,r,rint,stats] = regress(y,[ones(I,1) resultwSPCA2.T]);
rsq2008=stats(1)
yhat = [ones(I,1) resultwSPCA2.T]*b1;
MSEtrain = sum((y-yhat).^2)/length(y)
%predictive quality for the continuous outcome
T2007=DATA2007*pinv(resultwSPCA2.P')
y2007pred=[ones(I2007,1) T2007]*b1;
rsq2007=corr(y2007,y2007pred)^2
MSEtest = sum((y2007-y2007pred).^2)/length(y2007)

%***********************************
%     STABILITY SELECTION
%***********************************

%0. Initialize zero/non-zero status loading matrix with all zero values
SEL=zeros(J,R);
teller=1;
q=0;
%1. Sequence over CARD until nr of stable coefs > q
while q<qlambda
    SEL_old=SEL;%Usually we need the matrix BEFORE the while loop is interrupted,
    %this is SEL_old (in order to maintain the expected number of false
    %positives conservatively under control)
    %2. Calculate probability matrices for each of the CARD tuning
    CARD = round(qlambda*((100+25*teller)/100));
    [PrMAT]=subsamplingWSPCAcard(DATA,W.^(0.5),R,0,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'rational',nsamples,splitratio,CARD);
    SEL=SEL+(PrMAT>=pi_thr);
    q=sum(sum(SEL~=0))
    fname=['PrMAT_R1_N250card',num2str(teller),'.mat'];
    save([fname],'PrMAT')
    teller=teller+1;
end;
%obtain exactly qlambda selected variables: form the last run, those needed
%to attain this number are selected as those with highest prob.select.
k=find(SEL_old==0 & PrMAT>=pi_thr);
[x ind]=sort(PrMAT(k),'descend');
nrstillneeded = qlambda-sum(sum(SEL_old~=0));
SEL_old(k(ind(1:nrstillneeded)))=1;%because this gives the 
SEL=SEL_old;

% %*****
% %This code may be useful to restart somewhere in loop over lasso values or
% %to recalculate the final zero/non-zero status of the loadings
% SEL = zeros(J2007,R);
% for i=1:5
%     eval(['load PrMAT_N250card',num2str(i),'.mat'])
%     SEL=SEL+(PrMAT>=pi_thr);
%     q=sum(sum(SEL~=0));
%     teller=teller+1;
% end;

%3. Obtain estimates of non-zero coefficients
[Ti,Pi,ci,si,Bi,Lossi]=WSPCAfast(DATA,W.^(0.5),R,0,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],SEL,orth,'rational');
LOSS=Lossi;
resultwSPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si);
for nrstart=1:NRSTARTS;
    [Ti,Pi,ci,si,Bi,Lossi]=WSPCAfast(DATA,W.^(0.5),R,0,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],SEL,orth,START);
    if Lossi<LOSS
        LOSS=Lossi;
        resultwSPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si,'Loss',LOSS);
    end;
end;
XHAT=resultwSPCA.T*resultwSPCA.P';
RES=DATA-XHAT;
errspca=sum(sum(RES.^2))/sum(sum(DATA.^2));
corrcoef(resultwSPCA.T)

%4a. Calculate fit (2008 data) using continuous measure for vaccine
%efficacy
[b1,bint,r,rint,stats] = regress(y,[ones(I,1) resultwSPCA.T]);
rsq2008=stats(1)
yhat = [ones(I,1) resultwSPCA.T]*b1;
MSEtrain = sum((y-yhat).^2)/length(y)

%4b. Calculate prediction error
%estimate T for 2007 using the loadings from 2008 data
%load Rocke Lorenzato/Durbin variance parameters and the D0 and D3 data
load ../R/GSE29614/rawdata/GSE29614_RLparsD0.txt;
load ../R/GSE29614/rawdata/GSE29614_RLparsD3.txt;
load ../DATA/DATA_D0_2007
load ../DATA/DATA_D3_2007
%calculate weights according to Rocke Lorenzato model: this is difference
%of expression scores (D3-D0) so we use 4.20 of the Rocke & Durbin (2001) paper
%Because of background correction, we consider offsets to be equal to zero
epsilonD0=GSE29614_RLparsD0(2);
etaD0=GSE29614_RLparsD0(3);
epsilonD3=GSE29614_RLparsD3(2);
etaD3=GSE29614_RLparsD3(3);
varlnexpr=etaD0^2+etaD3^2+(epsilonD0^2*(DATA_D0).^(-1))+(epsilonD3^2*(DATA_D3).^(-1));
WEIGHTS=(1./varlnexpr);
WEIGHTS=WEIGHTS';
W2007=WEIGHTS;%weights for 2007 data yet maybe it makes more sense to use an unweighted approach
%see here below
for nrstart=1:NRSTARTS;
    [T,P,c,s,B,Loss]=WSPCAnoupdateP(DATA2007,W2007.^(0.5),R,1e-9,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],resultwSPCA.P(1:54675,:),[],orth,START);
    if nrstart==1
        LOSS=Loss;
        result2007=struct('T',T,'P',P,'c',c,'s',s);
    elseif Loss<LOSS
        LOSS=Loss;
        result2007=struct('T',T,'P',P,'c',c,'s',s,'Loss',LOSS);
    end;
end;
T2007=result2007.T;
T2007=DATA2007*pinv(resultwSPCA.P')

%predictive quality for the continuous outcome
y2007pred=[ones(I2007,1) T2007]*b1;
rsq2007=corr(y2007,y2007pred)^2
MSEtest = sum((y2007-y2007pred).^2)/length(y2007)
%predictive quality for the binary outcome
yhatbin2007 = predict(mdl, T2007)
1-sum((round(yhatbin2007)-ybin2007).^2)/length(ybin2007)

%SAVE LOADINGS
dlmwrite(['../DATA/wspcaloadings.txt'],resultwSPCA.P,'delimiter','\t')


%******************
%    ANNOTATION
%******************

%make rnk files / files with GENEIDS for GSEA / AmiGO
fid=fopen('../DATA/annotation.txt');
AFFYID_all=textscan(fid,'%s%s%s');
AFFYID=AFFYID_all{1};
GENEID=AFFYID_all{2};
ID=AFFYID_all{3};
fclose(fid);

%check gene names of selected genes and save their GENEID
for r=1:R
    geneids=GENEID(resultwSPCA.P(:,r)~=0);
    affies=AFFYID(W(:,r)~=0);
    fid=fopen(['../RESULTS/GENEIDS_WSPCA_pc',num2str(r),'.txt'],'w');
    for row = 1:length(geneids)
        %fprintf(fileID,formatSpec,C{row,:});
        fprintf(fid,'%s \n',geneids{row,:})
    end;
    fclose(fid)
end;

%Also make files for annotation of the loadings resulting from PMD
load '../RESULTS/PMDloadings.txt'
for r=1:R
    geneids=GENEID(PMDloadings(:,r)~=0);
    affies=AFFYID(W(:,r)~=0);
    fid=fopen(['../RESULTS/GENEIDS_PMD_pc',num2str(r),'.txt'],'w');
    for row = 1:length(geneids)
        %fprintf(fileID,formatSpec,C{row,:});
        fprintf(fid,'%s \n',geneids{row,:})
    end;
    fclose(fid)
end;
