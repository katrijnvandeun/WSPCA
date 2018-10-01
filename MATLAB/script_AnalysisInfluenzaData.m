%last update: 27th August 2018

%load DATA (difference scores)
load ../DATA/DATA2008
DATA=DATA_diffD3';
load ../DATA/DATA2007
DATA2007=DATA_diffD3';
[I J]=size(DATA);
[I2007 J2007]=size(DATA2007);

%load Rocke Lorenzato/Durbin variance parameters and D0 and D3 data
load ../R/GSE29617/rawdata/GSE29617_RLparsD0.txt;
load ../R/GSE29617/rawdata/GSE29617_RLparsD3.txt;
load ../DATA/DATA_D0
load ../DATA/DATA_D3
%Calculate WEIGHTS according to Rocke Lorenzato model: this is difference
%of expression scores (D3-D0) so we use expression (4.20) of the
%Rocke & Durbin (2001) paper. Because of background correction, we consider
%the offsets to be equal to zero
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

%The WSPCA assumes weights to be between zero and one
W=WEIGHTS/max(max(WEIGHTS));

%Histogram of weights, Creates Figure 2 in the paper.
%Export of pdf and sizing requires to download the fig.m and exportfig.m functions
%fig('units','centimeters','width',15,'height',15,'font','Arial')
hist(W(:),50)
axis([0 1 0 1.4*10^5])
xlabel('Weight')
ylabel('Frequency')
h=gcf();
%exportfig(h,'../PlotHistogram.pdf','Resolution',300)

%Load Criterion variable
load ../DATA/TIVtiter2008
y=log2(Y(:,4));%4th column is max of three titer types
load ../DATA/TIVtiter2007
y2007=log2(Y(:,4));%4th column is max of three titer types

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
LASSOt=0; %penalty on component scores
nLambda=60; %number of values for tuning parameter!
orth=1;%1 or 0 for orthogonal or oblique component scores respectively

%SELECTING NUMBER OF COMPONENTS: SCREE PLOT FOR WPCA
R = 10;
VAF=[];
SUMVAF=[];
for r=1:R
    [Ti,Pi,ci,si,Bi,Lossi]=WSPCA(DATA,W.^(0.5),r,0,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'rational',[]);
    LOSS=Lossi;
    resultwPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si);
    for nrstart=1:NRSTARTS;
        [Ti,Pi,ci,si,Bi,Lossi]=WSPCA(DATA,W.^(0.5),r,0,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'random',[]);
        if Lossi<LOSS
            LOSS=Lossi;
            resultwPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si,'Loss',LOSS);
        end;
    end;
    DATAhat_r=resultwPCA.T(:,r)*resultwPCA.P(:,r)';
    DATAhat_allr=resultwPCA.T(:,1:r)*resultwPCA.P(:,1:r)';
    VAF_r=1-sum(sum((DATA-DATAhat_r).^2))/(sum(sum(DATA.^2)));
    sumVAF_r=1-sum(sum((DATA-DATAhat_allr).^2))/(sum(sum(DATA.^2)));
    VAF=[VAF VAF_r];%note: components are not nested!
    SUMVAF=[SUMVAF sumVAF_r];
end;
subplot(1,2,1)
bar(SUMVAF-[0 SUMVAF(1:end-1)])%VAF per component
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

%***********************************
%     SET PARAMETERS SUBSAMPLING
%***********************************
%Based on Meinshausen & Buhlmann (2010)
nsamples=250;%number of resamples
splitratio=0.9;%proportional sample size for resamples
EV=1;%error rate: expected number of false positives
pi_thr=0.9;%probability trhreshold
%*****

%STABILITY SELECTION: Based on cardinality constraint
%1.Determine q, number of non-zero coefficients
J=size(DATA,2);
qlambda=R*sqrt(EV*(2*pi_thr-1)*J);
CARD=qlambda;
%2. Sequence over CARD until nr of stable coefs > q
%Initialize zero/non-zero status loading matrix with all zero values
SEL=zeros(J,R);
%2. Loop over lambda values, from qlambda until the number of
%non-zero coefficients is >= qlambda
teller=0;
q=0;
while q<qlambda
    SEL_old=SEL;%we need the matrix BEFORE the while loop is interrupted,
    %this is SEL_old (in order to maintain the expected number of false
    %positives conservatively under control)
    %2. Calculate probability matrices for each of the lasso tuning
    %parameter values
    %CARD = round(qlambda*((100+25*teller)/100))%CARD should be in range 1200-1400
    CARD = round(1100*((100+2*teller)/100))%CARD should start in range 1100-1300
    %at the end of the loop CARD was equal to 1298
    [PrMAT]=subsamplingWSPCAcard(DATA,W.^(0.5),R,0,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],[],orth,'rational',nsamples,splitratio,CARD);
    SEL=SEL+(PrMAT>=pi_thr);
    q=sum(sum(SEL~=0))
    teller=teller+1;
end;

%3. Obtain estimates of non-zero coefficients
[Ti,Pi,ci,si,Bi,Lossi]=WSPCA(DATA,W.^(0.5),R,0,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],SEL_old,orth,'rational',[]);
LOSS=Lossi;
resultwSPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si);
for nrstart=1:NRSTARTS;
    [Ti,Pi,ci,si,Bi,Lossi]=WSPCA(DATA,W.^(0.5),R,0,OFFSET,SCALING,LASSOTYPE,MAXITER,CONVERGENCE,HISTORY,LASSOt,[],SEL_old,orth,START,[]);
    if Lossi<LOSS
        LOSS=Lossi;
        resultwSPCA=struct('T',Ti,'P',Pi,'c',ci,'s',si,'Loss',LOSS);
    end;
end;
XHAT=resultwSPCA.T*resultwSPCA.P';
RES=DATA-XHAT;
errspca=sum(sum(RES.^2))/sum(sum(DATA.^2));
corrcoef(resultwPCA.T)

%4a. Calculate fit (2008 data) using continuous measure for vaccine
%efficacy
[b1,bint,r,rint,stats] = regress(y,[ones(I,1) resultwSPCA.T]);
rsq2008=stats(1)

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

%predictive quality for the continuous outcome
y2007pred=[ones(I2007,1) T2007]*b1;
rsq2007=corr(y2007,y2007pred)^2

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
        fprintf(fid,'%s \n',geneids{row,:})
    end;
    fclose(fid)
end;
