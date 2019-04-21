%load expression, titers, and pheno data
%expression data
load ../R/GSE29617/rawdata/GSE29617_rma.txt;
DATA=GSE29617_rma;

%antibody titers
load ../DATA/TIVtiters_2008.txt
TITERS=TIVtiters_2008;

%pheno data:
fid=fopen('../R/GSE29617/rawdata/GSE29617_pdata.txt');
allcovs=textscan(fid,'%s%s','Delimiter','\t','EndOfLine','\n');
subjID=allcovs{1};
timeID=allcovs{2};
fclose(fid);

% %make two blocks of data
%pre-vaccination D0
D0=strmatch('time point: D0',timeID);
Subj_D0=subjID(D0,:);
%three days after D3
D3=strmatch('time point: D3',timeID);
Subj_D3=subjID(D3,:);

%find subjects with data on D0 and D
matchD0=[];
matchD3=[];
for i=1:length(D3)
    match0=strmatch(Subj_D3(i,:),Subj_D0,'exact');
    matchD0=[matchD0;D0(match0,:)];
    matchD3=[matchD3;D3(i,:)];
end;
%CHECK!
[subjID(matchD0,:) subjID(matchD3,:)] % -> OK
[timeID(matchD0,:)  timeID(matchD3,:)] % -> OK

DATA_D0=DATA(:,matchD0);
DATA_D3=DATA(:,matchD3);
TITERS_D3=TITERS(:,matchD3);

%Obtain difference scores wrt baseline (D0)
DATA_diffD3=DATA_D3-DATA_D0;

%pre-processing: centering + scaling to unit SS per gene
DATA_diffD3_std=(STD(DATA_diffD3')');

%pre-processing of titers
m1=TITERS_D3(2,:)./TITERS_D3(1,:);
m2=TITERS_D3(4,:)./TITERS_D3(3,:);
m3=TITERS_D3(6,:)./TITERS_D3(5,:);
Y=[m1' m2' m3' max([m1' m2' m3'],[],2)]

%Write data to data folder
save ../DATA/DATA2008_std DATA_diffD3_std;
save ../DATA/DATA2008 DATA_diffD3;
save ../DATA/DATA2008_std.txt DATA_diffD3_std -double -ascii;
save ../DATA/DATA2008.txt DATA_diffD3 -double -ascii;
save ../DATA/DATA_D0 DATA_D0;
save ../DATA/DATA_D3 DATA_D3;
save ../DATA/TIVtiter2008 Y
save ../DATA/TIVtiter2008.txt Y -double -ascii;