%Script to generate data according to sparse bilinear model with additive and
%multiplicative error and with missing data

%set seed
rng(2010)

%Define factors and levels
propmissing_v=[0];% 0.2 0.4 0.6];
sigma_mult=[0.05 0.01 0];
sigma_add=sigma_mult;
propsparse_v=[0 0.5];
I_v=[100 1000];
J_v=[100 1000];
R_v=[2];
NRREPLICS=10;


RESULT=[];
teller=1;
for nrpc=1:length(R_v)
    R=R_v(nrpc);
    for i=1:length(I_v)
        I=I_v(i);
        for j=1:length(J_v)
            J=J_v(j);
            for pr=1:length(propmissing_v)
                propmissing=propmissing_v(pr);
                for prsp=1:length(propsparse_v)
                    propsparse=propsparse_v(prsp);
                    for noise_mult=1:length(sigma_mult)
                        noiselevel_mult=sigma_mult(noise_mult);
                        for noise_add=1:length(sigma_add)
                            noiselevel_add=sigma_add(noise_add);
                            for r=1:NRREPLICS
                                %Generation TRUE structure \mu=TP'
                                X=randn(I,J);
                                Xc=X-ones(I,1)*(mean(X));
                                [T s v]=svds(Xc,R);%orthogonal T
                                %T=abs(T);
                                P=sqrt(I)*rand(J,R);
                                P=abs(P);
                                %P=v*s;
                                v=randperm(R*J,round(R*J*propsparse));%P sparse
                                P(v)=0;
                                TRUEpca=T*P';
                                %rescale to standard normal distribution
                                ss=sum(sum(TRUEpca.^2));
                                TRUEpca=sqrt(I*J-1)*TRUEpca/(sqrt(ss));
                                TRUE=exp(TRUEpca);
                                vartrue=var(TRUE(:));
                                addNOISE=sqrt(vartrue*noiselevel_add)*randn(I,J);
                                multNOISE=sqrt(vartrue*noiselevel_mult)*randn(I,J);
                                X=(TRUE.*exp(multNOISE))+addNOISE;
                                %find values < eps (log transform is taken)
                                f=find(X<0);
                                X(f)=1;
                                DATA=log(X);
                                musq=TRUE.^2;
                                if noiselevel_mult~=0 | noiselevel_add~=0
                                    W=musq./((musq*noiselevel_mult*vartrue)+(noiselevel_add*vartrue));
                                else W=ones(I,J);
                                end;
                                %create missing values
                                %v=randperm(I*J,round(I*J*propmissing));
                                W(f)=0;
                                %W(v)=0;
                                f2=find(isinf(abs(DATA)));
                                DATA(f2)=0;
                                W(f2)=0;
                                nrnumeric = prod(size(f))+prod(size(f2));
                                loop=1;
                                while(nrnumeric>0.1*I*J)
                                    loop=loop+1;
                                    X=randn(I,J);
                                    Xc=X-ones(I,1)*(mean(X));
                                    [T s v]=svds(Xc,R);%orthogonal T
                                    %T=abs(T);
                                    P=sqrt(I)*randn(J,R);
                                    P=abs(P);
                                    %P=v*s;
                                    v=randperm(R*J,round(R*J*propsparse));%P sparse
                                    P(v)=0;
                                    TRUEpca=T*P';
                                    %rescale to standard normal distribution
                                    ss=sum(sum(TRUEpca.^2));
                                    TRUEpca=sqrt(I*J-1)*TRUEpca/(sqrt(ss));
                                    TRUE=exp(TRUEpca);
                                    vartrue=var(TRUE(:));
                                    addNOISE=sqrt(vartrue*noiselevel_add)*randn(I,J);
                                    multNOISE=sqrt(vartrue*noiselevel_mult)*randn(I,J);
                                    X=(TRUE.*exp(multNOISE))+addNOISE;
                                    %find values < eps (log transform is taken)
                                    f=find(X<0);
                                    X(f)=1;
                                    DATA=log(X);
                                    musq=TRUE.^2;
                                    if noiselevel_mult~=0 | noiselevel_add~=0
                                        W=musq./((musq*noiselevel_mult*vartrue)+(noiselevel_add*vartrue));
                                    else W=ones(I,J);
                                    end;
                                    W(f)=0;
                                    f2=find(isinf(abs(DATA)));
                                    DATA(f2)=0;
                                    W(f2)=0;
                                    nrnumeric = prod(size(f))+prod(size(f2));
                                    if loop>25
                                        loop
                                        teller
                                        nrnumeric=0;
                                    end;
                                end;
                                fname=['../DATA/DATAnomissc',num2str(teller),'.dat'];
                                dlmwrite(fname,DATA,'delimiter','\t')
                                fname=['../DATA/WEIGHTSnomissc',num2str(teller),'.dat'];
                                dlmwrite(fname,W,'delimiter','\t')
                                fname=['../DATA/TRUEnomissc',num2str(teller),'.dat'];
                                dlmwrite(fname,TRUEpca,'delimiter','\t')
                                fname=['../DATA/TRUEPnomissc',num2str(teller),'.dat'];
                                dlmwrite(fname,P,'delimiter','\t')
                                fname=['../DATA/TRUETnomissc',num2str(teller),'.dat'];
                                dlmwrite(fname,T,'delimiter','\t')
                                setting=struct('R',R,'PrSparse',propsparse,'AddNoise',noiselevel_add,'MultNoise',noiselevel_mult,...
                                    'PrMiss',propmissing,'AddNoiseVar',vartrue*noiselevel_add,'MultNoiseVar',vartrue*noiselevel_mult);
                                fname=['../DATA/STRUCTnomissc',num2str(teller),'.mat'];
                                save([fname],'setting');
                                teller=teller+1;
                            end;
                        end;
                    end;
                end;
            end;
        end;
    end;
end;