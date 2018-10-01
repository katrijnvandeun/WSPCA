function loss=WSPCALOSS(DATA,W,lasso,T,P,B,c,s,lassoT,B_T)
%WSPCALOSS calculates the value of the weighted sparse PCA function
%||W°((X-1c')S-TP')||²+lambda|P|_1
%K. Van Deun, Dept. Psychology, KU Leuven
%March 2013

[I,J]=size(DATA);

TDATA=(DATA-ones(I,1)*c').*(ones(I,1)*s');
Xhat=T*P';
PENALTY=sum(sum(B.*abs(P)));
PENALTYT=sum(sum(B_T.*abs(T)));
RES=TDATA-Xhat;
loss=sum(sum(((W.*RES).^2)))+lasso*PENALTY+lassoT*PENALTYT;