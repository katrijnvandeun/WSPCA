function RES=residual(DATA,c,s,T,P)
%RESIDUAL calculates the difference between the (transformed) observed and modeled scores
%function: [(DATA-1c')S-TP']
%K. Van Deun, Dept. Psychology, KU Leuven
%version 1: March 2013
[I,J]=size(DATA);

TDATA=(DATA-(ones(I,1)*c')).*(ones(I,1)*s');
XHAT=T*P';
RES=TDATA-XHAT;