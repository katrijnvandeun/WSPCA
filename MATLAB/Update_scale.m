function [scale P]=Update_scale(W,DATA,c,s,T,P)
%UPDATE_SCALE Calculates the update of the scale in WSPCA
%   Objective: minimize L(S)=||W°((X-1c')S-TP^T)|| s.t. S nonnegative and
%   diagonal.
%   The WLS problem is solved by solving the problem for a majorizing
%   function; the problem than becomes one of solving a constrained 
%   least-squares problem 
%   
%Author: K. Van Deun, APRIL 2013

[I,J]=size(DATA);
minspread=1/(10*J);

A=(DATA-ones(I,1)*c');
sumasq=sum(A.^2);
wmax=max(sumasq);
XHAT=(DATA-ones(I,1)*c')*diag(s);
RES=(T*P')-XHAT;
YSTAR=ystar(XHAT,RES,W);
dy=sum(A.*YSTAR);
%majorizing to obtain right metric; see Van Deun, 2005
q=s+(wmax^(-1))*(dy'-s.*sumasq');
q=q/norm(q);
k=find(q<0);
q(k)=abs(q(k));
P(k,:)=-1*P(k,:);
k=find(q<minspread);
q(k)=minspread;
scale=sqrt(J)*q;
end

