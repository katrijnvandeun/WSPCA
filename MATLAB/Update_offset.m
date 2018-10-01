function offsets =Update_offset(W,DATA,c,s,T,P)
%UPDATE_OFFSET Calculates the update of the offsets in WSPCA
%   Objective: minimize L(c)=||W°((X-1c')S-TP^T)||
%   The WLS problem is solved by solving the problem for a majorizing
%   function; the problem than becomes one of iteratively solving a
%   generalized (Penrose) regression problem with closed form solution

[I,J]=size(DATA);

RES=residual(DATA,c,s,T,P);
XHAT=ones(I,1)*(c'.*s');
YSTAR=ystar(XHAT,RES,W);
offsets=(sum(YSTAR).*(s.^(-1))'*(I^(-1)))';
offsets=(mean(YSTAR))';
end

