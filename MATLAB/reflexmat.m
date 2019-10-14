function MAT=reflexmat(m)
%REFLEXMAT function to construct matrix of reflections
MAT=ones(1,m);
for i=1:m-1
    B=nchoosek(1:m,i);
    for j=1:size(B,1)
        v=ones(1,m);
        v(B(j,:))=-1;
        MAT=[MAT; v];
    end;
end;