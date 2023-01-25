function ret = SimHPN_isSolution(v,A)
% returns 1 if the contPN system defined by the matrix A and the output
% given by the places in v is observable

S = zeros(length(v),size(A,1));
for i = 1 : length(v)
    S(i,v(i)) = 1;
end
ob = obsv(A,S);
if (rank(ob) == size(A,2))
    ret = 1;
else
    ret = 0;
end
return;