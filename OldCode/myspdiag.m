function [result] = myspdiag(subdia,diag,supdia)

K=length(diag)+2;
lineindex = [[1:K-2],[1:K-2],[1:K-2]];
colindex = [[1:K-2],[2:K-1],[3:K]];
data = [subdia,diag,supdia];

result=sparse(lineindex,colindex,data,K-2,K,3*K);