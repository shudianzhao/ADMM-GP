function [At, b, Ac] = diag_constr(n);
% generates At and b, such that At'(X)=b is diag(X)=e 
% input: n...dimension of X
% call: [At, b, Ac] = diag_constr(n);
 
% 08/06/04

Ac = cell(n,1);  
for i=1:n
  B = sparse(i,i,1,n,n);
  At(:,i) = B(:);
  Ac{i} = B;
end;

b = ones(n,1);
