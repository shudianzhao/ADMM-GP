function [A,B,b,f,C]=make_GPKC(Ac,a,W);    
% form the GPKC problem
% min <C,X> s.t. A(X)=b, B(X)=f;
% Ac: Adjacency matrix
% W: capacity weight
% a: vertices weights

  n = size(Ac,1);
  n2 = n*n;
  L =  diag(Ac*ones(n,1))-Ac;
  C = L/2;
  [At,b,Ec]= diag_constr(n);
  A = At';

  B = sparse(n,n2);
  f = -W*ones(n,1);
  for i=1:n
     Btemp = zeros(n,n);
     Btemp(i,:) = -a'/2;
     Btemp(:,i)=-a'/2;
     Btemp(i,i)=- a(i);
     B(i,:) = Btemp(:)';
  end
end

%%
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
end

%%
