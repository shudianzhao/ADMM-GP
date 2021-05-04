function [A,b,C] = kequi_form(Ac,k);
% return the SDP relaxation for k-equipartition
% min <C,X> s.t. A(X) = b
% Ac: adjacency matrix
% k: partition number

[n1,n1] = size(Ac);
L =  diag(Ac*ones(n1,1))-Ac;
C= L./2;

m = idivide(int16(n1), k,'ceil') ;
[Et,b1,Ec]= diag_constr(n1);
[Et2, b2, Ec2] = seperate_row(n1,m);
A = [Et , Et2 ];
b = [b1;b2];
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
function [At, b, Ac] = seperate_row(n,m);

% generates At and b, such that At'(X)=m*e_n i.e Xen=men (sum of row i = m)
% input: n...dimension of X
% call: [At, b, Ac] = seperate_row(n,m);

%31/01/19
% At = sparse(n*n,0);
Ac = cell(n,1);  
m = double(m);

for i=1:n
    B = sparse(i,1:n,1,n,n);
    B = (B + B')/2;
    At(:,i) = B(:);
    Ac{i} = B;
end

b = m*ones(n,1);
end