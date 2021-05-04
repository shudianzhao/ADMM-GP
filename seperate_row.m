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

