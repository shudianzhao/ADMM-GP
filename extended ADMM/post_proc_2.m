function [y,LB] = post_proc_2(Z,  At, C,b, S);
%PostPROCESSING: determine a feasible y, given Z and S

% dual SDP: max b'y s.t. Aty + Z + S = C;
%  S nonnegative

n = size(C,1);
n2 = n*n;
if n2==size(At,2);
    At= At';
end

eigtmp = eig(full(Z)); 
idx = find(eigtmp < 0); 
numneg = length(idx); 
if (numneg) 
   mineig = min(eigtmp(idx));
   fprintf('\n numneg = %3.0d,  mineigZnew = %- 3.2e',numneg,mineig);
    Zplus= project_W(Z);
else 
    fprintf('\n Z is PSD.');
    Zplus = Z;
end
%RHS
% M = C - Zplus;
M = C- Z;

if ~isempty(S);
    [y,LB1]= linprog(-b,At,M(:));
else
    [y,LB1]= linprog(-b,[],[],At,M(:));
end
LB= -LB1;

fprintf('LP safe lower bound LBnew = %10.9e  \n',LB);
