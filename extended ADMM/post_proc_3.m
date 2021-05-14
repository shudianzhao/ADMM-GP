function [y,LB] = post_proc_3(Z,At,Bt, C,b,f, S);
%PostPROCESSING: determine a feasible y, given Z and S from dual
% SDP with ineq.constraints.

% (P) min <C,X> s.t. A*X(:) = b;B*X(:) = s; X psd,non-negative, s>=f
% (D) max b'y s.t. C - A^t(y)-B^t(y_bar) = Z + S ; v=y_bar;Z psd, S
%     non-negative,

n = size(C,1);
n2 = n*n;
if n2==size(At,2);
    At= At';
end

if n2==size(Bt,2);
    Bt= Bt';
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
M = C - Zplus;

if ~isempty(S);
     [y,LB1]= linprog([-b;-f],[At,Bt],M(:),[],[],[-inf*ones(length(b),1);zeros(length(f),1)],[inf*ones(length([b;f]),1)]);
else
    [y,LB1]= linprog([-b;-f],[],[],[At,Bt],M(:),[-inf*ones(length(b),1);zeros(length(f),1)]);
end
LB= -LB1;

fprintf('LP safe lower bound LBnew = %10.9e  \n',LB);
