function [ Y,y,y_bar,Z,S,s,v, sigma,iter, secs, dual,primal, R,Rt] = mprw_ineq_general( A, B, b, C, f,L,U,max_iter, sigma, tol, Y, Z,v,s,R, Rt)
% simple "regularization method" : 
% (based on primal prox, or dual augmented Lagrangien)
% corresponds to boundary point method, and to block-coordinate descent

% solves: min <C,X> s.t. A*X(:) = b;B*X(:) = s; X psd,non-negative, s>=f
%         max b'y s.t. C - A^t(y)-B^t(y_bar) = Z + S ; v=y_bar;Z psd, S
%         non-negative,
% max_iter: max number of iterations
% sigma: penalty parameter for augmented Lagrangian 
%        a typical value (b and C normalized) is between .1 and 10
% tol: stopping condition
%      we stop once primal and dual relative infeasibility < tol
% 
% simplest call: [ X, y,y_bar, Z,S, s,v] = mprw_ineq_general( A, B, b, C ,f);
% call:[ X,y,y_bar,Z,S,s,v, sigma,R,Rt]   = mprw2_ineq( A, B, b, C, f, max_iter, sigma, tol, X, Z, S,v,s,R, Rt);
% optional input (for restart):  sigma, tol, X, Z, S

% no scaling

tstart = cputime; 
m=min(size(A));
p=min(size(B));
n2=max(size(A));
n = sqrt(n2);
if length(b)~=m; error('Size A mismatches b.'); end
if (min(size(C)) == 1)&&(max(size(C)) == n2); C=reshape(C,n,n); end
if (size(C,1)~=n)||(size(C,2)~=n); error('Size C mismatches A.'); end

if size(A,1) == m ;
  At = A';
else
  At = A; A = At';
end

if size(B,1) == p ;
  Bt = B';
else
  Bt = B; B = Bt';
end




% form Q' (and cholesky, if not part of input) 


AAT = A*At; 
BBT = B*Bt;

Q = [AAT, A*Bt ; B*At, BBT+speye(p)];


if nargin <= 13;
   [R,q] = chol( Q);
   Rt = R';    % only for speed-up
   secs = cputime - tstart;
   fprintf(' secs after chol:   %12.5f \n', secs);
   if q>0; 
         fprintf(' rows of A linearly dependent. q = %4.0d, m + p = %4.0d\n',q,m+p);
          error(' nothing done');
   end
end


% sinon
normb = norm(b);
normC = sqrt(C(:)'*C(:));

% initialize
if nargin == 8; max_iter=100; sigma=0.2; tol=1e-5; end
if nargin == 9; sigma = 0.2; tol = 1e-5; end
if nargin == 10; tol = 1e-5; end
if nargin == 11; Y = zeros(n,n); Z =Y; end
if nargin <= 13;  S= ones(n,n); s = zeros(p,1); v = s;Y = zeros(n,n); Z =Y; end;
  
% outer stopping condition
iter = 1;              % iteration count

% g2 = f - B* Y(:);
done = 0;              % stopping condition 
update_it = 1;        % update sigma every update_it iterations 

fprintf(' it     secs       dual         primal     lg(rrd)    lg(rrp)   sigma\n');


g1 = b - A*Y(:);          % primal residue 1
g2 = s - B*Y(:);


t_min = 1e-4;
t_max = 1e+4;
% start outer iterations
while done == 0;       % while stopping conditions not satisfies
  
  weight = 2^(-(iter-1)/100);
% w = (y, y_bar)
% given Y,Z,S  and sigma, solve for y and y_bar
  rhs1  = A*(C(:) - Z(:) - S(:)) + g1/sigma;
  rhs2  = B*(C(:) - Z(:) - S(:)) + g2/sigma + v;
  rhs = [rhs1 ; rhs2];
  
  w_tmp = Rt\rhs; w = R\w_tmp;
  
  y = w(1:m);
  y_bar = w((m+1):end);


% now form M0= M0(y,y_bar,Y,sigma,Z)
  Aty = reshape(At*y, n, n);    % A^T(y)
  Bty_bar = reshape(Bt*y_bar, n ,n);

  
% now form M1= M1(y,y_bar,Y,sigma,Z) & M2=M2(v,y,s,sigma)
  M0 = Aty + Bty_bar + Y/sigma + Z - C;
  M0 = (M0 + M0')/2;      % should be symmetric

  
  M1 =  -y_bar + s/sigma;
  
  S_hat = M0*sigma;
  
%   S_hat(S_hat<=0)=0; %projf(S,L,U)
   S_hat = projf(S_hat,L,U);
  
  S = S_hat/sigma - M0;
  
  s = sigma*M1;
  
  s = projf(s,f,Inf);
  
  v=  s/sigma - M1;
  
% now compute eigenvalue decomposition of M2
  M2 = Aty + Bty_bar + Y/sigma + S - C;
  M2 = (M2 + M2')/2;      % should be symmetric
   
  
 
  
  
    
  [ev, lam] = eig(M2); lam = diag(lam);

% compute projection Mp and Mn of M onto PSD and NSD
  I = find(lam> 0); j = length(I);
  if j < n/2;
     evp = zeros( n,j);
     for r=1:j;
         ic = I(r); evp(:,r) = ev(:,ic)*sqrt(lam(ic)); 
     end
     if j==0; evp = zeros(n,1); end
     evpt= evp';
     Mp = evp*evpt;     
     Mn = M2 - Mp;
  else
     I = find(lam < 0); j = length( I); % should be <= n/2
     evn = zeros( n,j);
     for r=1:j;
        ic = I(r); evn(:,r) = ev(:,ic)*sqrt( -lam(ic)); 
     end
     if j==0; evn = zeros(n,1); end
     evnt= evn';
     Mn = -evn*evnt;
     Mp = M2- Mn; 
  end

  Mn = (Mn+Mn')/2; Mp = (Mp+Mp')/2;  % these should be symm.

  Z = -Mn; X = sigma * Mp;
  
 

  

% update Y
%   Y = (1-1.6)*Y + 1.6*X;
  
  Y = X;  
  g1 = b- A * Y(:); 
  g2 = s- B * Y(:);
  G1 = -Z -S + C - Aty - Bty_bar;
  G2 = v - y_bar;

    
  
% some output

err_d = norm(G1,'fro');

err_b1 = norm(Y-projf(Y-S,0,inf));

err_b2 = norm(s-projf(s-v,f,inf));


dual = b'*y +s'*v + S(:)'*Y(:); 
secs = cputime-tstart;
err_p = norm(g1);   
primal = C(:)'*Y(:);

rel_err_p = err_p/(1+normb)+norm(g2)/(1+norm(s)); rel_err_d = err_d/(1+normC)+norm(G2)/(1+norm(v));


rel_err_b1 = err_b1/(1+norm(Y)+norm(S));
rel_err_b2 = err_b2/(1+norm(s)+norm(v));
iter = iter+ 1;

if (mod(iter,50)==0);
fprintf( '%3.0d %8.2f %13.5e %13.5e %8.3f %8.3f  %9.6f\n', iter, ...
   secs,  dual, primal,log10(rel_err_d), log10(rel_err_p), sigma );
end

 
 
 
% check stopping conditions

if max([rel_err_d, rel_err_p,rel_err_b1,rel_err_b2]) < tol || iter > max_iter ; 
   fprintf( '%3.0d %8.2f %13.5e %13.5e %8.3f %8.3f  %9.6f\n', iter, ...
   secs, dual, primal, log10(rel_err_d),log10(rel_err_p), sigma );
   fprintf('total time: %10.3f \n', secs); 
  if (iter>max_iter);  
    fprintf('max outer iterations reached. \n');
  end
    done=1;
end  





% check for reduction of sigma
ratio = rel_err_p/rel_err_d;  
const=1.2; % 1.2
% 
%adjust stepsize
%  ratio = norm(X(:))/norm(Z(:));
 % sigma = (1-weight)*sigma + weight*projf(ratio,t_min,t_max) ;



% %increase sigma if p-error is less ; decrease sigma if p-error is large
 if (ratio < 1) && (mod(iter,update_it)==0) %0.2
     sigma = min(1e3,const*sigma);
 elseif (ratio > 10)  % 5
    sigma = max(1e-2,sigma/const);
 end

end
