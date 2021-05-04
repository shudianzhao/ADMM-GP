function [ b, T, hash, g_new,brk] = separation3_x( x, b, T, hash, max_ineq, viol_tol);
%find new violated triangular constraints at x
% -x_ij - x_ik +x_jk >=-1  i, j, k \in [n],
%  x psd, x_ij \in [0,1]
% b,T,hash: data structure of existing inequalities
% T: Triangle constrains (index of nodes and the type of inequlities) b: rhs, constants 1
% gamma: violation value, gamma =  -(- xij - xik + xjk)-1, gamma = b-A(x)
% g_new: initial variables for gamma at new constraints (set =0)  
% new_ineq: minimum number of new inequalities is 

% call:  [ b, T, hash, gamma_new,brk] = separation3( x, b, T, hash,max_ienq,viol_tol);



  

  n = size(x,1);    
  brk = 0;
  


  y = reshape( x,n^2, 1);   % y vector for X
  h_tmp = sort( hash);        % hash function sorted

 [Tn, ~, hashn, g_new] = tri_sep_kc( y, n, max_ineq, h_tmp);

  m = length(hashn);         % number of new ineq.
  
  if m == 0;
      disp('nothing new violated.'); 
      g_new = []; brk=1; 
      return    % return with old data structure
  end;
  if m > 0;
    Tn = reshape( Tn, m, 4);
    bn = -ones( m, 1);      % right hand side must be -1
  end
% merge new with old constraints
  b = [ b; bn];             % update rhs vector
  hash = [hash; hashn];     % update hash function
  T = [T; Tn];
  
  if min(g_new) < 0; % b-A(X)<)
      disp('error separation.');
  end
  if max(g_new) < viol_tol;
      disp('violation is small');
      brk = 1;
  end 
  
  
%   ax = kc_a(n,m,Tn(:),x(:));  % A(X)>b is the tri. ineq
%   if size(T,1) <= 40
%    gamma_new = bn-ax;      % initialize new dual variables % gamma_new is violation A(X)-b 
%    if min(gamma_new)<0; disp('error separation.'); 
%       keyboard;
%    end
%   if max(gamma_new) < viol_tol;
%     disp('violation is small');
%     brk = 1;
%   end  
%   end
  