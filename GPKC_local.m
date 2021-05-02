function  [part_cell_final,part_final,newX_final,obj]= GPKC_local(L,W,a,part_cell,part,X,obj);

% Input: a kequiparition
% output: a new kequipairtion with better objective function value
% Objective function: min <L/2, X>, X {0,1}
% L: adjency matrix
% k: the partitioning number

fprintf('Best obj before local adjustment: %6.0d \n',X(:)'*L(:)/2);

part_cell_final = part_cell;
part_final = part;
newX_final = X;
% obj_cluster = obj;

n = size(L,1);
% m = n/k;
count = 0;
k = length(part_cell);
pair_perms = nchoosek(1:k,2);
if k==2;
    mlocal = 1;
else
    mlocal = 3;
end
% 

if isempty(obj);
    obj = L(:)'*X(:)/2;
end

for iter=1:mlocal;
    pairsets = randperm(nchoosek(k,2));
    for pair = pairsets;
        pair_temp = pair_perms(pair,:);
    % for iter=1:mlocal;
%         pair_temp =  randperm(k,2);
% 
%         if iter>1
%             while isempty(setdiff(pair_temp,pair_old))   
%                 pair_temp =  randperm(k,2);
%             end
%         end
        

        part_s_idx = pair_temp(1);
        part_t_idx = pair_temp(2);

        
%  start from the current best results
        part_cell_new = part_cell_final;   
        part_new = part_final;
        Xnew_temp = newX_final;
        
        part_s = reshape(part_cell_new{pair_temp(1)},1,[]);
        part_t = reshape(part_cell_new{pair_temp(2)},1,[]);
        
        
        ns_old = length(part_s);
        nt_old = length(part_t);
        
        
        
        
        nV = length(part_s) + length(part_t);
        subX = Xnew_temp([part_s,part_t],[part_s,part_t]);
        
        subL = L([part_s,part_t],[part_s,part_t]);
        subL(1:(nV+1):end)=0;
        subL(1:(nV+1):end) = -sum(subL);

        suba = a([part_s,part_t]);
        subx = 2*subX(1,:)-1; % transfer to {-1,1}

%       apply 2opt


        [subx_new,new_xlx,old_xlx] = GPKC_2opt(subL,W,suba,subx);
        
        if new_xlx < old_xlx
            subX_simplex = subx_new'*subx_new;
            Xnew_temp([part_s,part_t],[part_s,part_t]) = (subX_simplex +1)/2;
            part_s_new = union(part_s(subx_new(1:ns_old)==subx(1:ns_old)),part_t(subx_new((ns_old+1):end)~=subx((ns_old+1):end)));
            part_t_new = union(part_s(subx_new(1:ns_old)~=subx(1:ns_old)),part_t(subx_new((ns_old+1):end)==subx((ns_old+1):end)));
            part_cell_new{part_s_idx} = part_s_new;
            part_cell_new{part_t_idx} = part_t_new;
            part_new{part_s_idx} = num2str(part_cell_new{part_s_idx}'); 
            part_new{part_t_idx} = num2str(part_cell_new{part_t_idx}');
            obj_current = obj - old_xlx/4 + new_xlx/4;
             %obj = 1/2 xlx
            obj = obj_current;
            part_cell_final =part_cell_new;
            part_final = part_new;
            newX_final = Xnew_temp;
            if iter >1;
                count = 0;
            end
        elseif iter>1;
            count = count + 1;
            if count > nchoosek(k,2)/2
                return
            end
        end
%         pair_old = pair_temp;
    end
end
for i=1:n;
    for j = (i+1):n;
        if newX_final(i,j)~=0 & newX_final(i,j)~=1;
            error('no iterger!!');
            pause;
        end
    end
end

end
%%
function [xnew,cost,oldcost] = GPKC_2opt(L,W,a,x);
% minimization
% obj: xLx
% {-1,1}
% equipartition
% swap 2 vertices in two clusters
% find the best pair from all combinations

m=length(x)/2;

s_idx = find(x==1);
t_idx = find(x==-1);
ns = length(s_idx);
nt = length(t_idx);


Lx = L*x';       % auxiliary vector
d = diag(L);
oldcost = x*Lx;
cost = x*Lx;
% old_Lx = Lx;

delta = x'.*Lx-d;
improve_best = 0;
for s =1:ns;
    i = s_idx(s);
    for t =1:nt;
        j = t_idx(t);

        % check feasibility after swap
        Ws_new = sum(a(s_idx))-a(i) + a(j);
        Wt_new = sum(a(t_idx))-a(j) + a(i);

        if Wt_new > W || Ws_new >W;
            continue;
        else 
            improve_temp = delta(i) + delta(j) - 2*L(i,j)*x(i)*x(j);
            if improve_temp >  improve_best;
                improve_best = improve_temp;
                i_best = i;
                j_best = j;
            end
        end
    end
end


while improve_best > 0.0001;
    
    I1 = find(L(:,i_best));
    I2 = find(L(:,j_best));
   
    if x(i_best)>0;
       Lx(I1) = Lx(I1)  - 2 *L(I1,i_best);
       Lx(I2) = Lx(I2)  + 2 *L(I2,j_best);
    else
       Lx(I1) = Lx(I1)  + 2 *L(I1,i_best);
       Lx(I2) = Lx(I2)  - 2 *L(I2,j_best);
    end;
    
    x(i_best) = -x(i_best);
    x(j_best) = -x(j_best);
    
    
    cost = cost - 4*improve_best;
    delta = x'.*Lx-d;
    
    s_idx = find(x==1);
    t_idx = find(x==-1);
    improve_best = 0;
    
    
    for s =1:ns
        i = s_idx(s);
        for t =1:nt
            j = t_idx(t);
            improve_temp = delta(i) + delta(j) - 2*L(i,j)*x(i)*x(j);

            % check feasibility after swap
            Ws_new = sum(a(s_idx))-a(i) + a(j);
            Wt_new = sum(a(t_idx))-a(j) + a(i);

            if Wt_new > W || Ws_new >W;
                continue;
            else
                if improve_temp >  improve_best;
                    improve_best = improve_temp;
                    i_best = i;
                    j_best = j;
                end
            end
        end
    end
end
    

xnew = x;

if cost~= xnew*L*xnew';
    error('cost error！')
end
end
