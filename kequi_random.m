function [part, newX_final,obj,part_cell] = kequi_random(seed,X,k,C);


n = size(X,1);
m = n/k;

randiter = 1;
% specifiy the random number generator
rng(seed);
obj=inf;
% orthogonal vectors
% r0 = zeros(n,k);
% for i=1:k;
%   r0(i,i)=1;  
% end

% r0 = normrnd(0,1,n,k);



part = {};
part_cell = {};

vertex = 1:n;

% newX = X;



%%

% covert X to (-1/(k-1),1)

X = (k*X-1)/(k-1);

%generate the n veoctors from X
[V] = chol_psd(X);
V = V;

% R=r0;

for iter =1:randiter;
    R = normrnd(0,1,n,k);    
    vertex_finish = [];
    newX = eye(n,n);
    part_cell_temp = cell(k,1);
    for part_temp =1:k;
        ri  = R(:,part_temp);
        score_temp=ri'*V;
        score_temp(vertex_finish) = -100;
        [~,vertex_temp] =maxk(score_temp,m);
        part_cell_temp{part_temp} = vertex_temp;
        vertex_finish = [vertex_finish;vertex_temp];
    end
    for p = 1:k;
        subset_temp = part_cell_temp{p};
        for i = subset_temp;
            for j = 1:n;
                if ismember(j,subset_temp);
                    newX(i,j) = 1;
                else
                    newX(i,j) = 0;
                end
            end
        end
    end
    obj_temp = newX(:)'*C(:);
    if obj_temp < obj;
        obj = obj_temp;
        part_cell = part_cell_temp;
        newX_final = newX;
    end
    % % random rotate matrix
    % R_temp = -ones(n,k);
    % while min(min(R_temp)) <0
        % R_temp=R;
        % for riter=1:k;
        %     rM = eye(n,n);
        %     randgrad = 2*pi*rand;
        %     idxs = randperm(n,2);
        %     idxs = sort(idxs,'ascend');
        %     rM(idxs(1),idxs(1)) = cos(randgrad);
        %     rM(idxs(1),idxs(2)) = -sin(randgrad);
        %     rM(idxs(2),idxs(1)) = sin(randgrad);
        %     rM(idxs(2),idxs(2)) = cos(randgrad);
        %     R_temp = rM*R;
        % end
    % end
    % R = R_temp;
end

for part_id=1:k;
    part{end+1} = num2str(part_cell{part_id});
end

% feasible check
parts = [];
for i =1:k;
    if all(ismember(part_cell{i},parts)==0);
         parts = union(parts,part_cell{i},'stable');
    else
        error('infeasible partition');
    end
end

% for iter=1:randiter;
%     for r=1:n;
%         vi = V(:,r);
%         score = vi'*R;
%         [~,part_idx]= max(score);
%         part_cell{part_idx} = [part_cell{part_idx};r];
%     end
% % % random rotate matrix
%     randgrad = 2*pi*rand;
%     rM(n-2,n-1) = cos(randgrad);
%     rM(n-2,n)= -sin(randgrad);
%     rM(n-1,n-1) = sin(randgrad);
%     rM(n-1,n) = cos(randgrad);
%     
%     R = rM*R;
%     part_finish = [];
%     vertex_finish = [];num2str(subset_temp');
% %     for part_temp=1:k;
% %        if length(part_cell{part_temp})>m;
% %            vertex_idx = part_cell{part_temp};
% %            
% %            score_temp = V*R(:,part_temp);
% %            score(~vertex_idx) = -100;
% %            [~,vertex_new] = maxk(score_temp,m);
% %            part_cell{part_temp} = vertex_new;
% %            part_finish = [part_finish;part_temp];
% %            vertex_finish = union(vertex_finish,vertex_new);
% %        end
% %     end
%     
% end