function [part, newX_final,obj,part_cell] = kequi_random(seed,X,k,C);


n = size(X,1);
m = n/k;

randiter = 1;
% specifiy the random number generator
rng(seed);
obj=inf;




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
