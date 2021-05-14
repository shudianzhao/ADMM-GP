function [part_final, newX_final,obj,part_cell_final] = kequi_rounding(seed,X,k,C);
n = size(X,1);
m = n/k;
obj=inf;

%% vector clustering method
%generate the n veoctors from X
 

miter = 1; %100
% specifiy the random number generator
rng(seed);
%%
for round =1:miter;
    part = {};
    part_cell = {};
    vertex = 1:n;
    
    newX = X;
    for count =1:k;
    
        idx_rnd = randi(length(vertex));
        i = vertex(idx_rnd);
        vi = newX(i,:);
    
        VV=newX*vi';
        logical_vertex = ismember(1:n, vertex);
        VV(~logical_vertex)=-100;
        [~,subset_temp]=  maxk(VV,m);
        part{end+1} = num2str(subset_temp');
        part_cell{end+1} = subset_temp;
        for i = subset_temp;
            for j = 1:n;
                if ismember([j],subset_temp);
    %                 X(i,j)  = 1;
                    newX(i,j) = 1;
                else
    %                 X(i,j) = 0;
                    newX(i,j) = 0;
                end
            end
        end
        vertex = setdiff(vertex,subset_temp);
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
    obj_current = C(:)'*newX(:);

    if obj_current < obj;
        obj = obj_current;
        part_cell_final =part_cell;
        part_final = part;
        newX_final = newX;
    end

end
