function [part, newX,obj,part_cell] = GPKC_rounding_flex_rand(seed,X,C,W,a);
n = size(X,1);
part = {};
part_cell = {};
vertex = 1:n;
% subset_temp =[];
newX = X;

rng(seed);
while ~isempty(vertex);
   Xa = newX*a';
   sum0 = 0;
   logical_vertex = ismember(1:n, vertex);
%    Xa(~logical_vertex) = 0;
%    [~,i] = max(Xa); 

   idx_rnd = randi(length(vertex));
   i = vertex(idx_rnd);
   
   vi = newX(i,:);
   V = newX*vi';
   
   subset_temp = [];
   V(~logical_vertex) = 0;
   [~,idx] =sort(V,'descend');
   idx = intersect(idx,vertex,'stable'); % only keep vertices that are not assigned
   for count= 1:length(idx);
       j = idx(count);
       sum_new = sum0 + a(j);
       if sum_new <= W;
          sum0 = sum_new;
          subset_temp = [subset_temp;j];
       else
           continue;
       end
   end
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
for i =1:length(part_cell);
    if (all(ismember(part_cell{i},parts)==0)) && (sum(a(part_cell{i}) <=W));
         parts = union(parts,part_cell{i},'stable');
    else
        error('infeasible partition');
    end
end

obj = C(:)'*newX(:);
