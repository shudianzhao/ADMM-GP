The dataset consists of randomly generated graphs for k-equipartition problems and Graph partition with knapsack constraints (GPKC)


We generate the graphs in following way:

1. Choose edges of a complete graph randomly with probability $20\%$, $50\%$ and $80\%$.
2. The nonzero edge weights are integers between $(0,100]$. 
3. Choose the partition numbers as divisors of the graph


For k-equipartition problem, the instances are:

Keq_rand_p0.80.mat : Adjacency matrices for graphs of size n=100,200,300,400,500,600,700,800,900,1000. The density of graphs is 80%.
Keq_rand_p0.50.mat : Adjacency matrices for graphs of size n=100,200,300,400,500,600,700,800,900,1000. The density of graphs is 50%.
Keq_rand_p0.20.mat : Adjacency matrices for graphs of size n=100,200,300,400,500,600,700,800,900,1000. The density of graphs is 20%.


For GPKC problem, we form the adjacency matrix in the same way and we also randomly generate nonnegative integer weight between (0,1000] for each vertex.

GPKCrandRound_p80.mat:  Adjacency matrices, vetex weights for graphs of size n=100,200,300,400,500. The density of graphs is 80%.
GPKCrandRound_p50.mat:  Adjacency matrices, vetex weights for graphs of size n=100,200,300,400,500. The density of graphs is 50%.
GPKCrandRound_p20.mat:  Adjacency matrices, vetex weights for graphs of size n=100,200,300,400,500. The density of graphs is 20%.

The capacity weights for each graph are in:

GPKCrandRound_p80_W.mat
GPKCrandRound_p50_W.mat
GPKCrandRound_p20_W.mat

%%load the dataset in matlab

load('GPKCrandRound_p80.mat');
load('GPKCrandRound_p80_W.mat.mat');

% graph of size n=100
% adjacency matrix
A =GPKC_data100(1:(end-1),:);
% vertex weight
a = GPKC_data100(end,:);
% capacity weights
W = Wlist100;

