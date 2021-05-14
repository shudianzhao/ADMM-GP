

# ADMM-GP

This package is written based on the following papers:

Angelika Wiegele， Shudian Zhao: SDP-based bounds for graph partition via extended ADMM.


## Folders summary

* rand_data: benchmark instances used in numerical experiments;
* extended ADMM: codes for extended ADMM and post processing methods;
* form_model: codes to form SDP problems;
* heuristics: codes to build upper bounds for graph paritioning problems.


## Generate graphs via Rudy in terminal:

source ./rudy_code_randGraphs.sh 


## Matlab functions:

* Form SDP problems

  kequi_form.m : form k-equipartition SDP relaxation
  
  make_GPKC.m : form GPKC SDP relaxation

* Extended ADMM
aadmm_3b.m : extended ADMM for 3 blocks with adaptive method.

mprw_ineq_general.m : extended ADMM for SDP with inequalities and nonnegative constraints.

* Post-processing methods

post_proc_2.m: the post-processing method for aadmm_3b.m outputs.

post_proc_3.m: the post-processing nethod for mprw_ineq_general.m outputs.

* Heuristics

kequi_local.m: 2-opt method for k-equipartition.

kequi_rounding.m: Vector-clustering for k-equipartition.

kequi_random.m: Hyperplane rounding method for k-equipartition.

GPKC_local.m: 2-opt method for GPKC.

GPKC_rounding_flex_rand.m: Vertices-clustering for GPKC.

### Example script:
* k-equipartition
  1. Form the input data (A,b,C) for extended ADMM with the adjacency matrix Ac and the partition number k
  ```
  % matlab code
  % Form the SDP relaxation 
  % <C，X> s.t. A(X) =b, X >= 0, X is postive semidefinite.
  
  [A,b,C] = kequi_form(Ac,k);
  ```
  2. Call the extended ADMM fucntion to solve the SDP relaxation with nonnegative constraints

  ```
  % maximum number of iteration is 100000
  % X1 is the primal solution, (y1,Z1,S1) is the dual solution
  [ X1, y1, Z1,S1] =  aadmm_3b(A, b, C,10000);

  ```
  3. Get the safe lower bound
  ```
  % tuning dual variable y and S to get a feasible solution for the dual SDP problem
  % ynew is the new y variable, LB is the safe lower bound
   [ynew1,LB1] = post_proc_2(Z1, A', C,b, S1);
  ```

  4. Add violated triangle inequalities

  
    * Find violated triangle indices and form the corresponding inequalties <img src="https://render.githubusercontent.com/render/math?math=\mathcal{B}(X) \geq f">
    ```
    % initialize hash,f
    hash =[];
    f = []
    % max_ineq is the maximum number of violated triangle inequalties
    max_ineq = 100;
    [f0, T, hash,gamma,brk]=separation_kequi(X,f,T,hash,max_ineq);
    % form data for SDP with inqualtieis 
     [B,f]=formtri_kc(T,f0,n1);
    L = zeros(n1,n1);
    U = Inf*ones(n1,n1);

    ```
    * Solve the new SDP relaxation with inequalities with extended ADMM  and get the safe lower bound
    ``` 
    %intialize sigma
    sigma= 0.1;
    % tolarence for infeasibility
    tol_err = 1e-5;
    % initalize X, Z with preivious ADMM results X1, Z1

    [ X2, y2,y2_bar, Z2,S2,s2,v2] = mprw_ineq_general(A, B, b, C ,f,L,U,10000,sigma,tol_err,X1,Z1);

    [ynew2,LB2] = post_proc_3(Z2,A',B',C,b,f,S2); 
    ```
   5. Build upper bounds from the SDP solution by heuristics
     * Vector clustering method for the solution of DNN relaxtion
     ```
     % rndseed: random seed
      [part1, newX1,ub1,partcell1] = kequi_rounding(rndseed,X1,k,C);
     ```
     * Hyperplane rounding method
    ```
    % rndseed: random seed
    [part1, newX1,ub1,partcell1] = kequi_random(rndseed,X1,k,C);
    
    ```
     * 2-opt method on output from vector clustering method/ hyperplane rounding method
     ```
       [partcell1_1,part1_1,newX1_1,ub1_1]= kequi_local(L,k,partcell1,part1,newX1,ub1);
     ```
   * GPKC
   
    1. Form the GPKC problem with adjacency matrix Ac, vertex weights a and capacity weight W
        ```
        [A,B,b,f,C]=make_GPKC(Ac,a,W);
        ```
    2. Solve the SDP relaxation with nonnegative constraints
    
      ```
      n1 = size(Ac,1);
      L = zeros(n1,n1);
      U = Inf*ones(n1,n1);
      [ X1 y1, y1_bar, Z1,S1,s1,v1] = mprw_ineq_general(A, B, b, C ,f,L,U,10000)；
      ```
    3. Get the safe lower bound
      
      ```
        [ynew1,LB1] = post_proc_3(Z1,A',B', C,b,f, S1); 
      ```
    4. Build upper bounds from the SDP solution
    
     * vector clustering method
     ```
     [part1, newX1,ub1,partcell1] = GPKC_rounding_flex_rand(rnd_seed,X1,C,W,a);
     ```
    
     * 2-opt method
    ```
    [partcell1_1,part1_1,newX1_1,ub1_1]= GPKC_local(Ac,W0,a,partcell1,part1,newX1,ub1);
    ```
