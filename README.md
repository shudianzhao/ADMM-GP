# ADMM-GP

This package is written based on the following papers:

Angelika Wiegele，Shudian Zhao: SDP-based bounds for graph partition via extended ADMM.


**gerate instances via Rudy in terminal 

source ./rudy_code_randGraphs.sh 


**matlab functions:

aadmm_3b.m : extended ADMM for 3 blocks with adaptive method.

mprw_ineq_general.m : extended ADMM for SDP with inequalities and nonnegative constraints.

post_proc_2.m: the post-processing method for aadmm_3b.m outputs.

post_proc_3.m: the post-processing nethod for mprw_ineq_general.m outputs.

kequi_local.m: 2-opt method for k-equipartition.

kequi_rounding.m: Vertices-clustering for k-equipartition.

kequi_random.m: Hyperplane rounding method for k-equipartition.

GPKC_local.m: 2-opt method for GPKC.

GPKC_rounding_flex_rand.m: Vertices-clustering for GPKC.




