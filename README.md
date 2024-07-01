# CNI-SC
CNI-SC is proposed on a local basis to complete the construction of global causal network.

The "SCI.app" and the "SCI.h" is the C++ code which is more efficient, the usage is to open it in Matlab and add it to the search path, then it can works.

We provide an alarm network as the test data which is "alarm.txt". The users can directly open the main function and run it with one click. More algorithmic process details can be found in the manuscript.

<*contents* >

"Gen_ADJ.M" is the generation method of the candidate adjacent node set of the target node.

"RePC.M" is the network structure learning method based on SCI.

"SCI.h" "SCI.app" "SCI.Mexw64" is the  Conditional independence criterion based on stochastic complexity (SCI).

"TXYCO.M" is the causality orienting algorithm.

"cauculate_lamda.m" is the function cauculate lamda.

"neighbors.m" is the function find neighbors of the target node.

"parents.m"  is the function find parents of the target node.

"setdiff.m"  is the function find the difference of the node sequence.

"subsets.m" is the function find the sub-sets of a specific from a given set.

"union.m" is the function union of two sets of positive integers.

"alarm_500.txt" is the Alarm causal network with sample numbers of 500.

<*Notice* >

Please make sure the MATLAB is the newer version. If the Matlab version or computer hardware configuration is inappropriate, the Matlab software may crash, the user can first step run several times and then press the continue button to complete the program, or reset Matlab and try step run again to solve the issue. 
