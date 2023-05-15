
## [Results of the Experiments](out/Results/)

## [Precomputed cores and core information](out/Results/CoreInfo/)

## Build Code:

> 1. Clone this repository and navigate into *GCoreApproximation* folder: <br> ```git clone https://github.com/fseiffarth/GCoreApproximation.git && cd GCoreApproximation```
> 2. Create folders *GraphData* and *ExternalLibraries* with the Snap-6.0 library in the parent folder of *GCoreApproximation*: <br> ```mkdir ../GraphData && mkdir ../ExternalLibraries && unzip Snap-6.0.zip -d ../ExternalLibraries/``` <br> (the original version can be found [here](http://snap.stanford.edu/releases/Snap-6.0.zip))
> 3. Create build inside *GCoreApproximation* folder and compile:
   <br> ```mkdir build && cd build && cmake -DCMAKE_BUILD_TYPE=Release .. && make -j 4```
> 4. Download the graphs from http://snap.stanford.edu/data/index.html and save unpacked *.txt* in the *GraphData* folder
> 5. Convert the graphs to the used format: <br> ```./ConvertGraphs```
   

### Run Core Experiment

> 1. Move precomputed cores to the right place: <br> ```cp -a out/Results/CoreInfo/. ../GraphData/```
> 2. Run the executables:
>   1. > ```./ExpSamplingRuntime``` for sampling runtime experiment (output is written to *out/Results/Closure/*)
>   2. > ```./ExpSamplingQuality``` for sampling quality experiment (output is written to *out/Results/Sampling/*)
>   3. > ```./ExpClosureRuntime``` for closure runtime experiment   (output is written to *out/Results/Sampling/*)
>   4. > ```./ExpApproxCore``` for approximate core computation (either you have to compute the cores first (see v.) or use the precomputed cores). To run the experiments from the paper use the below commands (output is written to *out/Results/Approximation/*)
>       > 
>      > *Grid Search*
>      >  - ```./ExpApproxCore -i ../../GraphData/ --generators 5 10 100 500 1000 --threshold 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.15 --threads 2 --outer_loop 0 --small_graphs --no_tree --file_name grid_search_small```
>      >  - ```./ExpApproxCore -i ../../GraphData/ --generators 5 10 100 500 1000 --threshold 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.15 --threads 3 --outer_loop 0 --large_graphs --no_tree --file_name grid_search_large```
>      >
>      > 
>      >*Iterative Approx*
>      >  - ```./ExpApproxCore -i ../../GraphData/ --generators 10 --threshold 0.01 --threads 2 --outer_loop 1 --no_tree --small_graphs --file_name iteration_small```
>      >  - ```./ExpApproxCore -i ../../GraphData/ --generators 10 --threshold 0.01 --threads 3 --outer_loop 1 --no_tree --large_graphs --file_name iteration_large```
>      >
>      > 
>      >*Very large Graphs* (compute approx core of the graphs without knowing the exact core)
>      >  - ```./ExpApproxCore -i ../../GraphData/ --generators 5 10 --threshold 0.01 0.02 --threads 2 --outer_loop 0 --save_load_samples --no_tree --no_core --file_name large_graphs```
>      >
>      >| Optional Arguments | ```-i```  | ```-o```  | ```--threads```  | ```--generators``` | ```--generator_seed``` | ```--threshold``` | ```--core_iterations```  | ```--samples``` | ```--sample_seed```  | ```--max_nodes``` | ```--max_edges``` |
>      >| :---:   | :-: | :-: | :-: | :------------: | :-----------------: | :------------------: | :------------------: | :------------: | :------------: | :------------: | :------------: |
>      >| Seconds | input path | output path | thread num | generator sizes | generator seed | threshold sizes | iterations of the core | number of samples | sample seed | max graph size | max graph edges |
      
>   5. > ```./ExpExactCore``` for exact core computation to recalculate the cores for the graphs use: <br> ```./ExpExactCore -i ../../GraphData/ --threads 8 --recalculate```
>      >
>      >   | Optional Arguments | ```-i```  | ```--threads```  | ```--generators``` | ```--generator_seed``` | ```--core_iterations``` | ```--max_nodes``` | ```--max_edges``` |
>      >   | :---:   | :-: | :-: | :------------: | :-----------------: | :------------------: | :------------: | :------------: |
>      >   | Seconds | input path | thread num | generator number | generator seed | iterations of the core | max graph size | max graph edges |
       

### Run Custom Core Approximation Experiment

>#### With Exact Core Computation:
>>1.  Compute the Core of your graph(s) with <br> ```./ExpExactCore -i ../../GraphData/ --threads 4 --recalculate```
>>2.  Compute the Approximate Cores with <br> ```./ExpApproxCore -i ../../GraphData/ --generators 5 --threshold 0.01 --threads 4 --outer_loop 0 --no_tree --file_name custom_experiment```
>>- Add the argument ```--save_load_samples```, if the samples are too large to fit in the RAM, use ```--no_core``` **only** if no core was computed before
>#### Without Exact Core Computattion:
>>1.  Compute the Approximate Cores with <br> ```./ExpApproxCore -i ../../GraphData/ --generators 5 --threshold 0.01 --threads 4 --outer_loop 0 --no_tree --no_core --file_name custom_experiment```
>>- Add the argument ```--save_load_samples```, if the samples are too large to fit in the RAM, use ```--no_core``` **only** if no core was computed before
>#### Getting the results
>>- The _*.core_ file in _../../GraphData/_ contains the node ids of the core nodes *(relabeled from 0 to n-1)*. 
>>- The _*.core_info_ contains some information about the size, runtime etc.
>>- The _*.approx_core_ file in _../../GraphData/_ contains the node ids of the nodes in the approximated core *(relabeled from 0 to n-1)*.
>>- The _custom_experiment.csv_ and *custom_experiment_summary.csv* files in */out/Approximation/* contain all resp. a summary of the most important results (e.g., runtime, similarity).


### Run Geodesic Closure Operator

>1.  Run the executable to compute the geodesic closure of set of graph vertices:
>    ```./GraphClosure -i path/to/graph.txt -o path/to/output -ids path/to/input_ids.txt (one id per line in the file) [-dist threshold for the geodesic closure based on the distance in the graph (optional)]```


       
