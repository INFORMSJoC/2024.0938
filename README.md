[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# A column generation algorithm with dynamic constraint aggregation for minimum sum-of-squares clustering

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

**CG+DCA** is an exact algorithm, based of the branch-and-bound technique, which combines Column Generation (CG)
with Dynamic Constraint Aggregation (DCA) the solving the minimum sum-of-squares clustering problem (MSSC) as described in the paper ["A column generation algorithm with dynamic constraint aggregation for minimum sum-of-squares clustering"](). This repository contains the source code written in C++.

**Important: This code is being developed on an on-going basis at https://github.com/antoniosudoso/mssc-dca. Please go there if you would like to get a more recent version or would like support**

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2024.0938

https://doi.org/10.1287/ijoc.2024.0938.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{SudosoAloise2025,
  author =        {A. M. Sudoso, D. Aloise},
  publisher =     {INFORMS Journal on Computing},
  title =         {A column generation algorithm with dynamic constraint aggregation for minimum sum-of-squares clustering},
  year =          {2025},
  doi =           {10.1287/ijoc.2024.0938.cd},
  url =           {https://github.com/INFORMSJoC/2024.0938},
  note =          {Available for download at https://github.com/INFORMSJoC/2024.0938},
}  
```

## Installation

**CG+DCA** is entirely open-source with the exception of the commercial code [Gurobi](https://www.gurobi.com/) used for solving linear programming problems.

Ubuntu instructions:

1) Install Gurobi (>= 10.0) and compile its C++ libraries

2) Open the makefile `clustering_cpp/Makefile`
    - Set the variable `gurobi_path` with your Gurobi installation folder.

3) Compile the code:

```
cd clustering_cpp/
make
```

This code has been tested under Ubuntu 22.04 LTS with Gurobi 11.0.

## Source File Description

Folder `clustering_cpp`:

- `main_cpp` contains the routine that reads the parameters from the configuration file `config.txt` and executes the column generation algorithm.
- `JobQueue.cpp` and `JobQueue.h` provide the implementation of the queue that stores the nodes of the branch-and-bound tree.
- `Node.h` contains the data structures used to represent the branch-and-bound nodes.
- `ThreadPool.cpp` and `ThreadPool.h` implement a configurable pool of POSIX threads to perform the branch-and-bound search in parallel.
- `cg_util.cpp` and `cg_util.h` provide functions for generating aggregated columns, checking their compatibility, and handling the aggregation partition.
- `Kmeans2D.cpp` and `Kmeans2D.h` contain the implementation of the k-means algorithm for 2D instances.
- `CGStruct.cpp` and `CGStruct.h` define structs for storing data related to column generation and branch-and-bound algorithms.
- `cg_pricing.cpp` and `cg_pricing.h` implement the geometric approach used for the pricing problem.
- `cg_pricing_multi.cpp` and `cg_pricing_multi.h` contain the implementation of the multiple pricing strategy.
- `cg_branch_and_price.cpp` and `cg_branch_and_price.h` implement the branch-and-price algorithm.
- `cg_branching.cpp` and `cg_branching.h` implement the pricing problem in the nodes of the branch-and-bound tree.
- `cg_reduction.cpp` and `cg_reduction.h` implement the dynamic constraint aggregation routine with the aggregated restricted master problem.
- `cg_stabilization.cpp` and `cg_stabilization.h` implement the stabilized (aggregated) restricted master problem.


## Configuration

Various parameters used in **CG+DCA** can be modified in the configuration file `clustering_cpp/config.txt`:

- `KMEANS_INIT` - Initial assignment: use a pre-computed assignment in the `initial_assignment` folder (0); run k-means with several initializations and save the result in the `initial_assignment` folder (>= 1)
- `CG_TIME_LIMIT` - CG time limit (int)
- `CG_MAX_ITER` - Maximum number of CG iterations (int)
- `CG_LOG_STEP` - Display CG summary statistics every CG_LOG_STEP iterations (int)
- `CG_DUAL_DISAGGREGATION` - Strategy for disaggregating the dual variables: average (0); sparse (1); complementary (2)
- `CG_PARTITION_UPDATE` - Update Q with the column having the smallest number of incompatibilities (0); update Q with the minimum reduced cost column (1)
- `CG_DISAGGREGATION_THRESHOLD` - Force the disaggregation of Q when the ratio between the minimum RC compatible and incompatible columns is smaller than the threshold (float)
- `CG_MAX_COMP_COLS` - Maximum number of compatible columns to add in the RMP at each CG iteration (int)
- `CG_MAX_COLS` - Maximum number of columns to generate at each CG iteration (int)
- `CG_MIP_HEURISTIC` - Run MIP-based RMP heuristic to get a feasible clustering at the end of CG algorithm: yes (1); no (0)
- `CG_VERBOSE` - Do not display CG log (0); display CG log (1)
- `CG_INHERIT_PARTITION_CHILD` - Inherit the aggregating partition Q from the parent node: yes (1); no (0)
- `CG_HEURISTIC_PRICING_CHILD` - Run heuristic pricing in children nodes before exact pricing: yes (1); no (0)
- `CG_MAX_COMP_COLS_CHILD` - Maximum number of compatible columns to add in the RMP at each CG iteration in a branching node (int)
- `BB_TOL` - B&B optimality tolerance (float)
- `BB_THREADS` - Single thread (1); multi-thread (> 1 number of threads)
- `BB_MAX_NODES` - Maximum number of B&B nodes to explore (int)
- `BB_VISITING_STRATEGY` - B&B visiting strategy: best first (0); depth first (1); breadth first (2)
- `BB_VERBOSE` - Do not display B&B log (0); display B&B log (1)

## Usage

```
cd clustering_cpp/
./bb <DATASET_NAME> <K> <AGGREGATION_LEVEL> <CONFIG_FILE>
```
- `DATASET_NAME` - dataset name (see the dataset names in the `instances` folder)
- `K` - number of clusters (int)
- `AGGREGATION_LEVEL` - desired aggregation level (initial number of covering constraints) (int)
- `CONFIG_FILE` - name of the configuration file

The dataset file must contain 2D coordinates and include a header line with the number of data points `n` and columns `2`:

```
n 2
x_11 x_12
x_21 x_22
...
...
x_n1 x_n2
```

## Log

If `CG_VERBOSE = 1` the log reports the progress of the CG algorithm at each iteration:

- `Iter` - Iteration number
- `m` - Number of covering constraints (aggregated data points) in agRMP
- `Col ARMP` - Number of columns in the agRMP
- `Opt ARMP` - Optimal objective value of agRMP
- `Time ARMP` - Computational time in seconds for solving agRMP
- `Comp` - Number of compatible columns added to the agRMP
- `Best Inc` - Minimum reduced cost column value
- `Best Comp` - Minimum reduced cost compatible column value
- `Ratio` - Ratio between `Best Comp` and `Best Inc`
- `Time PP` - Computational time in seconds for solving the pricing problem
- `LB Inc` - Valid lower bound on the optimal MSSC value
- `CG Gap [%]` - Percentage gap between `Opt ARMP` and `LB Inc`
- `IP Gap [%]` - Percentage integrality gap
- `Time [s]` - Cumulative computational time in seconds

If `BB_VERBOSE = 1` the log reports the progress of the B&B algorithm at each node:

- `ID` - ID of the node
- `TYPE` - Node type (RT - root node; ML - must-link node; CL - cannot-link node)
- `I J` - Indices of branching decision
- `LB` - Lower bound of the node
- `ITER` - Number of CG iterations of the node
- `N START` - Initial number of covering constraints (aggregated data points) at the node
- `N END` - Final number of covering constraints (aggregated data points) at the node
- `TIME [s]` - Running time in seconds of the current node
- `UB` - Upper bound of the current node
- `GUB` - Global upper bound
- `GAP [%]` - Optimality gap at the current node
- `GGAP [%]` - Global optimality gap
- `QUEUE` - Number of open nodes in the queue

## Example

```
cd clustering_cpp/
./bb ruspini 8 8 config
```

CG execution at the root node (`BB_MAX_NODES = 1`) with summary statistics and metrics

```

K-MEANS INIT: 1000

DATA: ruspini
K: 8
N_COMPONENTS: 8

CG_TIME_LIMIT: 86400
CG_MAX_ITER: 1000000
CG_LOG_STEP: 10
CG_DUAL_DISAGGREGATION: 0
CG_PARTITION_UPDATE: 1
CG_DISAGGREGATION_THRESHOLD: 0.25
CG_MAX_COMP_COLS: 10
CG_MAX_COLS: 100
CG_MIP_HEURISTIC: 0
CG_VERBOSE: 1
CG_INHERIT_PARTITION_CHILD: 1
CG_HEURISTIC_PRICING_CHILD: 1
CG_MAX_COMP_COLS_CHILD: 10

BB_TOL: 0.0001
BB_THREADS: 1
BB_MAX_NODES: 1
BB_VISITING_STRATEGY: 0
BB_VERBOSE: 1

python run_kmeans.py ./instances/ruspini.txt 8 1000 ./initial_assignment/ruspini_8.txt
Input data shape: (75, 2)
Using kmeans (python) solution with value: 6149.64

    ------------------------------------------------------------------------  COLUMN GENERATION WITH DCA FOR THE MSSC PROBLEM  --------------------------------------------------------------------
    Iter          m        Col ARMP        Opt ARMP    Time ARMP     Comp        Best Inc       Best Comp       Ratio      Time PP          LB Inc      CG Gap [%]     IP Gap [%]       Time [s]
       0          8               8         6149.64         0.00        0         -539.82            0.00        0.00         0.00         1831.10          70.224          70.224            0.00
      10         24              29         6149.64         0.00        6         -234.85          -61.05        0.26         0.00         4270.82          24.167          24.167            0.01
      20         36              71         6146.52         0.00       10          -60.59          -31.20        0.52         0.00         5661.82           7.886           7.933            0.01
      30         43             121         6146.52         0.00        3          -28.84           -1.16        0.04         0.00         5915.80           3.754           3.803            0.02
      40         50             163         6146.52         0.00        5          -81.31          -81.31        1.00         0.00         5496.06           2.549           2.599            0.02
      50         55             203         6146.52         0.00        3           -7.63           -3.84        0.50         0.00         6085.50           0.993           1.043            0.03
      55         57             211         6146.52         0.00        0            0.00            0.00        0.00         0.00         6146.52           0.000           0.051            0.03
      Active LB: 1      Active UB: 2
      57         57             214         6148.59         0.00        0            0.00            0.00        0.00         0.00         6148.59           0.000           0.017            0.03
      Stop: problem solved!
      Optimal solution is not integral
    -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


Number of points: 75
Number of clusters: 8
Number of components: 8

Iter: 58
Initial upper bound: 6149.64
Final lower bound: 6148.59
Final upper bound: 6149.64
Initial IP gap: 0.02 [%]
Final IP gap: 0.02 [%]

Initial number of covering constraints: 8
Final number of covering constraints: 57
Number of times the partition has been updated: 31
Count active bounds on lambda: 1

Total time master (sec) = 0.01
Total time pricing (sec) = 0.02
Total time (sec) = 0.03

Average time master (sec): 0.00
Average time pricing (sec): 0.00
Average incompatibilities (best d-incompatible column) = 1.58
Average incompatibilities (best incompatible column) = 1.58
Average covering = 39.10
Average number of generated columns = 54.38
Average number of compatible columns = 4.33

```

B&B execution with summary statistics and metrics

```

K-MEANS INIT: 1000

DATA: ruspini
K: 8
N_COMPONENTS: 8

CG_TIME_LIMIT: 86400
CG_MAX_ITER: 1000000
CG_LOG_STEP: 100
CG_DUAL_DISAGGREGATION: 0
CG_PARTITION_UPDATE: 1
CG_DISAGGREGATION_THRESHOLD: 0.25
CG_MAX_COMP_COLS: 10
CG_MAX_COLS: 100
CG_MIP_HEURISTIC: 0
CG_VERBOSE: 0
CG_INHERIT_PARTITION_CHILD: 1
CG_HEURISTIC_PRICING_CHILD: 1
CG_MAX_COMP_COLS_CHILD: 10

BB_TOL: 0.0001
BB_THREADS: 1
BB_MAX_NODES: 100
BB_VISITING_STRATEGY: 0
BB_VERBOSE: 1

python run_kmeans.py ./instances/ruspini.txt 8 1000 ./initial_assignment/ruspini_8.txt
Input data shape: (75, 2)
Using kmeans (python) solution with value: 6149.64

      ------------------------------------------------------  BRANCH AND BOUND WITH COLUMN GENERATION AND DCA FOR THE MSSC PROBLEM  --------------------------------------------------
      ID     TYPE        I        J              LB       ITER       COLS    N START      N END       TIME              UB             GUB         GAP [%]        GGAP [%]     QUEUE
       0       RT       -1       -1         6148.59         58        214          8         57       0.03         6149.64         6149.64           0.017           0.017          0
      Running exact pricing...
      Exact pricing reduced cost: -0.000
       1       ML        3        9        6149.639         30        358         57         65       0.66         6149.64         6149.64          -0.000           0.017          1
      Running exact pricing...
      Exact pricing reduced cost: -0.000
       2       CL        3        9        6156.207         34        425         57         65       0.35         6149.64         6149.64          -0.107           0.017          0
      --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Final MSSC solution: 6149.639
Gap [%]: 0.000
Nodes: 3
Time [s]: 1.042

```

The assignment vector corresponding to the optimal MSSC solution is saved in the `optimal_assignment` folder.

