ALGORITHM ENGINEERING SUMMER TERM 2020  
SUBMISSION FOR EXERCISE SHEET 5 BY TEAM 1  
	Johannes Berthold  
	Christopher Denker  
	Luka Stärk  

The program consists of the following files:
* branch_and_bound.py
* vc_graph.py
* reduction.py
* sat_solver.py
* lower_bound.py
* hopcroft_karp.py
* vectorized_func.py
* read_stdin.py
* clique_cover.pyx
* packing_constraints.pyx
* utils.pyx

1. set the stack size limit to 64MB, install dependencies and compiles Cython files by running the script: `. ./config.sh`
(IMPORTANT: execute config inside your bash with exactly `. ./config.sh` and NOT `./config.sh`)

2. install *graph-tool* with a package-manager	
3. Run the program with:  
    `python3 branch_and_bound.py pre=1 branching2=1 red_grp=2 greedy=1`
    
 Dependencies: numpy>=1.18.1, scipy>=1.4.1, graph_tool, python-sat==0.1.5.dev14 and Cython==0.29.20.

The data format for the undirected graph is as specified by the text of exercise sheet 1.
Input is read from standard input and output is written to standard output.

* **branching1** in [0,1], default: 0
    - Splits graph into connected components and solves each for fixed k start with a lower bound until k equals the solution size. 

* **branching2** in [0,1], default: 0
    - Constrained branching with upper bound, two packing constraints and clique cover lower bound updates to cut branches. 

* **greedy** in [0,1], default: 1
  - Apply greedy algorithm. Returns a vertex cover that is not minimal.

* **local_search** in [0,1], default: 1
  - Apply local search on result(s) of greedy algorithm.

* **pre** in [0,1], default: 1
  - Apply preprocessing according to setting of red_grp (reduction group).

* **sat** in [0,1], default: 0
  - Use SAT solver to obtain minimum vertex cover using 'glucose3'.

* **red_grp** in [0,8], default: 2
  - Apply a combination of reduction rules according to index. See `reduction.py`

* **red_frq** in [0,3], default: 0
  - Apply reduction rules (compare red_grp) with varying frequencies according to index. See `reduction.py`.
    