# Social Ranking Manipulability
## Tahar Allouche, Bruno Escoffier, Meltem Ozturk, Stefano Moretti

# IJCAI-2020

**Paper Abstract ([paper available here](https://www.ijcai.org/proceedings/2020/0003.pdf)):**\
We investigate the issue of manipulability for social ranking rules, where the goal is to rank individuals given the ranking of coalitions formed by them and each individual prefers to reach the highest positions in the social ranking. This problem lies at the intersection of computational social choice and the algorithmic theory of power indices. Different social ranking rules have been recently proposed and studied from an axiomatic point of view. In this paper, we focus on rules representing three classical approaches in social choice theory: the marginal contribution approach, the lexicographic approach and the (ceteris paribus) majority one. We first consider some particular members of these families analysing their resistance to a malicious behaviour of individuals. Then, we analyze the computational complexity of manipulation, and complete our theoretical results with simulations in order to analyse the manipulation frequencies and to assess the effects of manipulations.


**Repository:**\
This repository contains the julia code used in the simulations. It uses JuMP with the [CBC optimizer](https://github.com/jump-dev/Cbc.jl) for linear optimization.

A related R package for usual social ranking operations can be found on [CRAN](https://cran.r-project.org/web/packages/socialranking/index.html).


**Structure:**\
- source folder: contains the following scripts:
    - `generation.jl`: For generating random Power Relations.
	- `helpers.jl`: Contains some useful functions to process power relations.
	- `banzhaf.jl`: Contains the linear optimization program to manipulate Banzhaf social rankings.
	- `copeland.jl`: Contains the linear optimization program to manipulate Copeland social rankings.
	- `kramer_simpson.jl`: Contains the linear optimization program to manipulate Kramer-Simpson social rankings.
	- `experiments.jl`: This is where tests run through whole datasets to get some statistics.


**Code:**\
Here's how to run simulations for each rule:
1. Run julia on terminal: `$ julia`
2. Include the `experiments.jl` script: `include("experiments.jl")`.
3. Execute the function: `SolveDataSet()`, passing the number of agents, the number of samples, and choose the social ranking method.
5. The final results will be displayed at the end of the compilation.

