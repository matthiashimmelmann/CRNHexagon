# CRNHexagon.jl
 
This repository is a supplementary material for the article "Empirically Exploring the Space of Monostationarity in Dual Phosphorylation" by May Cai, Matthias Himmelmann and Birte Ostermann. Through this Julia package, we provide the data used in the article and the code with which the data was generated. To install the package, you first need to clone this repository to your local machine. Then, open the Julia REPL and navigate to the folder "CRNHexagon" that was just created via the command `cd("CRNHexagon")` (this depends on the location of your Julia installation). In the subfolder "src" (reachable with `cd("src")`), where the source code is located. There, use the following line of commands to 

```julia
julia> ]
(@v1.9) pkg> activate cleanCRNenv
(cleanCRNenv) pkg> add HomotopyContinuation, Polyhedra, Plots, Colors, GLMakie, LinearAlgebra
julia> include("CRNHexagon.jl")
```

These commands open the module `CRNHexagon` to your local Julia installation. Afterwards, the package exports four functions. The command

```julia
julia> CRNHexagon.runTest()
```

generates the data from the relative comparison with the cover used in the article "The kinetic space of multistationarity in dual phosphorylation" by Feliu et al. It is referred to as CC(9) in our article. The parameter `boxsize=1` can be used to alter the size of the hypercube `[0,boxsize]^12` from which the points are sampled. `numberOfSamplingRuns=150` determines the amounts of times that 1,000,000 samples are drawn. `prefix="michaelismentontest"` and `suffix="NEW"` allows us to change the name of the data file this method outputs. The data from this test run can be accessed using

```julia
julia> CRNHexagon.computeCoverInvariants()
```

which has the parameters `boxsizes=["0.1","1","10","100"]`, `prefix="michaelismentontest"` and `suffix="NEW"`. This method outputs the LaTex code for Table 1 and Table 2 used in the article associated with this repository (see Section 3). Moreover,

```julia
julia> CRNHexagon.printValues()
```

with the parameters `boxsizes=["0.1","1","10","100"]`, `prefix="michaelismentontest"` and `suffix="NEW"` can be used to print the relative containment information of the covers used to create the Hasse diagram in Figure 6 of the article. Analogously, 

```julia
julia> CRNHexagon.runTest_weighted([ [[1,3,6],[2,5,7],[8,9],[4,10]], [[3,5,8],[1,4,7],[9,10],[2,6]], [[1,5],[7,3],[8,9],[6,2],[4,10]] ])
```

provides a simplicial homotopy between the covers CC(4), CC(9) and CC(15) (see the code in the main file `CRNHexagon` or the article for details). In particular, all equally weighted covers involving these three pure covers are considered by this method. Again, the parameters are `boxsize=1`, `numberOfSamplingRuns=100`, `prefix="linearweight"` and `suffix="4,9"`. However, there is an additional option to set the discretization of the homotopy via `discretization=10`, which sets the step size to `1/discretization`.
