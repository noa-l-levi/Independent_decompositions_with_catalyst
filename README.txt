*****************************************************************

	Independent Decomposition of Chemical Reaction Networks

*****************************************************************

This code is used to partition the reactions of a chemcial reaction network into algebraically independent
subnetworks. 

This code was adapted for the Julia language from the function indepDecomp within the package 
COMPILES (see reference below), which uses Matlab

Hernandez BS, Lubenia PVN, Johnston MD, Kim JK (2023) 
A framework for deriving analytic steady states of biochemical reaction networks. 
PLOS Computational Biology 19(4): e1011039. https://doi.org/10.1371/journal.pcbi.1011039

The function independent_decomposition requires the following Julia packages to be installed:
RowEchelon, Catalyst, Graphs, Combinatorics, MetaGraphs

INPUT : Catalyst-defined reaction network
OUTPUT : The finest possible partition of the network reactions into independent subnetworks

Line 148 contains a test chemical reaction network taken from:

Anderson, D.F., Cappelletti, D. 
Discrepancies between extinction events and boundary equilibria in reaction networks. 
J. Math. Biol. 79, 1253–1277 (2019). https://doi.org/10.1007/s00285-019-01394-9

Note that reactions are numbered from in the order they appear in the catalyst defined reaction network.

Lines 161 - 165 contain the skeleton for defining your own reaction network and running the independent
decomposition. See Catalyst.jl for more details on how to define your network.
