*****************************************************************

	Independent Decomposition of Chemical Reaction Networks

*****************************************************************

This code is used to partition the reactions of a chemcial reaction network into algebraically independent
subnetworks. 

This code was adapted for Julia from the MATLAB function indepDecomp within the package 
COMPILES (see reference below)

Hernandez BS, Lubenia PVN, Johnston MD, Kim JK (2023) 
A framework for deriving analytic steady states of biochemical reaction networks. 
PLOS Computational Biology 19(4): e1011039. https://doi.org/10.1371/journal.pcbi.1011039

The function independent_decomposition requires the following Julia packages to be installed:
RowEchelon, Catalyst, Graphs, Combinatorics, MetaGraphs

INPUT : Julia file containing a Catalyst-defined reaction network (See Catalyst.jl for details on how to define your network).
The name of the file must be "-name-of-your-reaction-network.jl"
To run the independent decomposition, add the name of the reaction network where prompted in lines 150 & 155
OUTPUT : The finest possible partition of the network reactions into independent subnetworks. 
*Note that reactions are numbered in the order they appear in the Catalyst reaction network.

We have included an example network - "Retina.jl", which is a reaction network for cellular cholesterol regulation in the retina. 

