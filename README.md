# NetworkGenerators
Scripts for generating random phylogenetic networks according to different models.

For help, just run either of the python files with an -h argument.

The HeathNetwork script produces networks by an extention of the Heath model for trees. 
The extention adds hybridization with rates depending on the distance between taxa, or HGT with rates implemented like speciation and extinction in the Heath model.

The BetaSplittingNetwork script is based on the Beta splitting model for trees.
It first produces a tree according to this model, and then adds edges.
Adding the edges can either be done uniformly, or locally.
