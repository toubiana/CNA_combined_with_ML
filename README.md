# Source code: CNA_combined_with_ML
Code to compute network properties (features) for subgraphs of a correlation-based network
## DESCRIPTION:
The file CNA_combind_with_ML.R found under URL: contains all functions used to compute network features as described in the
corresponding publication by Toubiana et al. (Communications Biology, 2019).
For each pathway/subgraph all feature functions need to be executed in order
to obtain feature values for all metabolic pathways.
A container function 'network_features' is provided, which allows all feature functions
to be triggered at once.
The function receives the following parameters:

      graph: which represents the entire network (igraph object)
      
      subgraph: which represents the subgraph corresponding to a pathway
                or a set of random nodes (igraph object)
                
      subgraph.node.index = provides the node indices of the subgraph
                           by default 1
                           
      subgraph.edge.index = provides the link indices of the subgraph
                           by default 1
                           
      community.detection = a logical parameter. Shall community detection
                           be performed? As this is a time consuming process
                           it is adviced to run this process only once.
                           The communities are coded as global variables

 Feature functions where the mean, standard deviation and two central moments are computed are
 structured as the following:
 
 One HEAD FUNCTION, which returns the mean, sd, skewness, and kurtosis as one variable
 seperate functions, which return on of central moments only, e.g. see functions 2-5

 ## PREREQUISITES:
 The following code is a correlation based network (CN), where nodes in the network
 represent metabolites and links between them the significant correlations.
 Prior to constructing the network, spurious links ought to be removed.
 For instructions on how to construct a correlation-based network (CN), we refer the user to
 Toubiana et al, 2013, Network analysis: tackling complex data to study plant metabolism,
 Trends in Biotechnology, Vol. 31, Issue 1, Pages 29-36.
 The CN should be converted into an igraph object.
 Also, pathways/subgraphs should be supplied as igraph objects.
 NOTE: a CN is not supplied here.
 The libraries 'igraph' and 'e1071' need to be preinstalled in order for the corresponding
 source code to run.
