## ------------------------------DISCLAIMER-------------------------------------##
## This SOFTWARE PRODUCT is provided by THE PROVIDER "as is" and "with all faults."
## THE PROVIDER makes no representations or warranties of any kind concerning the safety,
## suitability, lack of viruses, inaccuracies, typographical errors, or other harmful
## components of this SOFTWARE PRODUCT. There are inherent dangers in the use of any software,
## and you are solely responsible for determining whether this SOFTWARE PRODUCT is compatible
## with your equipment and other software installed on your equipment.
## You are also solely responsible for the protection of your equipment and backup of your data,
## and THE PROVIDER will not be liable for any damages you may suffer in connection with using,
## modifying, or distributing this SOFTWARE PRODUCT.


## ------------------------------AUTHOR: DAVID TOUBIANA----------------------------------------##

## ---------------------------------DESCRIPTION------------------------------------------------##
## this file contains all functions used to compute network features as described in the
## corresponding publication by Toubiana et al. (Communications Biology, 2019)
## for each pathway/subgraph all feature functions need to be executed in order
## to obtain feature values for all metabolic pathways
## A container function 'network_features' is provided, which allows all feature functions
## to be triggered at once
## the function receives the following parameters:
##      graph: which represents the entire network (igraph object)
##      subgraph: which represents the subgraph corresponding to a pathway
##                or a set of random nodes (igraph object)
##      subgraph.node.index: provides the node indices of the subgraph
##                           by default 1
##      subgraph.edge.index: provides the link indices of the subgraph
##                           by default 1
##      community.detection: logical parameter. Shall community detection
##                           be performed? As this is a time consuming process
##                           it is adviced to run this process only once.
##                           the communities are coded as global variables
##
## feature functions where the mean, sd and two central moments are computed are
## structured as the following:
## One HEAD FUNCTION, which returns the mean, sd, skewness, and kurtosis as one variable
## seperate functions, which return on of central moments only, e.g. see functions 2-5
## ---------------------------------DESCRIPTION END--------------------------------------------##

## ------------------------------PREREQUISITES-------------------------------------------------##
## the following code is based on a correlation based network, where nodes in the network
## represent metabolites and links between them the significant correlations
## prior to constructing the network, spurious links ought to be removed
## for instructions on how to construct a correlation-based network (CN), we refer the user to
## Toubiana et al, 2013, Network analysis: tackling complex data to study plant metabolism,
## Trends in Biotechnology, Vol. 31, Issue 1, Pages 29-36
## the CN should be converted into an igraph object
## also, pathways/subgraphs should be supplied as igraph objects
## NOTE: a CN is not supplied here
## the libraries 'igraph' and 'e1071' need to be preinstalled in order for the following
## source code to run
## ------------------------------PREREQUISITES END-------------------------------------------##

## --------------------------------LICENSE---------------------------------------------------##

# License
# 
# The following license governs the use of the source code for CNA_combined_with_ML in academic and educational environments. 
# Commercial use requires a commercial license from the author directly: 
#     David Toubiana
# email: dtoubiana@ucdavis.edu or david.toubiana@gmail.com
# 
# ACADEMIC PUBLIC LICENSE
# 
# Original license written by Andras Varga and modified by David Toubiana (license text is in public domain)
# 
# Preamble
# 
# This license contains the terms and conditions of using the source code for CNA_combined_with_ML in noncommercial settings: at academic
# institutions for teaching and research use, and at non-profit research organizations. You will find that this license provides 
# noncommercial users of the source code for CNA_combined_with_M with rights that are similar to the well-known GNU General Public License, 
# yet it retains the possibility for the source code for CNA_combined_with_M authors to financially support the development by selling 
#commercial licenses. In fact, if you intend to use the source code for CNA_combined_with_M in a “for-profit” environment, where research
# is conducted to develop or enhance a product, is used in a commercial service offering, or when a commercial company uses the source 
# code for CNA_combined_with_M to participate in a research project (for example government-funded or EU-funded research projects), 
# then you need to obtain a commercial license for the source code for CNA_combined_with_M. In that case, please contact the Author to
# inquire about commercial licenses.
# 
# What are the rights given to noncommercial users? Similarly to GPL, you have the right to use the software, to distribute copies, to 
# receive source code, to change the software and distribute your modifications or the modified software. Also similarly to the GPL, 
# if you distribute verbatim or modified copies of this software, they must be distributed under this license.
# 
# By modeling the GPL, this license guarantees that you’re safe when using the source code for CNA_combined_with_M in your work, for 
# teaching or research. This license guarantees that the source code for CNA_combined_with_M will remain available free of charge for 
# nonprofit use. You can modify the source code for CNA_combined_with_M to your purposes, and you can also share your modifications. 
# Even in the unlikely case of the authors abandoning the source code for CNA_combined_with_M entirely, this license permits anyone to 
# continue developing it from the last release, and to create further releases under this license.
# 
# We believe that the combination of noncommercial open-source and commercial licensing will be beneficial for the whole user community, 
# because income from commercial licenses will enable faster development and a higher level of software quality, while further enjoying 
# the informal, open communication and collaboration channels of open source development.
# 
# The precise terms and conditions for using, copying, distribution and modification follow.
# 
# TERMS AND CONDITIONS FOR USE, COPYING, DISTRIBUTION AND MODIFICATION
# Definitions
# 
# “Program” means a copy of the source code for CNA_combined_with_M, which is said to be distributed under this Academic Public License.
# 
# “Work based on the Program” means either the Program or any derivative work under copyright law: that is to say, a work containing the 
# Program or a portion of it, either verbatim or with modifications and/or translated into another language. (Hereinafter, translation is 
# included without limitation in the term “modification”.)
# 
# “Using the Program” means any act of creating executables that contain or directly use libraries that are part of the Program, running 
# any of the tools that are part of the Program, or creating works based on the Program.
# 
# Each licensee is addressed as “you”.
# 
# §1. Permission is hereby granted to use the Program free of charge for any noncommercial purpose, including teaching and research at 
# universities, colleges and other educational institutions, research at non-profit research institutions, and personal non-profit 
# purposes. For using the Program for commercial purposes, including but not restricted to consulting activities, design of commercial
# hardware or software networking products, and a commercial entity participating in research projects, you have to contact the Author 
# for an appropriate license. Permission is also granted to use the Program for a reasonably limited period of time for the purpose of
# evaluating its usefulness for a particular purpose.
# 
# §2. You may copy and distribute verbatim copies of the Program’s source code as you receive it, in any medium, provided that you 
# conspicuously and appropriately publish on each copy an appropriate copyright notice and disclaimer of warranty; keep intact all the
# notices that refer to this License and to the absence of any warranty; and give any other recipients of the Program a copy of this 
# License along with the Program.
# 
# §3. You may modify your copy or copies of the Program or any portion of it, thus forming a work based on the Program, and copy and 
# distribute such modifications or work under the terms of Section 2 above, provided that you also meet all of these conditions:
#     
# a) You must cause the modified files to carry prominent notices stating that you changed the files and the date of any change.
# 
# b) You must cause any work that you distribute or publish, that in whole or in part contains or is derived from the Program or any part 
#    thereof, to be licensed as a whole at no charge to all third parties under the terms of this License.
# 
# These requirements apply to the modified work as a whole. If identifiable sections of that work are not derived from the Program, and
# can be reasonably considered independent and separate works in themselves, then this License, and its terms, do not apply to those 
# sections when you distribute them as separate works. But when you distribute the same sections as part of a whole which is a work based 
# on the Program, the distribution of the whole must be on the terms of this License, whose regulations for other licensees extend to the 
# entire whole, and thus to each and every part regardless of who wrote it. (If the same, independent sections are distributed as part of 
# a package that is otherwise reliant on, or is based on the Program, then the distribution of the whole package, including but not 
# restricted to the independent section, must be on the unmodified terms of this License, regadless of who the author of the included 
# sections was.)
# 
# Thus, it is not the intent of this section to claim rights or contest your rights to work written entirely by you; rather, the intent
# is to exercise the right to control the distribution of derivative or collective works based or reliant on the Program.
# 
# In addition, mere aggregation of another work not based on the Program with the Program (or with a work based on the Program) on a 
# volume of storage or distribution medium does not bring the other work under the scope of this License.
# 
# §4. You may copy and distribute the Program (or a work based on it, under Section 3) in object code or executable form under the terms
#     of Sections 2 and 3 above provided that you also do one of the following:
#     
# a) Accompany it with the complete corresponding machine-readable source code, which must be distributed under the terms of Sections 2
#    and 3 above on a medium customarily used for software interchange; or,
# 
# b) Accompany it with a written offer, valid for at least three years, to give any third party, for a charge no more than your cost of
#    physically performing source distribution, a complete machine-readable copy of the corresponding source code, to be distributed 
#    under the terms of Sections 2 and 3 above on a medium customarily used for software interchange; or,
# 
# c) Accompany it with the information you received as to the offer to distribute corresponding source code. (This alternative is allowed 
#    only for noncommercial distribution and only if you received the program in object code or executable form with such an offer, in
#    accord with Subsection b) above.)
# 
# The source code for a work means the preferred form of the work for making modifications to it. For an executable work, complete source
# code means all the source code for all modules it contains, plus any associated interface definition files, plus the scripts used to
# control compilation and installation of the executable. However, as a special exception, the source code distributed need not include 
# anything that is normally distributed (in either source or binary form) with the major components (compiler, kernel, and so on) of the
# operating system on which the executable runs, unless that component itself accompanies the executable.
# 
# If distribution of executable or object code is made by offering access to copy from a designated place, then offering equivalent access
# to copy the source code from the same place counts as distribution of the source code, even though third parties are not compelled to 
# copy the source along with the object code.
# 
# §5. You may not copy, modify, sublicense, or distribute the Program except as expressly provided under this License. Any attempt 
# otherwise to copy, modify, sublicense or distribute the Program is void, and will automatically terminate your rights under this License. 
# However, parties who have received copies, or rights, from you under this License will not have their licenses terminated so long as 
# such parties remain in full compliance.
# 
# §6. You are not required to accept this License, since you have not signed it. Nothing else grants you permission to modify or 
# distribute the Program or its derivative works; law prohibits these actions if you do not accept this License. Therefore, by modifying
# or distributing the Program (or any work based on the Program), you indicate your acceptance of this License and all its terms and 
# conditions for copying, distributing or modifying the Program or works based on it, to do so.
# 
# §7. Each time you redistribute the Program (or any work based on the Program), the recipient automatically receives a license from the
# original licensor to copy, distribute or modify the Program subject to these terms and conditions. You may not impose any further 
# restrictions on the recipients’ exercise of the rights granted herein. You are not responsible for enforcing compliance by third parties
# to this License.
# 
# §8. If, as a consequence of a court judgment or allegation of patent infringement or for any other reason (not limited to patent 
# issues), conditions are imposed on you (whether by court order, agreement or otherwise) that contradict the conditions of this License,
# they do not excuse you from the conditions of this License. If you cannot distribute so as to satisfy simultaneously your obligations 
# under this License and any other pertinent obligations, then as a consequence you may not distribute the Program at all. For example, 
# if a patent license would not permit royalty-free redistribution of the Program by all those who receive copies directly or indirectly 
# through you, then the only way you could satisfy both it and this License would be to refrain entirely from distribution of the Program.
# 
# If any portion of this section is held invalid or unenforceable under any particular circumstance, the balance of the section is 
# intended to apply and the section as a whole is intended to apply in other circumstances.
# 
# §9. If the distribution and/or use of the Program are restricted in certain countries either by patents or by copyrighted interfaces,
# the original copyright holder who places the Program under this License may add an explicit geographical distribution limitation 
# excluding those countries, so that distribution is permitted only in or among countries not thus excluded. In such case, this License
# incorporates the limitation as if written in the body of this License.
# 
# NO WARRANTY
# 
# §10. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. 
# EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM “AS IS” WITHOUT WARRANTY OF ANY 
# KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
# PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME 
# THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.
# 
# §11. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED ON IN WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY 
# AND/OR REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR 
# CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
# RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN 
# IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
# 
# END OF TERMS AND CONDITIONS
# 
# In case the above text differs from the license file in the source distribution, the latter is the valid one.

## --------------------------------LICENSE END-----------------------------------------------##

## ------------------------------FUNCTIONS---------------------------------------------------##
## function network_features
network_features  <- function(graph, subgraph, subgraph.node.index=1,subgraph.edge.index=1,
                              community.detection=F){
    # this is the function that will be called, which will make calls to all other functions
    # using parameters graph and subgraph to compute features
    # libraries to be imported: e1071, igraph
    # library e1071 has to be imported for skewness and kurtosis measures to run
    library(e1071); library(igraph)
    
    # variable to be returned
    #------------------------------------------------------------------------------------------------
    # precalculating community structures if parameter community.detection=T - by default it is FALSE
    # only want to do it once for the graph
    if(community.detection){
        
        ## <<-  defining global variables
        gr <<- graph
        E(gr)$weight <<- abs(E(gr)$weight)
        # fast greedy, non-weighted and weighted
        fgc <<- fastgreedy.community(gr,weights=NULL)
        fgc.w <<- fastgreedy.community(gr,weights=E(gr)$weight)
        # walktrap, non-weighted and weighted
        wtc <<- walktrap.community(gr,weights=NULL)
        wtc.w <<- walktrap.community(gr,weights=E(gr)$weight)
        # edge betweenness community
        # doesn't work for some peculiar reason
        #ebc <<- edge.betweenness.community(gr,weights=NULL)
        #ebc.w <<- edge.betweenness.community(gr,weights=E(gr)$weight)
        # label propagation community
        lpc <<- label.propagation.community(gr,weights=NULL)
        lpc.w <<- label.propagation.community(gr,weights=E(gr)$weight)
        # leading eigenvector community
        lec <<- leading.eigenvector.community(gr,weights=NULL)
        lec.w <<- leading.eigenvector.community(gr,weights=E(gr)$weight)
        # multi level community
        mlc <<- multilevel.community(gr,weights=NULL)
        mlc.w <<- multilevel.community(gr,weights=E(gr)$weight)
    }
    #------------------------------------------------------------------------------------------------
    return.variable = NULL
    
    return.variable <- cbind(return.variable,subgraph.num.edges(subgraph))                    # no. 1   ok
    return.variable <- cbind(return.variable,t(subgraph.degree(subgraph)))                    # no. 2-5   ok
    return.variable <- cbind(return.variable,t(subgraph.weighted.degree(subgraph)))           # no. 6-10   ok
    return.variable <- cbind(return.variable,t(subgraph.abs.weighted.degree(subgraph)))       # no. 11-15   ok
    print(15)
    return.variable <- cbind(return.variable,geodesic.distance.subgraph(subgraph))            # no. 16 ok
    return.variable <- cbind(return.variable,t(geodesic.distance.subgraph.moments(subgraph)))    # no 17-20 
    return.variable <- cbind(return.variable,t(geodesic.distance.weighted.subgraph.moments(subgraph)))    # no 21-25
    #return.variable <- cbind(return.variable,group.betweenness.subgraph(graph,subgraph.node.index)) # no 26 - uses lists - slower
    # only for testing purposes - uncomment when done
    # is too slow put NAs instead
    return.variable <- cbind(return.variable,NA)
    return.variable <- cbind(return.variable,NA)
    # return.variable <- cbind(return.variable,group.betweenness.subgraph(graph,subgraph.node.index)) # no 26 - uses loops - faster
    # return.variable <- cbind(return.variable,weighted.group.betweenness.subgraph(graph,subgraph.node.index)) # no 27 - uses lists - faster
    # return.variable <- cbind(return.variable,weighted.group.betweenness.subgraph.2(graph,subgraph.node.index)) # no 27 - uses loops - slower
    return.variable <- cbind(return.variable,t(closeness.centrality.subgraph(subgraph)))    # no 28-31
    (print(31))
    return.variable <- cbind(return.variable,t(weighted.closeness.centrality.subgraph(subgraph)))    # no 32-35
    return.variable <- cbind(return.variable,t(closeness.centrality.graph.subgraph(graph,subgraph.node.index)))   # no 36-39  
    return.variable <- cbind(return.variable,t(weighted.closeness.centrality.graph.subgraph(graph,subgraph.node.index)))  # no 40-43
    return.variable <- cbind(return.variable,diameter.subgraph(subgraph))  # no 44
    return.variable <- cbind(return.variable,weighted.diameter.subgraph(subgraph)) # no 45
    print(45)
    # only for testing purposes - uncomment when done
    # is too slow put NAs instead
    return.variable <- cbind(return.variable,NA)
    return.variable <- cbind(return.variable,NA)
    # return.variable <- cbind(return.variable,diameter.centrality(graph,subgraph.node.index))  # no 46
    # return.variable <- cbind(return.variable,weighted.diameter.centrality(graph,subgraph.node.index)) # no 47
    return.variable <- cbind(return.variable,t(clustering.coefficient.subgraph(subgraph)))  # no 48-52
    return.variable <- cbind(return.variable,t(weighted.clustering.coefficient.subgraph(subgraph)))  # no 53-57
    return.variable <- cbind(return.variable,
                             t(clustering.coefficient.subgraph.graph(graph,subgraph.node.index))) # no 58-61
    print(61)
    return.variable <- cbind(return.variable,
                             t(weighted.clustering.coefficient.subgraph.graph(graph,subgraph.node.index))) # no 62-65
    return.variable <- cbind(return.variable,t(node.betweenness.subgraph(graph,subgraph.node.index))) # no 66-69
    return.variable <- cbind(return.variable,t(weighted.node.betweenness.subgraph(graph,subgraph.node.index))) # no 70-73
    return.variable <- cbind(return.variable,t(edge.betweenness.subgraph(graph,subgraph.edge.index))) # no 74-77
    print(77)
    return.variable <- cbind(return.variable,t(weighted.edge.betweenness.subgraph(graph,subgraph.edge.index))) # no 78-81
    return.variable <- cbind(return.variable,assortativity.subgraph(subgraph)) # no 82
    return.variable <- cbind(return.variable,density.subgraph(subgraph)) # no 83
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,fgc)) # no 84
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,fgc.w)) # no 85
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,wtc)) # no 86
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,wtc.w)) # no 87
    # edge betweeness community detection doesn't work for some reason, putting NAs instead
    return.variable <- cbind(return.variable,NA)
    return.variable <- cbind(return.variable,NA)
    #return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,ebc)) # no 88
    #return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,ebc.w)) # no 89
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,lpc)) # no 90
    print(90)
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,lpc.w)) # no 91
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,lec)) # no 92
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,lec.w)) # no 93
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,mlc)) # no 94
    return.variable <- cbind(return.variable,greatest.community.subgraph(graph,subgraph.node.index,mlc.w)) # no 95
    # friends measures: common friends, total friend, distinct friends, mixed friends, and
    # preferential attachment score
    return.variable <- cbind(return.variable,t(friend.measures(graph,subgraph.node.index))) # no 96 - 100
    print(100)
    #return(return.variable)
    # the jaccard coefficient is the commonfriends measure / totalfriends measure
    # features 101-105
    if(length(subgraph.node.index)<2){
        sni <- unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))
        return.variable <- cbind(return.variable,t(jaccard.coefficient(graph,sni,
                                                                       return.variable[length(return.variable)-4],
                                                                       return.variable[length(return.variable)-3])))
    } else{
        return.variable <- cbind(return.variable,t(jaccard.coefficient(graph,subgraph.node.index,
                                                                       return.variable[length(return.variable)-4],
                                                                       return.variable[length(return.variable)-3])))
    }
    print(105)
    return.variable <- cbind(return.variable,t(friends.measure.subgraph(subgraph)))  # no 106-110
    return.variable <- cbind(return.variable,t(friends.measure.graph(graph,subgraph.node.index)))  # no 111-115
    print(115)
    # whole network features
    return.variable <- cbind(return.variable,num.edges.graph(graph,subgraph.node.index)) # no 116
    #print(116)
    return.variable <- cbind(return.variable,t(degree.graph(graph,subgraph.node.index))) # no 117-120
    #print(120)
    return.variable <- cbind(return.variable,t(weighted.degree.graph(graph,subgraph.node.index))) # no 121-125
    print(125)
    return.variable <- cbind(return.variable,t(abs.weighted.degree.graph(graph,subgraph.node.index))) # no 126-130
    #print(130)
    return.variable <- cbind(return.variable,t(geodesic.distances.graph(graph,subgraph.node.index))) # no 131-135
    return.variable <- cbind(return.variable,t(geodesic.distances.weighted.graph(graph,
                                                                                 subgraph.node.index))) # no 136-140
    # is too slow put NAs instead
    return.variable <- cbind(return.variable,t(rep(NA,8)))
    #return.variable <- cbind(return.variable,NA)
    # return.variable <- cbind(return.variable,t(stress.centrality(graph,subgraph.node.index))) # no 141-144
    # return.variable <- cbind(return.variable,t(weighted.stress.centrality(graph,subgraph.node.index))) # no 145-148
    print(148)
    return(return.variable)
} # end function


#=============================================================================================
#=============================================================================================
#  feature functions
#=============================================================================================
#=============================================================================================
# feature no. 1: number of edges
subgraph.num.edges <- function(subgraph){
    return(ecount(subgraph))
} # end function
# one function to calculate all
#=============================================================================================
# head feature: degree / calculating all mathematical moments
subgraph.degree <- function(subgraph){
    temp <- degree(subgraph)
    return(c(mean(temp),sd(temp),skewness(temp),kurtosis(temp)))
}
# feature no. 2: average degree
subgraph.average.degree <- function(subgraph){
    return(mean(degree(subgraph)))
} # end function
# feature no. 3: standard deviation degree
subgraph.sd.degree <- function(subgraph){
    return(sd(degree(subgraph)))
} # end function
# feature no. 4: skewness degree
subgraph.skewness.degree <- function(subgraph){
    return(skewness(degree(subgraph)))
} # end function
# feature no. 5: kurtosis degree
subgraph.kurtosis.degree <- function(subgraph){
    return(kurtosis(degree(subgraph)))
} # end function
#=============================================================================================
# head feature: weighted degree / calculating all mathematical moments
subgraph.weighted.degree <- function(subgraph){
    if(is.weighted(subgraph)){
        temp <- E(subgraph)$weight
        return(c(sum(temp),mean(temp),sd(temp),skewness(temp),kurtosis(temp)))
    }
    else{
        print('A non-weighted graph supplied! Returning NA')
        return(rep(NA,5))
    }
}
# feature no. 6 total weighted degree
subgraph.total.weighted.degree  <- function(subgraph){
    if(is.weighted(subgraph)){
        return(sum(E(subgraph)$weight))
    }
    else{
        print('A non-weighted graph supplied! Returning NA')
        return(NA)
    }
} # end function
# feature no. 7 average weighted degree
subgraph.average.weighted.degree <- function(subgraph){
    if(is.weighted(subgraph)){
        return(mean(E(subgraph)$weight))
    }
    else{
        print('A non-weighted graph supplied! Returning NA')
        return(NA)
    }
} # end function
# feature no. 8 standard deviation weighted degree
subgraph.sd.weighted.degree  <- function(subgraph){
    if(is.weighted(subgraph)){
        return(sd(E(subgraph)$weight))
    }
    else{
        print('A non-weighted graph supplied! Returning NA')
        return(NA)
    }
} # end function
# feature no. 9 skewness weighted degree
subgraph.skewness.weighted.degree <- function(subgraph){
    if(is.weighted(subgraph)){
        return(skewness(E(subgraph)$weight))
    }
    else{
        print('A non-weighted graph supplied!Returning NA')
        return(NA)
    }
} # end function
# feature no. 10 kurtosis weighted degree
subgraph.kurtosis.weighted.degree <- function(subgraph){
    if(is.weighted(subgraph)){
        return(kurtosis(E(subgraph)$weight))
    }
    else{
        print('A non-weighted graph supplied!Returning NA')
        return(NA)
    }
} # end function
#=============================================================================================
# head feature: weighted degree / calculating all mathematical moments
subgraph.abs.weighted.degree <- function(subgraph){
    if(is.weighted(subgraph)){
        temp <- abs(E(subgraph)$weight)
        return(c(sum(temp),mean(temp),sd(temp),skewness(temp),kurtosis(temp)))
    }
    else{
        print('A non-weighted graph supplied! Returning NA')
        return(rep(NA,5))
    }
}
# feature no.11 total weighted degree
subgraph.abs.total.weighted.degree  <- function(subgraph){
    if(is.weighted(subgraph)){
        return(sum(abs(E(subgraph)$weight)))
    }
    else{
        print('A non-weighted graph supplied! Returning NA')
        return(NA)
    }
} # end function
# feature no. 12 average weighted degree
subgraph.abs.average.weighted.degree <- function(subgraph){
    if(is.weighted(subgraph)){
        return(mean(abs(E(subgraph)$weight)))
    }
    else{
        print('A non-weighted graph supplied! Returning NA')
        return(NA)
    }
} # end function
# feature no. 13 standard deviation weighted degree
subgraph.abs.sd.weighted.degree  <- function(subgraph){
    if(is.weighted(subgraph)){
        return(sd(abs(E(subgraph)$weight)))
    }
    else{
        print('A non-weighted graph supplied! Returning NA')
        return(NA)
    }
} # end function
# feature no. 14 skewness weighted degree
subgraph.abs.skewness.weighted.degree <- function(subgraph){
    if(is.weighted(subgraph)){
        return(skewness(abs(E(subgraph)$weight)))
    }
    else{
        print('A non-weighted graph supplied!Returning NA')
        return(NA)
    }
} # end function
# feature no. 15 kurtosis weighted degree
subgraph.abs.kurtosis.weighted.degree <- function(subgraph){
    if(is.weighted(subgraph)){
        return(kurtosis(abs(E(subgraph)$weight)))
    }
    else{
        print('A non-weighted graph supplied!Returning NA')
        return(NA)
    }
} # end function
#=============================================================================================
# feature no. 16 geodesic distance of the subgraph - nonweigthed
geodesic.distance.subgraph <- function(subgraph){
    # return NA if single node subgraph
    if(vcount(subgraph)<2){
        print('No shortest path can be calculated for a single node graph. Returning NA.')
        return(NA)
    }
    # make sure to run a non-weighted graph
    t <- shortest.paths(subgraph,weights=NA)
    t <- triangle(t)
    # exclude all 0 and infinity length connections
    x <- which(t==0)
    x <- c(x,which(t==Inf))
    if(length(x)==0){ # no zero or infinit connection found
        return(min(t))
    } # end if
    else if(length(x)==(dim(t)[1]*dim(t)[2]))
    {
        return(Inf)
    } # end else if
    else{
        return(min(t[-x]))
    }
} # end function
# head feature: calculating 4 mathematical moments
geodesic.distance.subgraph.moments <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No shortest path can be calculated for a single node graph. Returning NA.')
        return(rep(NA,4))
    }
    # make sure to run a non-weighted graph
    t <- shortest.paths(subgraph,weights=NA)
    t <- triangle(t)
    x <- which(t==Inf)
    if(length(x)!=0)
        t <- t[-x]
    return(c(mean(t),sd(t),skewness(t),kurtosis(t)))
} # end function
# feature no. 17 average nonweighted path length of subgraph
average.geodesic.distance.subgraph <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No shortest path can be calculated for a single node graph. Returning NA.')
        return(NA)
    }
    # make sure to run a non-weigted graph
    return(average.path.length(subgraph,directed=F))
} # end function
# feature no. 18 sd nonweighted path length of subgraph
sd.geodesic.distance.subgraph <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No shortest path can be calculated for a single node graph. Returning NA.')
        return(NA)
    }
     # make sure to run a non-weighted graph
    t <- shortest.paths(subgraph,weights=NA)
    t <- triangle(t)
    x <- which(t==Inf)
    if(length(x)!=0)
        t <- t[-x]
    return(sd(t))
} # end function
# feature no. 19 skewness nonweighte path length of subgraph
skewness.geodesic.distance.subgraph <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No shortest path can be calculated for a single node graph. Returning NA.')
        return(NA)
    }
    # make sure to run a non-weighted graph
    t <- shortest.paths(subgraph,weights=NA)
    t <- triangle(t)
    x <- which(t==Inf)
    if(length(x)!=0)
        t <- t[-x]
    return(skewness(t))
} # end function
# feature no. 20 skewness nonweighte path length of subgraph
kurtosis.geodesic.distance.subgraph <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No shortest path can be calculated for a single node graph. Returning NA.')
        return(NA)
    }
    # make sure to run a non-weighted graph
    t <- shortest.paths(subgraph,weights=NA)
    t <- triangle(t)
    x <- which(t==Inf)
    if(length(x)!=0)
        t <- t[-x]
    return(kurtosis(t))
} # end function
#=============================================================================================
# head feature: calculating 4 mathematical moments
geodesic.distance.weighted.subgraph.moments <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No shortest path can be calculated for a single node graph. Returning NA.')
        return(rep(NA,5))
    }
    gr  <- subgraph
    E(gr)$weight  <- abs(E(gr)$weight)
    t <- shortest.paths(gr,weights=NULL)
    # only upper triangle
    t <- triangle(t)
    x  <- which(t==Inf)
    # if only non-connections found
    if(length(x)==length(t)){
        return(rep(Inf,4))
    }
    else if (length(x)==0){
        return(c(min(t),mean(t),sd(t),skewness(t),kurtosis(t)))
    } 
    else{
        return(c(min(t[-x]),mean(t[-x]),sd(t[-x]),skewness(t[-x]),kurtosis(t[-x])))
    }
} # end function
# feature no. 21 weighted geodesic distance of subgraph
geodesic.distance.weighted.subgraph <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No shortest path can be calculated for a single node graph. Returning NA.')
        return(NA)
    }
    gr  <- subgraph
    E(gr)$weight  <-  abs(E(gr)$weight)
    t <- shortest.paths(gr,weights = NULL)
    t <- triangle(t)
    x <- which(t==0)
    x <- c(x,which(t==Inf))
    if(length(x)==0){ # no zero or infinit connection found
        return(min(t))
    } # end if
    else if(length(x)==(dim(t)[1]*dim(t)[2]))
    {
        return(Inf)
    } # end else if
    else{
        return(min(t[-x]))
    }
} # end function

# feature no. 22 average weighted path length of subgraph
average.weighted.geodesic.subgraph <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No shortest path can be calculated for a single node graph. Returning NA.')
        return(NA)
    }
    gr  <- subgraph
    E(gr)$weight  <- abs(E(gr)$weight)
    t <- shortest.paths(gr,weights=NULL)
    # only upper triangle
    t <- triangle(t)
    x  <- which(t==Inf)
    # if only non-connections found
    if(length(x)==length(t)){
        return(Inf)
    }
    else if (length(x)==0){
        return(mean(t))
    } 
    else{
        return(mean(t[-x]))
    }
} # end function
# feature no. 23 sd weighted path length of subgraph
sd.weighted.geodesic.subgraph <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No shortest path can be calculated for a single node graph. Returning NA.')
        return(NA)
    }
    gr  <- subgraph
    E(gr)$weight  <- abs(E(gr)$weight)
    t <- shortest.paths(gr,weights=NULL)
    # only upper triangle
    t <- triangle(t)
    x  <- which(t==Inf)
    # if only non-connections found
    if(length(x)==length(t)){
        return(Inf)
    }
    else if (length(x)==0){
        return(sd(t))
    } 
    else{
        return(sd(t[-x]))
    }
} # end function
# feature no. 24 skewness weighted path length of subgraph
skewness.weighted.geodesic.subgraph <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No shortest path can be calculated for a single node graph. Returning NA.')
        return(NA)
    }
    gr  <- subgraph
    E(gr)$weight  <- abs(E(gr)$weight)
    t <- shortest.paths(gr,weights=NULL)
    # only upper triangle
    t <- triangle(t)
    x  <- which(t==Inf)
    # if only non-connections found
    if(length(x)==length(t)){
        return(Inf)
    }
    else if (length(x)==0){
        return(skewness(t))
    } 
    else{
        return(skewness(t[-x]))
    }
} # end function
# feature no. 25 sd weighted path length of subgraph
kurtosis.weighted.geodesic.subgraph <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No shortest path can be calculated for a single node graph. Returning NA.')
        return(NA)
    }
    gr  <- subgraph
    E(gr)$weight  <- abs(E(gr)$weight)
    t <- shortest.paths(gr,weights=NULL)
    # only upper triangle
    t <- triangle(t)
    x  <- which(t==Inf)
    # if only non-connections found
    if(length(x)==length(t)){
        return(Inf)
    }
    else if (length(x)==0){
        return(kurtosis(t))
    } 
    else{
        return(kurtosis(t[-x]))
    }
} # end function
#=============================================================================================
# feature no 26: compute how many geodesic distances (non-weighted) of the entire graph run through the subgraph
group.betweenness.subgraph.2 <- function(graph,subgraph.node.index){
   # print(paste(Sys.time(),"nonweighted: using loops",sep=" : "))
    node.index <- c(1:vcount(graph))
    node.index <- node.index[-subgraph.node.index]
    # loop over all nodes and find geodesic distances between them
    count.distances <- 0 # vector("list",0)
    for(i in 1:(length(node.index)-1)){
        #print(node.index[i])
        for(j in (i+1):(length(node.index))){
            # print(node.index[j])
            distances <- get.all.shortest.paths(graph,node.index[i],node.index[j],weights=NA,mode="all")$res
            # only if there is a connection between the two nodes
            if(length(distances)>0){
                for(k in 1:length(distances)){
                    m  <- match(subgraph.node.index,distances[[k]])
                    if (sum(is.na(m))!=length(m))  # that means it passes through one node of the subgraph
                        count.distances  <- count.distances + 1
                } # end nested nested for loop
            } # end if
        } # end nested loop
    } # end for loop
   # print(Sys.time())
    return(count.distances)
} # end function
group.betweenness.subgraph <- function(graph,subgraph.node.index){
   # print(paste(Sys.time(),"nonweighted: using lists",sep=" : "))
    node.index <- c(1:vcount(graph))
    node.index <- node.index[-subgraph.node.index]
    # loop over all nodes and find geodesic distances between them
    count.distances <- 0 # vector("list",0)
   
    for(i in 1:(length(node.index)-1)){
            distances <- get.all.shortest.paths(graph,node.index[i],node.index[-(1:i)],weights=NA,mode="all")$res
            # only if there is a connection between the two nodes
            if(length(distances)>0){
                # find all shortest paths
                m <- lapply(distances,match,subgraph.node.index)
                # calculate the length of each cell in the list 'm'
                length.m <- sapply(m,length)
                # which cells of the vector of each cell in the list are NAs
                # non-NAs indicate that the shortest paths goes through one of the nodes of the subgraph
                isna <- lapply(m,is.na)
                # sum up the number of NAs in each vector of the list
                sum.na <- sapply(isna,sum)
                # now count the instances where the length of the vectors is not equal to the number of NAs within
                count.distances <- count.distances + length(which((sum.na==length.m)==F))
            }
    } # end for loop
  #  print(Sys.time())
    return(count.distances)
} # end function
# feature no 27: compute how many weighted geodesic distances of the entire graph run through the subgraph
weighted.group.betweenness.subgraph.2 <- function(graph,subgraph.node.index){
   # print(paste(Sys.time(),"weighted: using loops",sep=" : "))
    options(warn = -1)
    node.index <- c(1:vcount(graph))
    node.index <- node.index[-subgraph.node.index]
    # set all weights to positive values
    E(graph)$weight  <- abs(E(graph)$weight)
    # loop over all nodes and find geodesic distances between them
    count.distances <- 0 # vector("list",0)
    for(i in 1:(length(node.index)-1)){
        #print(node.index[i])
        for(j in (i+1):(length(node.index))){
            # print(node.index[j])
            distances <- get.all.shortest.paths(graph,node.index[i],node.index[j],weights=NULL,mode="all")$res
            # only if there is a connection between the two nodes
            if(length(distances)>0){
                for(k in 1:length(distances)){
                    m  <- match(subgraph.node.index,distances[[k]])
                    if (sum(is.na(m))!=length(m))  # that means it passes through one node of the subgraph
                        count.distances  <- count.distances + 1
                } # end nested nested for loop
            } # end if
        } # end nested loop
    } # end for loop
    options(warn = 0)
   # print(Sys.time())
    return(count.distances)
} # end function
weighted.group.betweenness.subgraph <- function(graph,subgraph.node.index){
   # print(paste(Sys.time(),"weighted: using lists",sep=" : "))
    node.index <- c(1:vcount(graph))
    node.index <- node.index[-subgraph.node.index]
    E(graph)$weight  <- abs(E(graph)$weight)
    # loop over all nodes and find geodesic distances between them
    count.distances <- 0 # vector("list",0)
    
    for(i in 1:(length(node.index)-1)){
        distances <- get.all.shortest.paths(graph,node.index[i],node.index[-(1:i)],weights=NULL,mode="all")$res
        # only if there is a connection between the two nodes
        if(length(distances)>0){
            # find all shortest paths
            m <- lapply(distances,match,subgraph.node.index)
            # calculate the length of each cell in the list 'm'
            length.m <- sapply(m,length)
            # which cells of the vector of each cell in the list are NAs
            # non-NAs indicate that the shortest paths goes through one of the nodes of the subgraph
            isna <- lapply(m,is.na)
            # sum up the number of NAs in each vector of the list
            sum.na <- sapply(isna,sum)
            # now count the instances where the length of the vectors is not equal to the number of NAs within
            count.distances <- count.distances + length(which((sum.na==length.m)==F))
        }
    } # end for loop
   # print(Sys.time())
    return(count.distances)
} # end function
#=============================================================================================
# head feautre closeness centrality - calculates four mathematical moments
closeness.centrality.subgraph <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No closeness centrality can be calculated for a single node graph. Returning NA.')
        return(rep(NA,4))
    }
    cc <- closeness(subgraph,weights=NA)
    return(c(mean(cc),sd(cc),skewness(cc),kurtosis(cc)))
} # end function
# feature no. 28 average closeness centrality nonweighted
average.closeness.centrality.subgraph <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No closeness centrality can be calculated for a single node graph. Returning NA.')
        return(NA)
    }
    return(mean(closeness(subgraph,weights=NA)))
}
# feature no. 29 sd closeness centrality nonweighted
sd.closeness.centrality.subgraph <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No closeness centrality can be calculated for a single node graph. Returning NA.')
        return(NA)
    }
    return(sd(closeness(subgraph,weights=NA)))
}
# feature no. 30 skewness closeness centrality nonweighted
skewness.closeness.centrality.subgraph <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No closeness centrality can be calculated for a single node graph. Returning NA.')
        return(NA)
    }
    return(skewness(closeness(subgraph,weights=NA)))
}
# feature no. 31 kurtosis closeness centrality nonweighted
kurtosis.closeness.centrality.subgraph <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No closeness centrality can be calculated for a single node graph. Returning NA.')
        return(NA)
    }
    return(kurtosis(closeness(subgraph,weights=NA)))
}
# head feautre closeness centrality - calculates four mathematical moments
weighted.closeness.centrality.subgraph <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No closeness centrality can be calculated for a single node graph. Returning NA.')
        return(rep(NA,4))
    }
    E(subgraph)$weight  <- abs(E(subgraph)$weight)
    cc <- closeness(subgraph,weights=NULL)
    return(c(mean(cc),sd(cc),skewness(cc),kurtosis(cc)))
} # end function
# feature no. 32 average closeness centrality nonweighted
weighted.average.closeness.centrality.subgraph <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No closeness centrality can be calculated for a single node graph. Returning NA.')
        return(NA)
    }
    E(subgraph)$weight  <- abs(E(subgraph)$weight)
    return(mean(closeness(subgraph,weights=NULL)))
}
# feature no. 33 sd closeness centrality nonweighted
weighted.sd.closeness.centrality.subgraph <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No closeness centrality can be calculated for a single node graph. Returning NA.')
        return(NA)
    }
    E(subgraph)$weight  <- abs(E(subgraph)$weight)
    return(sd(closeness(subgraph,weights=NULL)))
}
# feature no. 34 skewness closeness centrality nonweighted
weighted.skewness.closeness.centrality.subgraph <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No closeness centrality can be calculated for a single node graph. Returning NA.')
        return(NA)
    }
    E(subgraph)$weight  <- abs(E(subgraph)$weight)
    return(skewness(closeness(subgraph,weights=NULL)))
}
# feature no. 35 kurtosis closeness centrality nonweighted
weighted.kurtosis.closeness.centrality.subgraph <- function(subgraph){
    if(vcount(subgraph)<2){
        print('No closeness centrality can be calculated for a single node graph. Returning NA.')
        return(NA)
    }
    E(subgraph)$weight  <- abs(E(subgraph)$weight)
    return(kurtosis(closeness(subgraph,weights=NULL)))
}
#=============================================================================================
# head feature closeness centrality subgraph within graph
closeness.centrality.graph.subgraph <- function(graph,subgraph.node.index){
    cc <- closeness(graph,subgraph.node.index,weights=NA)
    return(c(mean(cc),sd(cc),skewness(cc),kurtosis(cc)))
} # end function
# feature no. 36 average closeness centrality nonweighted
average.closeness.centrality.graph.subgraph <- function(graph,subgraph.node.index){
    return(mean(closeness(graph,subgraph.node.index,weights=NA)))
}
# feature no. 37 sd closeness centrality nonweighted
sd.closeness.centrality.graph.subgraph <- function(graph,subgraph.node.index){
    return(sd(closeness(graph,subgraph.node.index,weights=NA)))
}
# feature no. 38 skewness closeness centrality nonweighted
skewness.closeness.centrality.graph.subgraph <- function(graph,subgraph.node.index){
    return(skewness(closeness(graph,subgraph.node.index,weights=NA)))
}
# feature no. 39 kurtosis closeness centrality nonweighted
kurtosis.closeness.centrality.graph.subgraph <- function(graph,subgraph.node.index){
    return(kurtosis(closeness(graph,subgraph.node.index,weights=NA)))
}

# head feature closeness centrality subgraph within graph
weighted.closeness.centrality.graph.subgraph <- function(graph,subgraph.node.index){
    E(graph)$weight  <- abs(E(graph)$weight)
    cc <- closeness(graph,subgraph.node.index,weights=NULL)
    return(c(mean(cc),sd(cc),skewness(cc),kurtosis(cc)))
} # end function
# feature no. 40 average closeness centrality weighted from graph to subgraph
weighted.average.closeness.centrality.graph.subgraph <- function(graph,subgraph.node.index){
    E(graph)$weight  <- abs(E(graph)$weight)
    return(mean(closeness(graph,subgraph.node.index,weights=NULL)))
}
# feature no. 41 sd closeness centrality weighted from graph to subgraph
weighted.sd.closeness.centrality.graph.subgraph <- function(graph,subgraph.node.index){
    E(graph)$weight  <- abs(E(graph)$weight)
    return(sd(closeness(graph,subgraph.node.index,weights=NULL)))
}
# feature no. 42 skewness closeness centrality weighted from graph to subgraph
weighted.skewness.closeness.centrality.graph.subgraph <- function(graph,subgraph.node.index){
    E(graph)$weight  <- abs(E(graph)$weight)
    return(skewness(closeness(graph,subgraph.node.index,weights=NULL)))
}
# feature no. 43 kurtosis closeness centrality weighted from graph to subgraph
weighted.kurtosis.closeness.centrality.graph.subgraph <- function(graph,subgraph.node.index){
    E(graph)$weight  <- abs(E(graph)$weight)
    return(kurtosis(closeness(graph,subgraph.node.index,weights=NULL)))
}

#=============================================================================================
# feature no. 44 diameter of subgraph nonweighted
diameter.subgraph <- function(subgraph){
    if(vcount(subgraph)<2)
    {
        print("WARNING! Single-noded graph. Results = 0.")
    }
    return(diameter(subgraph,weights=NA,unconnected=TRUE))
} # end function
# feature no. 45 diameter of subgraph weighted
weighted.diameter.subgraph <- function(subgraph){
    if(vcount(subgraph)<2)
    {
        print("WARNING! Single-noded graph. Results = 0.")
    }
    E(subgraph)$weight  <- abs(E(subgraph)$weight)
    return(diameter(subgraph,weights=NULL,unconnected=TRUE))
} # end function
#=============================================================================================
# feature no. 46 diameter of graph through subgraph nonweighted = diameter centrality
diameter.centrality <- function(graph,subgraph.node.index){
    # first find all shortest paths between all nodex except the subgraph
    nodex <- c(1:vcount(graph))
    shortest.paths <- vector("list",0)
    
    # find all shortest paths
    for(i in 1:(length(nodex)-1)){
        sp <- get.all.shortest.paths(graph,nodex[i],nodex[-(1:i)],weights=NA,mode="all")$res
        if(length(sp)>0)
            shortest.paths[(length(shortest.paths)+1):(length(shortest.paths)+length(sp))] <- sp
    } # end for

    # now find the longest ones = diameter
    length.sps <- unlist(lapply(shortest.paths,length))
    diameters <- which(length.sps==max(length.sps))
    
    # now look over the actual path and determine wether it runs through the subgraph
    count.diameters <- 0
    if(length(diameters)>0){
       for(i in diameters){
           m  <- match(subgraph.node.index,shortest.paths[[i]])
           if (sum(is.na(m))!=length(m))  # that means it passes through one node of the subgraph
               count.diameters  <- count.diameters + 1
       } # end for loop
    } # end if
    return(count.diameters)
} # end function
# feature no. 47 diameter of graph through subgraph nweighted - diameter centrality
weighted.diameter.centrality <- function(graph,subgraph.node.index){
    options(warn=-1)
    # first find all shortest paths between all nodex except the subgraph
    nodex <- c(1:vcount(graph))
    shortest.paths <- vector("list",0)
    E(graph)$weight  <- abs(E(graph)$weight)
    diam <- diameter(graph,weights=NULL)
    # counting diameter variable
    count.diameters <- 0
    # find all shorest.paths
    for(i in 1:(length(nodex)-1)){
        for(j in (i+1):length(nodex)){
            sp <- get.all.shortest.paths(graph,nodex[i],nodex[j],weights=NULL,mode="all")$res
            dist <- distances(graph,nodex[i],nodex[j],weights=NULL,mode="all")
            if(dist!=diam){  # not equal in length
                next
            } # end if
            else{ # check if it runs through the subgraph
                m  <- match(subgraph.node.index,sp)
                if (sum(is.na(m))!=length(m)) { # that means it passes through one node of the subgraph
                    count.diameters  <- count.diameters + 1  
                } # end of 
            } # end else
        } # end nested for
    } # end for
    
    options(warn=0)
    return(count.diameters)
} # end function
#=============================================================================================
# head feature global clustering coefficent 48-52
clustering.coefficient.subgraph <- function(subgraph){
    gl.cc <- transitivity(subgraph,type="global",isolates="zero")
    loc.cc <- transitivity(subgraph,type="local",isolates="zero")
    # first one is global, the other local
    return(c(gl.cc,mean(loc.cc),sd(loc.cc),skewness(loc.cc),kurtosis(loc.cc)))
}
# feature no. 48 subgraph global clustering coefficient
global.clust.coef.subgraph <- function(subgraph){
    return(transitivity(subgraph,type="global",isolates="zero"))
} # end function   - ok
# feature no. 49 average local clustering coefficient of subgraph
average.local.clust.coef.subgraph <- function(subgraph){
    return(mean(transitivity(subgraph,type="local",isolates="zero")))
} # end function - ok
# feature no. 50 sd local clustering coefficient of subgraph
sd.local.clust.coef.subgraph <- function(subgraph){
    return(sd(transitivity(subgraph,type="local",isolates="zero",weights=NULL)))
} # end function - ok
# feature no. 51 skewness local clustering coefficient of subgraph
skewness.local.clust.coef.subgraph <- function(subgraph){
    return(skewness(transitivity(subgraph,type="local",isolates="zero")))
}# end function - ok
# feature no. 52 kurtosis local clustering coefficient of subgraph
kurtosis.local.clust.coef.subgraph <- function(subgraph){
    return(kurtosis(transitivity(subgraph,type="local",isolates="zero")))
}# end function - ok
# head feature weighted global clustering coefficient coefficient 53-57
weighted.clustering.coefficient.subgraph <- function(subgraph){
    gl.cc <- transitivity(subgraph,type="global",isolates="zero",weights=NULL)
    loc.cc <- transitivity(subgraph,type="local",isolates="zero",weights=NULL)
    # first one is global, the other local
    return(c(gl.cc,mean(loc.cc),sd(loc.cc),skewness(loc.cc),kurtosis(loc.cc)))
}
# feature no. 53 weighted subgraph global clustering coefficient
weighted.global.clust.coef.subgraph <- function(subgraph){
    return(transitivity(subgraph,type="global",isolates="zero",weights=NULL))
} # end function - ok
# feature no. 54 weighted average local clustering coefficient of subgraph
weighted.average.local.clust.coef.subgraph <- function(subgraph){
    return(mean(transitivity(subgraph,type="weighted",isolates="zero",weights=NULL)))
} # end function - ok
# feature no. 55 weighted sd local clustering coefficient of subgraph
weighted.sd.local.clust.coef.subgraph <- function(subgraph){
    return(sd(transitivity(subgraph,type="weighted",isolates="zero",weights=NULL)))
} # end function - ok
# feature no. 56 weighted skewness local clustering coefficient of subgraph
weighted.skewness.local.clust.coef.subgraph <- function(subgraph){
    return(skewness(transitivity(subgraph,type="weighted",isolates="zero",weights=NULL)))
} # end function - ok
# feature no. 57 kurtosis local clustering coefficient of subgraph
weighted.kurtosis.local.clust.coef.subgraph <- function(subgraph){
    return(kurtosis(transitivity(subgraph,type="weighted",isolates="zero",weights=NULL)))
} # end function - ok
#=============================================================================================
# head feature clustering coefficent 58-61
clustering.coefficient.subgraph.graph <- function(graph,subgraph.node.index){
    loc.cc <- transitivity(graph,type="local",isolates="zero",vids=subgraph.node.index)
    # f the other local
    return(c(mean(loc.cc),sd(loc.cc),skewness(loc.cc),kurtosis(loc.cc)))
}
# feature no. 58 average local clustering coefficient of subgraph within graph
average.local.clust.coef.subgraph.graph <- function(graph,subgraph.node.index){
    return(mean(transitivity(graph,type="local",isolates="zero",vids=subgraph.node.index)))
} # end function - ok
# feature no. 59 sd local clustering coefficient of subgraph within graph
sd.local.clust.coef.subgraph.graph <- function(graph,subgraph.node.index){
    return(sd(transitivity(graph,type="local",isolates="zero",vids=subgraph.node.index)))
} # end function - ok
# feature no. 60 skewness local clustering coefficient of subgraph within graph
skewness.local.clust.coef.subgraph.graph <- function(graph,subgraph.node.index){
    return(skewness(transitivity(graph,type="local",isolates="zero",vids=subgraph.node.index)))
} # end function - ok
# feature no. 61 kurtosis local clustering coefficient of subgraph within graph
kurtosis.local.clust.coef.subgraph.graph <- function(graph,subgraph.node.index){
    return(kurtosis(transitivity(graph,type="local",isolates="zero",vids=subgraph.node.index)))
} # end function - ok
# head feature weighted clustering coefficent 62-65
weighted.clustering.coefficient.subgraph.graph <- function(graph,subgraph.node.index){
    loc.cc <- transitivity(graph,type="local",isolates="zero",vids=subgraph.node.index,
                           weights=NULL)
    return(c(mean(loc.cc),sd(loc.cc),skewness(loc.cc),kurtosis(loc.cc)))
}
# feature no. 62 weighted average local clustering coefficient of subgraph within graph
weighted.average.local.clust.coef.subgraph.graph <- function(graph,subgraph.node.index){
    return(mean(transitivity(graph,type="weighted",isolates="zero",vids=subgraph.node.index, weights=NULL)))
} # end function - ok
# feature no. 63 weighted sd local clustering coefficient of subgraph within graph
weighted.sd.local.clust.coef.subgraph.graph <- function(graph,subgraph.node.index){
    return(sd(transitivity(graph,type="weighted",isolates="zero",vids=subgraph.node.index, weights=NULL)))
} # end function - ok
# feature no. 64 weighted skewness local clustering coefficient of subgraph within graph
weighted.skewness.local.clust.coef.subgraph.graph <- function(graph,subgraph.node.index){
    return(skewness(transitivity(graph,type="weighted",isolates="zero",vids=subgraph.node.index, weights=NULL)))
} # end function - ok
# feature no. 65 weighted kurtosis local clustering coefficient of subgraph within graph
weighted.kurtosis.local.clust.coef.subgraph.graph <- function(graph,subgraph.node.index){
    return(kurtosis(transitivity(graph,type="weighted",isolates="zero",vids=subgraph.node.index, weights=NULL)))
} # end function - ok
#=============================================================================================
# head feature node betweeness 66-69
node.betweenness.subgraph <- function(graph,subgraph.node.index){
    nb <- betweenness(graph,v=subgraph.node.index,weights=NA)
    return(c(mean(nb),sd(nb),skewness(nb),kurtosis(nb)))
} # end function - ok
# feature no. 66 average node betweenness of subgraph nodes within graph
average.node.betweenness.subgraph <- function(graph,subgraph.node.index){
    return(mean(betweenness(graph,v=subgraph.node.index,weights=NA)))
} # end function - ok
# feature no. 67 sd node betweenness of subgraph nodes within graph
sd.node.betweenness.subgraph <- function(graph,subgraph.node.index){
    return(sd(betweenness(graph,v=subgraph.node.index,weights=NA)))
} # end function - ok
# feature no. 68 skewness node betweenness of subgraph nodes within graph
skewness.node.betweenness.subgraph <- function(graph,subgraph.node.index){
    return(skewness(betweenness(graph,v=subgraph.node.index,weights=NA)))
} # end function - ok
# feature no. 69 kurtosis node betweenness of subgraph nodes within graph
kurtosis.node.betweenness.subgraph <- function(graph,subgraph.node.index){
    return(kurtosis(betweenness(graph,v=subgraph.node.index,weights=NA)))
} # end function - ok
# head feature node betweeness 70-73
weighted.node.betweenness.subgraph <- function(graph,subgraph.node.index){
    E(graph)$weight <- abs(E(graph)$weight)
    nb <- betweenness(graph,v=subgraph.node.index,weights=NULL)
    return(c(mean(nb),sd(nb),skewness(nb),kurtosis(nb)))
} # end function - ok
# feature no. 70 weighted average node betweenness of subgraph nodes within graph
weighted.average.node.betweenness.subgraph <- function(graph,subgraph.node.index){
    E(graph)$weight <- abs(E(graph)$weight)
    return(mean(betweenness(graph,v=subgraph.node.index,weights=NULL)))
} # end function - ok
# feature no. 71 weighted sd node betweenness of subgraph nodes within graph
weighted.sd.node.betweenness.subgraph <- function(graph,subgraph.node.index){
    E(graph)$weight <- abs(E(graph)$weight)
    return(sd(betweenness(graph,v=subgraph.node.index,weights=NULL)))
} # end function - ok
# feature no. 72 weighted skewness node betweenness of subgraph nodes within graph
weighted.skewness.node.betweenness.subgraph <- function(graph,subgraph.node.index){
    E(graph)$weight <- abs(E(graph)$weight)
    return(skewness(betweenness(graph,v=subgraph.node.index,weights=NULL)))
} # end function - ok
# feature no. 73 kurtosis node betweenness of subgraph nodes within graph
weighted.kurtosis.node.betweenness.subgraph <- function(graph,subgraph.node.index){
    E(graph)$weight <- abs(E(graph)$weight)
    return(kurtosis(betweenness(graph,v=subgraph.node.index,weights=NULL)))
} # end function - ok
#=============================================================================================
# head feature edge betweeness 74-77
edge.betweenness.subgraph <- function(graph,subgraph.edge.index){
    eb <- edge.betweenness(graph,e=subgraph.edge.index,weights=NA)
    return(c(mean(eb),sd(eb),skewness(eb),kurtosis(eb)))
} # end function  -  ok
# feature no. 74 average edge betweenness of subgraph edges within graph
average.edge.betweenness.subgraph <- function(graph,subgraph.edge.index){
    eb <- edge.betweenness(graph,e=subgraph.edge.index,weights=NA)
    if(length(eb)==0) return(NA)
    
    return(mean(eb))
    #return(mean(edge.betweenness(graph,e=subgraph.edge.index,weights=NA)))
    #  return(sd(edge.betweenness(graph,e=subgraph.edge.index,weights=NA)))
} # end function  - ok
# feature no. 75 sd edge betweenness of subgraph edges within graph
sd.edge.betweenness.subgraph <- function(graph,subgraph.edge.index){
    return(sd(edge.betweenness(graph,e=subgraph.edge.index,weights=NA)))
} # end function  -  ok
# feature no. 76 skewness edge betweenness of subgraph edges within graph
skewness.edge.betweenness.subgraph <- function(graph,subgraph.edge.index){
    return(skewness(edge.betweenness(graph,e=subgraph.edge.index,weights=NA)))
} # end function  -  ok
# feature no. 77 kurtosis edge betweenness of subgraph edges within graph
kurtosis.edge.betweenness.subgraph <- function(graph,subgraph.edge.index){
    return(kurtosis(edge.betweenness(graph,e=subgraph.edge.index,weights=NA)))
} # end function  -  ok
# head feature weighted edge betweeness 78-81
weighted.edge.betweenness.subgraph <- function(graph,subgraph.edge.index){
    E(graph)$weight <- abs(E(graph)$weight)
    eb <- edge.betweenness(graph,e=subgraph.edge.index,weights=NULL)
    return(c(mean(eb),sd(eb),skewness(eb),kurtosis(eb)))
} # end function  -  ok
# feature no. 78 weighted average edge betweenness of subgraph edges within graph
weighted.average.edge.betweenness.subgraph <- function(graph,subgraph.edge.index){
    E(graph)$weight <- abs(E(graph)$weight)
    return(mean(edge.betweenness(graph,e=subgraph.edge.index,weights=NULL)))
} # end function  -  ok
# feature no. 79 weighted sd edge betweenness of subgraph edges within graph
weighted.sd.edge.betweenness.subgraph <- function(graph,subgraph.edge.index){
    E(graph)$weight <- abs(E(graph)$weight)
    return(sd(edge.betweenness(graph,e=subgraph.edge.index,weights=NULL)))
} # end function  -  ok
# feature no. 80 weighted skewness node betweenness of subgraph nodes within graph
weighted.skewness.edge.betweenness.subgraph <- function(graph,subgraph.edge.index){
    E(graph)$weight <- abs(E(graph)$weight)
    return(skewness(edge.betweenness(graph,e=subgraph.edge.index,weights=NULL)))
} # end function  -  ok
# feature no. 81 kurtosis node betweenness of subgraph nodes within graph
weighted.kurtosis.edge.betweenness.subgraph <- function(graph,subgraph.edge.index){
    E(graph)$weight <- abs(E(graph)$weight)
    return(kurtosis(edge.betweenness(graph,e=subgraph.edge.index,weights=NULL)))
} # end function  -  ok
#=============================================================================================
# feature no. 82 assortativty of subgraph
assortativity.subgraph <- function(subgraph){
    if(vcount(subgraph)<2)
    {
        print("WARNING! Single-noded graph. Returning NA.")
        return(NA)
    }
    return(assortativity.degree(subgraph))
} # end function  -  ok
#=============================================================================================
# feature no. 83 density of subgraph
density.subgraph <- function(subgraph){
    if(vcount(subgraph)<2)
    {
        print("WARNING! Single-noded graph. Returning NA.")
        return(NA)
    }
    return(graph.density(simplify(subgraph),loops=F))
} # end function  -  ok
#=============================================================================================
# the following functions test whether the subgraph fall into the same community as computed by different community detecting algorithms
# that is, it returns the per centage (%) of the greatest portion of the subgraph belonging to the same community 

# head feature - general function to compute the greatest community
# this function takes in the results of the community detection and computes the greatest community
# as the community structure will remain the same for the entire graph and only the subgraph is
# subject to change
greatest.community.subgraph <- function(graph,subgraph.node.index,community){
    if(length(subgraph.node.index)<2){
        print("WARNING! Single-noded graph. Greatest community will be of length 1.")
        return(1)
    }
    sgc <- community$membership[subgraph.node.index]
    factors <- unique(sgc)
    
    greatest.community <- NULL
    for(i in 1:length(factors)){
        x <- length(which(sgc==factors[i]))
        if(i==1)
            greatest.community  <- x
        else if(x>greatest.community)
            greatest.community  <- x
    } # end for 
    return(greatest.community/length(subgraph.node.index))
} # end function
# function no. 84 fastgreedy.community nonweighted
fgc.subgraph <- function(graph, subgraph.node.index){
    fgc <- fastgreedy.community(graph,weights=NULL)
    # finding the membership of the subgraph community - sgc
    sgc <- fgc$membership[subgraph.node.index]
    factors <- unique(sgc)
    
    greatest.community <- NULL
    for(i in 1:length(factors)){
        x <- length(which(sgc==factors[i]))
        if(i==1)
            greatest.community  <- x
        else if(x>greatest.community)
            greatest.community  <- x
    } # end for 
    return(greatest.community/length(subgraph.node.index))
} # end function  -  ok
# function no. 85 fastgreedy.community weighted
fgc.weighted.subgraph <- function(graph, subgraph.node.index){
    fgc <- fastgreedy.community(graph,weights=abs(E(graph)$weight))
    # finding the membership of the subgraph community - sgc
    sgc <- fgc$membership[subgraph.node.index]
    factors <- unique(sgc)
    
    greatest.community <- NULL
    for(i in 1:length(factors)){
        x <- length(which(sgc==factors[i]))
        if(i==1)
            greatest.community  <- x
        else if(x>greatest.community)
            greatest.community  <- x
    } # end for 
    return(greatest.community/length(subgraph.node.index))
} # end function  -  ok
# function no. 86 walktrap.community nonweighted
wtc.subgraph <- function(graph, subgraph.node.index){
    wtc <- walktrap.community(graph,weights=NULL)
    # finding the membership of the subgraph community - sgc
    sgc <- wtc$membership[subgraph.node.index]
    factors <- unique(sgc)
    
    greatest.community <- NULL
    for(i in 1:length(factors)){
        x <- length(which(sgc==factors[i]))
        if(i==1)
            greatest.community  <- x
        else if(x>greatest.community)
            greatest.community  <- x
    } # end for 
    return(greatest.community/length(subgraph.node.index))
} # end function  -  ok
# function no. 87 walktrap.community weighted
wtc.weighted.subgraph <- function(graph, subgraph.node.index){
    wtc <- walktrap.community(graph,weights=E(graph)$weight)
    # finding the membership of the subgraph community - sgc
    sgc <- wtc$membership[subgraph.node.index]
    factors <- unique(sgc)
    
    greatest.community <- NULL
    for(i in 1:length(factors)){
        x <- length(which(sgc==factors[i]))
        if(i==1)
            greatest.community  <- x
        else if(x>greatest.community)
            greatest.community  <- x
    } # end for  
    return(greatest.community/length(subgraph.node.index))
} # end function  - ok
# function no. 88 edge.betweenness.community nonweighted
ebc.subgraph <- function(graph, subgraph.node.index){
    ebc <- edge.betweenness.community(graph,weights=NULL)
    # finding the membership of the subgraph community - sgc
    sgc <- ebc$membership[subgraph.node.index]
    factors <- unique(sgc)
    
    greatest.community <- NULL
    for(i in 1:length(factors)){
        x <- length(which(sgc==factors[i]))
        if(i==1)
            greatest.community  <- x
        else if(x>greatest.community)
            greatest.community  <- x
    } # end for 
    return(greatest.community/length(subgraph.node.index))
} # end function  -  ok
# function no. 89 edge.betweenness.community weighted
ebc.weighted.subgraph <- function(graph, subgraph.node.index){
    ebc <- edge.betweenness.community(graph,weights=abs(E(graph)$weight))
    # finding the membership of the subgraph community - sgc
    sgc <- ebc$membership[subgraph.node.index]
    factors <- unique(sgc)
    
    greatest.community <- NULL
    for(i in 1:length(factors)){
        x <- length(which(sgc==factors[i]))
        if(i==1)
            greatest.community  <- x
        else if(x>greatest.community)
            greatest.community  <- x
    } # end for 
    return(greatest.community/length(subgraph.node.index))
} # end function
# function no. 90 label.propagation.community nonweighted
lpc.subgraph <- function(graph, subgraph.node.index){
    lpc <- label.propagation.community(graph,weights=NA)
    # finding the membership of the subgraph community - sgc
    sgc <- lpc$membership[subgraph.node.index]
    factors <- unique(sgc)
    
    greatest.community <- NULL
    for(i in 1:length(factors)){
        x <- length(which(sgc==factors[i]))
        if(i==1)
            greatest.community  <- x
        else if(x>greatest.community)
            greatest.community  <- x
    } # end for 
    return(greatest.community/length(subgraph.node.index))
} # end function - ok
# function no. 91 label.propagation.community weighted
lpc.weighted.subgraph <- function(graph, subgraph.node.index){
    lpc <- label.propagation.community(graph,weights=abs(E(graph)$weight))
    # finding the membership of the subgraph community - sgc
    sgc <- lpc$membership[subgraph.node.index]
    factors <- unique(sgc)
    
    greatest.community <- NULL
    for(i in 1:length(factors)){
        x <- length(which(sgc==factors[i]))
        if(i==1)
            greatest.community  <- x
        else if(x>greatest.community)
            greatest.community  <- x
    } # end for 
    return(greatest.community/length(subgraph.node.index))
} # end function  -  ok
# function no. 92 leading.eigenvector.community nonweighted
lec.subgraph <- function(graph, subgraph.node.index){
    lec <- leading.eigenvector.community(graph,weights=NA)
    # finding the membership of the subgraph community - sgc
    sgc <- lec$membership[subgraph.node.index]
    factors <- unique(sgc)
    
    greatest.community <- NULL
    for(i in 1:length(factors)){
        x <- length(which(sgc==factors[i]))
        if(i==1)
            greatest.community  <- x
        else if(x>greatest.community)
            greatest.community  <- x
    } # end for 
    return(greatest.community/length(subgraph.node.index))
} # end function  -  ok
# function no. 93 leading.eigenvector.community weighted
lec.weighted.subgraph <- function(graph, subgraph.node.index){
    lec <- leading.eigenvector.community(graph,weights=abs(E(graph)$weight))
    # finding the membership of the subgraph community - sgc
    sgc <- lec$membership[subgraph.node.index]
    factors <- unique(sgc)
    
    greatest.community <- NULL
    for(i in 1:length(factors)){
        x <- length(which(sgc==factors[i]))
        if(i==1)
            greatest.community  <- x
        else if(x>greatest.community)
            greatest.community  <- x
    } # end for 
    return(greatest.community/length(subgraph.node.index))
} # end function  -  ok
# function no. 94 multilevel.community nonweighted
mlc.subgraph <- function(graph, subgraph.node.index){
    mlc <- multilevel.community(graph,weights=NA)
    # finding the membership of the subgraph community - sgc
    sgc <- mlc$membership[subgraph.node.index]
    factors <- unique(sgc)
    
    greatest.community <- NULL
    for(i in 1:length(factors)){
        x <- length(which(sgc==factors[i]))
        if(i==1)
            greatest.community  <- x
        else if(x>greatest.community)
            greatest.community  <- x
    } # end for 
    return(greatest.community/length(subgraph.node.index))
} # end function  -  ok
# function no. 95 multilevel.community weighted
mlc.weighted.subgraph <- function(graph, subgraph.node.index){
    mlc <- multilevel.community(graph,weights=abs(E(graph)$weight))
    # finding the membership of the subgraph community - sgc
    sgc <- mlc$membership[subgraph.node.index]
    factors <- unique(sgc)
    
    greatest.community <- NULL
    for(i in 1:length(factors)){
        x <- length(which(sgc==factors[i]))
        if(i==1)
            greatest.community  <- x
        else if(x>greatest.community)
            greatest.community  <- x
    } # end for 
    return(greatest.community/length(subgraph.node.index))
} # end function  - ok
#=============================================================================================
# friend measures
# feature no. 96-100 find all common friends of order 1 of nodes in subgraph, excluding friendship to each other
# if a subgraph with a single node is submitted, the nodes of order of the node will be considered
# as the subgraph and then nodes order will be calculated
# head feature friends = common,total,distinct, mixed friends, and preferential attachment score
friend.measures <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)>1){
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index) # returns a list
    }
    else{
        print("WARNING! Single-noded graph! Nodes of order 1 will be considered as subgraph")
        # first detect subgraph - order 1
        subgraph.node.index <- unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))
        # now detect order 1 neighborhood of new subgraph
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index)
    } 
    # common friends
    common.friends  <- 0
    # to calculate the preferential attachment score
    number.of.friends <- NULL
    if (length(nh)>0){
        # eliminate friendship to self and other nodes in the subgraph
        for(i in 1:length(nh)){
            nas <- match(nh[[i]],subgraph.node.index)
            idx <- which(!is.na(nas))
            nh[[i]] <- nh[[i]][-idx]
            number.of.friends <- rbind(number.of.friends,length(nh[[i]]))
        } # end for
        nh <- unlist(nh)
        common <- unique(nh)
        if (length(common)>0){
            for(j in 1:length(common)){
                x <- length(which(nh==common[j]))
                if(x/length(subgraph.node.index)==1){
                    # that means that the node here is common to all nodes in the subgraph
                    common.friends <- common.friends + 1
                }  # end if
            } # end for 
        } # end if
        # total friends
        common <- unique(nh)
        total.friends <- length(common)
        # distinct friends
        dup <- which(duplicated(nh)==T)
        m <- match(nh,nh[dup])
        distinct.friends <- length(which(is.na(m)))
        # mixed friends = total friends - common friends - distinct friends
        mixed.friends <- total.friends - common.friends - distinct.friends
        # preferential attachment score
        pas <- cumprod(number.of.friends)
        return(c(common.friends,total.friends,distinct.friends,mixed.friends,pas[length(pas)]))
    }
    else{
        return(rep(NA,5))
    }
} # end function
# feature 96
common.friends <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)>1){
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index) # returns a list
    }
    else{
        print("WARNING! Single-noded graph! Nodes of order 1 will be considered as subgraph")
        # first detect subgraph - order 1
        subgraph.node.index <- unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))
        # now detect order 1 neighborhood of new subgraph
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index)
    }
    common.friends  <- 0
    if (length(nh)>0){
        # eliminate friendship to self and other nodes in the subgraph
        for(i in 1:length(nh)){
            nas <- match(nh[[i]],subgraph.node.index)
            idx <- which(!is.na(nas))
            nh[[i]] <- nh[[i]][-idx]
        } # end for
        nh <- unlist(nh)
        # common friends
        common <- unique(nh)
        if (length(common)>0){
            for(j in 1:length(common)){
                x <- length(which(nh==common[j]))
                if(x/length(subgraph.node.index)==1){
                    # that means that the node here is common to all nodes in the subgraph
                    common.friends <- common.friends + 1
                }  # end if
            } # end for 
        } # end if
        # distinct friends
        return(common.friends)
    } # end if
    else
        return(NA)
} # end function  - ok
# feature no. 97 find all friends of order 1 of nodes in subgraph, excluding friendship to each other
# if a subgraph with a single node is submitted, the nodes of order of the node will be considered
# as the subgraph and then nodes order will be calculated
total.friends <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)>1){
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index) # returns a list
    }
    else{
        print("WARNING! Single-noded graph! Nodes of order 1 will be considered as subgraph")
        # first detect subgraph - order 1
        subgraph.node.index <- unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))
        # now detect order 1 neighborhood of new subgraph
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index)
    }

    if (length(nh)>0){
        # eliminate friendship to self and other nodes in the subgraph
        for(i in 1:length(nh)){
            nas <- match(nh[[i]],subgraph.node.index)
            idx <- which(!is.na(nas))
            nh[[i]] <- nh[[i]][-idx]
        } # end for
        nh <- unlist(nh)
        common <- unique(nh)
        return(length(common))
    } # end if
    else
        return(NA)
} # end function  -  ok
# feature no. 98 find all friends of order 1 of nodes in subgraph that are not in common, excluding friendship to each other
# if a subgraph with a single node is submitted, the nodes of order of the node will be considered
# as the subgraph and then nodes order will be calculated
distinct.friends <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)>1){
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index) # returns a list
    }
    else{
        print("WARNING! Single-noded graph! Nodes of order 1 will be considered as subgraph")
        # first detect subgraph - order 1
        subgraph.node.index <- unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))
        # now detect order 1 neighborhood of new subgraph
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index)
    }
    if (length(nh)>0){
        # eliminate friendship to self and other nodes in the subgraph
        for(i in 1:length(nh)){
            nas <- match(nh[[i]],subgraph.node.index)
            idx <- which(!is.na(nas))
            nh[[i]] <- nh[[i]][-idx]
        } # end for
        nh <- unlist(nh)
        dup <- which(duplicated(nh)==T)
        m <- match(nh,nh[dup])
        return(length(which(is.na(m))))
    } # end if
    else
        return(NA)
} # end function  - ok
# feature no. 99 mixed friends are common to 2 or more nodes in the subgraph, but not to all
# total friends - common friends - distinct friends
mixed.friends <- function(graph,subgraph.node.index){
    return(total.friends(graph,subgraph.node.index)-common.friends(graph,subgraph.node.index)-
               distinct.friends(graph,subgraph.node.index))
} # end function
# feature no. 100: preferential attachment score
# the preferential attachment score is the multiplication of all first neighbors for each node in the subgraph
preferential.attachment.score <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)>1){
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index) # returns a list
    }
    else{
        print("WARNING! Single-noded graph! Nodes of order 1 will be considered as subgraph")
        # first detect subgraph - order 1
        subgraph.node.index <- unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))
        # now detect order 1 neighborhood of new subgraph
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index)
    }
    
    number.of.friends <- NULL
    if (length(nh)>0){
        # eliminate friendship to self and other nodes in the subgraph
        for(i in 1:length(nh)){
            nas <- match(nh[[i]],subgraph.node.index)
            idx <- which(!is.na(nas))
            nh[[i]] <- nh[[i]][-idx]
            number.of.friends <- rbind(number.of.friends,length(nh[[i]]))
            #print(number.of.friends)
        } # end for
        prod <- cumprod(number.of.friends)
        return(prod[length(prod)])
    } # end if
    else
        return(NA)
} # end function - ok
# feature no. 101-105  -  ok
# the jaccard coefficient is the commonfriends measure / totalfriends measure
# this is usually applied to two nodes, but has been extended here to the entire subgraph
# head feature jaccard coef calculates jaccard coeff and four mathematical moments
# and is based on the friend measure which have been calculated before
jaccard.coefficient <- function(graph,subgraph.node.index,common.friends,total.friends){
    if(length(subgraph.node.index)==1){  # that means no neighborhood
        print("Node has no neighborhood. Returning NAs.")
        return(rep(NA,5))
        }
    jaccard.coef <- common.friends/total.friends
    jac.coef <- NULL
    for(i in 1:(length(subgraph.node.index)-1)){
        for(j in (i+1):length(subgraph.node.index)){
            jac.coef <- c(jac.coef,(common.friends(graph,c(i,j))/total.friends(graph,c(i,j))))
        } # end for
    } # end for
    return(c(jaccard.coef,mean(jac.coef),sd(jac.coef),skewness(jac.coef),kurtosis(jac.coef)))
} # end function
# feature 101
jaccard.coef <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)==1){  # that means no neighborhood
        print("Node has no neighborhood. Returning NAs.")
        return(NA)
    }
    return(common.friends(graph,subgraph.node.index)/total.friends(graph,subgraph.node.index))
} # end function  - ok
# feature no. 102  -  ok
# calculate the average jaccard coefficient for all posible pairs within the subgraph
average.jaccard.coef <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)==1){  # that means no neighborhood
        print("Node has no neighborhood. Returning NAs.")
        return(NA)
    }
    jac.coef <- NULL
    for(i in 1:(length(subgraph.node.index)-1)){
        for(j in (i+1):length(subgraph.node.index)){
            jac.coef <- c(jac.coef,(common.friends(graph,c(i,j))/total.friends(graph,c(i,j))))
        } # end for
    } # end for
    return(mean(jac.coef))
} # end function  - ok
# feature no. 103  -  ok
# calculate the standard deviation jaccard coefficient for all posible pairs within the subgraph
sd.jaccard.coef <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)==1){  # that means no neighborhood
        print("Node has no neighborhood. Returning NAs.")
        return(NA)
    }
    jac.coef <- NULL
    for(i in 1:(length(subgraph.node.index)-1)){
        for(j in (i+1):length(subgraph.node.index)){
            jac.coef <- c(jac.coef,(common.friends(graph,c(i,j))/total.friends(graph,c(i,j))))
        } # end for
    } # end for
    return(sd(jac.coef))
} # end function  - ok
# feature no. 104  -  ok
# calculate the skewness jaccard coefficient for all posible pairs within the subgraph
skewness.jaccard.coef <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)==1){  # that means no neighborhood
        print("Node has no neighborhood. Returning NAs.")
        return(NA)
    }
    jac.coef <- NULL
    for(i in 1:(length(subgraph.node.index)-1)){
        for(j in (i+1):length(subgraph.node.index)){
            jac.coef <- c(jac.coef,(common.friends(graph,c(i,j))/total.friends(graph,c(i,j))))
        } # end for
    } # end for
    return(skewness(jac.coef))
} # end function  - ok
# feature no. 105  -  ok
# calculate the kurtosis jaccard coefficient for all posible pairs within the subgraph
kurtosis.jaccard.coef <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)==1){  # that means no neighborhood
        print("Node has no neighborhood. Returning NAs.")
        return(NA)
    }
    jac.coef <- NULL
    for(i in 1:(length(subgraph.node.index)-1)){
        for(j in (i+1):length(subgraph.node.index)){
            jac.coef <- c(jac.coef,(common.friends(graph,c(i,j))/total.friends(graph,c(i,j))))
        } # end for
    } # end for
    return(kurtosis(jac.coef))
} # end function  - ok

#=============================================================================================
# the friends measure counts the number of nodes two nodes have in common incident on them
# regardless of the edge between them, if any at all
# since we are looking here at cliques, we are calculating average, sd, curtosis, and skewness
# for the subgraraph  
# head feature friends measure not be confused with friend measureS (see above)
# for subgraph
# features 106-110
friends.measure.subgraph <- function(subgraph){
    # check if there are more than one nodes in the subgraph
    if(vcount(subgraph)<2){
        print("WARNING! Single noded subgraph - friends measure cannot be calculated. Returning NA")
        return(rep(NA,5))
    }
    nh <- neighborhood(subgraph,order=1)#,mode <- "allnodes=subgraph.node.index)
    count <- NULL
    if (length(nh)>0){
        # eliminate friendship to self and other nodes in the subgraph
        for(i in 1:(length(nh)-1)){
            for(j in (i+1):length(nh)){
                if(length(nh[[i]][-1])>0 & length(nh[[j]][-1]>0)){
                    x <- match(nh[[i]][-1],nh[[j]][-1]) 
                    #print(x)
                    count <- c(count,length(which(!is.na(x))))
                    #count <- count + length(which(!is.na(x)))
                }
                else
                    #count <- count+0
                    count <- c(count,0)
            } # end nested loop
        } # end for
    } # end if
    return(c(sum(count),mean(count),sd(count),skewness(count),kurtosis(count)))
} # end function
# feature 106
friends.measure.subgraph.sum <- function(subgraph){
    if(vcount(subgraph)<2){
        print("WARNING! Single noded subgraph - friends measure cannot be calculated. Returning NA")
        return(rep(NA,1))
    }
    nh <- neighborhood(subgraph,order=1)#,mode <- "allnodes=subgraph.node.index)
    count <- NULL
    if (length(nh)>0){
        # eliminate friendship to self and other nodes in the subgraph
        for(i in 1:(length(nh)-1)){
            for(j in (i+1):length(nh)){
                if(length(nh[[i]][-1])>0 & length(nh[[j]][-1]>0)){
                    x <- match(nh[[i]][-1],nh[[j]][-1]) 
                    #print(x)
                    count <- c(count,length(which(!is.na(x))))
                    #count <- count + length(which(!is.na(x)))
                }
                else
                    #count <- count+0
                    count <- c(count,0)
            } # end nested loop
        } # end for
    } # end if
    return(sum(count))
}# end function  - ok
# feature 107
friends.measure.subgraph.avg <- function(subgraph.node){
    if(vcount(subgraph)<2){
        print("WARNING! Single noded subgraph - friends measure cannot be calculated. Returning NA")
        return(rep(NA,1))
    }
    nh <- neighborhood(subgraph,order=1)#,mode <- "allnodes=subgraph.node.index)
    count <- NULL
    if (length(nh)>0){
        # eliminate friendship to self and other nodes in the subgraph
        for(i in 1:(length(nh)-1)){
            for(j in (i+1):length(nh)){
                if(length(nh[[i]][-1])>0 & length(nh[[j]][-1]>0)){
                    x <- match(nh[[i]][-1],nh[[j]][-1]) 
                    # print(x)
                    count <- c(count,length(which(!is.na(x))))
                    #count <- count + length(which(!is.na(x)))
                }
                else
                    #count <- count+0
                    count <- c(count,0)
            } # end nested loop
        } # end for
    } # end if
    return(mean(count))
}# end function  - ok
# feature 108
friends.measure.subgraph.sd <- function(subgraph.node){
    if(vcount(subgraph)<2){
        print("WARNING! Single noded subgraph - friends measure cannot be calculated. Returning NA")
        return(rep(NA,1))
    }
    nh <- neighborhood(subgraph,order=1)#,mode <- "allnodes=subgraph.node.index)
    count <- NULL
    if (length(nh)>0){
        # eliminate friendship to self and other nodes in the subgraph
        for(i in 1:(length(nh)-1)){
            for(j in (i+1):length(nh)){
                if(length(nh[[i]][-1])>0 & length(nh[[j]][-1]>0)){
                    x <- match(nh[[i]][-1],nh[[j]][-1]) 
                    #print(x)
                    count <- c(count,length(which(!is.na(x))))
                    #count <- count + length(which(!is.na(x)))
                }
                else
                    #count <- count+0
                    count <- c(count,0)
            } # end nested loop
        } # end for
    } # end if
    return(sd(count))
}# end function  - ok
# feature 109
friends.measure.subgraph.skewness <- function(subgraph.node){
    if(vcount(subgraph)<2){
        print("WARNING! Single noded subgraph - friends measure cannot be calculated. Returning NA")
        return(rep(NA,1))
    }
    nh <- neighborhood(subgraph,order=1)#,mode <- "allnodes=subgraph.node.index)
    count <- NULL
    if (length(nh)>0){
        # eliminate friendship to self and other nodes in the subgraph
        for(i in 1:(length(nh)-1)){
            for(j in (i+1):length(nh)){
                if(length(nh[[i]][-1])>0 & length(nh[[j]][-1]>0)){
                    x <- match(nh[[i]][-1],nh[[j]][-1]) 
                    #print(x)
                    count <- c(count,length(which(!is.na(x))))
                    #count <- count + length(which(!is.na(x)))
                }
                else
                    #count <- count+0
                    count <- c(count,0)
            } # end nested loop
        } # end for
    } # end if
    return(skewness(count))
}# end function  - ok
# feature 110
friends.measure.subgraph.kurtosis <- function(subgraph.node){
    if(vcount(subgraph)<2){
        print("WARNING! Single noded subgraph - friends measure cannot be calculated. Returning NA")
        return(rep(NA,1))
    }
    nh <- neighborhood(subgraph,order=1)#,mode <- "allnodes=subgraph.node.index)
    count <- NULL
    if (length(nh)>0){
        # eliminate friendship to self and other nodes in the subgraph
        for(i in 1:(length(nh)-1)){
            for(j in (i+1):length(nh)){
                if(length(nh[[i]][-1])>0 & length(nh[[j]][-1]>0)){
                    x <- match(nh[[i]][-1],nh[[j]][-1]) 
                    #print(x)
                    count <- c(count,length(which(!is.na(x))))
                    #count <- count + length(which(!is.na(x)))
                }
                else
                    #count <- count+0
                    count <- c(count,0)
            } # end nested loop
        } # end for
    } # end if
    return(kurtosis(count))
}# end function  - ok
# feature 111-115
# for the entire graph
# head feature for friends measure of subgraph within the graph
friends.measure.graph <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)>1){
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index) # returns a list
    }
    else{
        print("WARNING! Single-noded graph! Nodes of order 1 will be considered as subgraph")
        # first detect subgraph - order 1
        subgraph.node.index <- unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))
        # now detect order 1 neighborhood of new subgraph
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index)
    }
    count <- NULL
    if (length(nh)>1){
        # eliminate friendship to self and other nodes in the subgraph
        for(i in 1:(length(nh)-1)){
            for(j in (i+1):length(nh)){
                if(length(nh[[i]][-1])>0 & length(nh[[j]][-1]>0)){
                    x <- match(nh[[i]][-1],nh[[j]][-1]) 
                    #print(x)
                    count <- c(count,length(which(!is.na(x))))
                    #count <- count + length(which(!is.na(x)))
                }
                else
                    #count <- count+0
                    count <- c(count,0)
            } # end nested loop
        } # end for
        return(c(sum(count),mean(count),sd(count),skewness(count),kurtosis(count)))
    } # end if else
    return(rep(0,5))
} # end function
# feature 111
friends.measure.graph.sum <- function(graph, subgraph.node.index){
    if(length(subgraph.node.index)>1){
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index) # returns a list
    }
    else{
        print("WARNING! Single-noded graph! Nodes of order 1 will be considered as subgraph")
        # first detect subgraph - order 1
        subgraph.node.index <- unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))
        # now detect order 1 neighborhood of new subgraph
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index)
    }
    count <- NULL
    if (length(nh)>1){
        # eliminate friendship to self and other nodes in the subgraph
        for(i in 1:(length(nh)-1)){
            for(j in (i+1):length(nh)){
                if(length(nh[[i]][-1])>0 & length(nh[[j]][-1]>0)){
                    x <- match(nh[[i]][-1],nh[[j]][-1]) 
                    #print(x)
                    count <- c(count,length(which(!is.na(x))))
                    #count <- count + length(which(!is.na(x)))
                }
                else
                    #count <- count+0
                    count <- c(count,0)
            } # end nested loop
        } # end for
        return(sum(count))
    } # end if
    return(0)
}# end function  - ok
# feature 112
friends.measure.graph.avg <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)>1){
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index) # returns a list
    }
    else{
        print("WARNING! Single-noded graph. Nodes of order 1 will be considered as subgraph")
        # first detect subgraph - order 1
        subgraph.node.index <- unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))
        # now detect order 1 neighborhood of new subgraph
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index)
    }
    count <- NULL
    if (length(nh)>1){
        # eliminate friendship to self and other nodes in the subgraph
        for(i in 1:(length(nh)-1)){
            for(j in (i+1):length(nh)){
                if(length(nh[[i]][-1])>0 & length(nh[[j]][-1]>0)){
                    x <- match(nh[[i]][-1],nh[[j]][-1]) 
                    # print(x)
                    count <- c(count,length(which(!is.na(x))))
                    #count <- count + length(which(!is.na(x)))
                }
                else
                    #count <- count+0
                    count <- c(count,0)
            } # end nested loop
        } # end for
        return(mean(count))
    } # end if
    return(0)
}# end function  - ok
# feature 113
friends.measure.graph.sd <- function(graph, subgraph.node.index){
    if(length(subgraph.node.index)>1){
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index) # returns a list
    }
    else{
        print("WARNING! Single-noded graph! Nodes of order 1 will be considered as subgraph")
        # first detect subgraph - order 1
        subgraph.node.index <- unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))
        # now detect order 1 neighborhood of new subgraph
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index)
    }
    count <- NULL
    if (length(nh)>1){
        # eliminate friendship to self and other nodes in the subgraph
        for(i in 1:(length(nh)-1)){
            for(j in (i+1):length(nh)){
                if(length(nh[[i]][-1])>0 & length(nh[[j]][-1]>0)){
                    x <- match(nh[[i]][-1],nh[[j]][-1]) 
                    #print(x)
                    count <- c(count,length(which(!is.na(x))))
                    #count <- count + length(which(!is.na(x)))
                }
                else
                    #count <- count+0
                    count <- c(count,0)
            } # end nested loop
        } # end for
        return(sd(count))
    } # end if
    return(0)
}# end function  - ok
# feature 114
friends.measure.graph.skewness <- function(graph, subgraph.node.index){
    if(length(subgraph.node.index)>1){
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index) # returns a list
    }
    else{
        print("WARNING! Single-noded graph! Nodes of order 1 will be considered as subgraph")
        # first detect subgraph - order 1
        subgraph.node.index <- unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))
        # now detect order 1 neighborhood of new subgraph
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index)
    }
    count <- NULL
    if (length(nh)>1){
        # eliminate friendship to self and other nodes in the subgraph
        for(i in 1:(length(nh)-1)){
            for(j in (i+1):length(nh)){
                if(length(nh[[i]][-1])>0 & length(nh[[j]][-1]>0)){
                    x <- match(nh[[i]][-1],nh[[j]][-1]) 
                    #print(x)
                    count <- c(count,length(which(!is.na(x))))
                    #count <- count + length(which(!is.na(x)))
                }
                else
                    #count <- count+0
                    count <- c(count,0)
            } # end nested loop
        } # end for
        return(skewness(count))
    } # end if
    return(0)
}# end function  - ok
# feature 115
friends.measure.graph.kurtosis <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)>1){
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index) # returns a list
    }
    else{
        print("WARNING! Single-noded graph! Nodes of order 1 will be considered as subgraph")
        # first detect subgraph - order 1
        subgraph.node.index <- unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))
        # now detect order 1 neighborhood of new subgraph
        nh <- neighborhood(graph,order=1,nodes=subgraph.node.index)
    }
    count <- NULL
    if (length(nh)>1){
        # eliminate friendship to self and other nodes in the subgraph
        for(i in 1:(length(nh)-1)){
            for(j in (i+1):length(nh)){
                if(length(nh[[i]][-1])>0 & length(nh[[j]][-1]>0)){
                    x <- match(nh[[i]][-1],nh[[j]][-1]) 
                    #print(x)
                    count <- c(count,length(which(!is.na(x))))
                    #count <- count + length(which(!is.na(x)))
                }
                else
                    #count <- count+0
                    count <- c(count,0)
            } # end nested loop
        } # end for
        return(kurtosis(count))
    } # end if
    return(0)
}# end function  - ok
#=============================================================================================
# feature no. 116: number of edges of all edges incident on nodes in the subgraph
num.edges.graph <- function(graph,subgraph.node.index){
    neighborhood <- sort(unique(unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))))
    if(length(neighborhood)==1){
        print("Node has no neighborhood. Thus, no edges incident on it. Returning 0")
        return(0)
    }
    sni <- NULL
    for(j in 1:(length(neighborhood)-1)){
        for(k in (j+1):length(neighborhood)){
            sni <- c(sni,neighborhood[j],neighborhood[k])
        } # end for
    } # end for
    nonzero <- which(get.edge.ids(graph,sni)!=0)
    eids <- get.edge.ids(graph, sni)[nonzero]
    gr <- subgraph.edges(graph,eids,delete.vertices=T)
    return(ecount(gr))
} # end function
# head feature degree
degree.graph <- function(graph,subgraph.node.index){
    deg <- degree(graph)[subgraph.node.index]
    return(c(mean(deg),sd(deg),skewness(deg),kurtosis(deg)))
}
# feature no. 117: average degree
average.degree.graph <- function(graph,subgraph.node.index){
    return(mean(degree(graph)[subgraph.node.index]))
} # end function
# feature no.118: standard deviation degree
sd.degree.graph <- function(graph,subgraph.node.index){
    return(sd(degree(graph)[subgraph.node.index]))
} # end function
# feature no. 119: skewness degree
skewness.degree.graph <- function(graph,subgraph.node.index){
    return(skewness(degree(graph)[subgraph.node.index]))
} # end function
# feature no. 120: kurtosis degree
kurtosis.degree.graph <- function(graph,subgraph.node.index){
    return(kurtosis(degree(graph)[subgraph.node.index]))
} # end function
#=============================================================================================
# head feature : weighted degree graph
weighted.degree.graph <- function(graph,subgraph.node.index){
    neighborhood <- sort(unique(unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))))
    if(length(neighborhood)==1){
        print("Node has no neighborhood. Thus, no edges incident on it. Returning 0")
        return(rep(0,5))
    }
    sni <- NULL
    for(j in 1:(length(neighborhood)-1)){
        for(k in (j+1):length(neighborhood)){
            sni <- c(sni,neighborhood[j],neighborhood[k])
        } # end for
    } # end for
    nonzero <- which(get.edge.ids(graph,sni)!=0)
    eids <- get.edge.ids(graph, sni)[nonzero]
    gr <- subgraph.edges(graph,eids,delete.vertices=T)
    return(c(sum(E(gr)$weight),mean(E(gr)$weight),sd(E(gr)$weight),
             skewness(E(gr)$weight),kurtosis(E(gr)$weight)))
} # end function
# feature no. 121 total weighted degree
total.weighted.degree.graph  <- function(graph,subgraph.node.index){
    neighborhood <- sort(unique(unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))))
    sni <- NULL
    for(j in 1:(length(neighborhood)-1)){
        for(k in (j+1):length(neighborhood)){
            sni <- c(sni,neighborhood[j],neighborhood[k])
        } # end for
    } # end for
    nonzero <- which(get.edge.ids(graph,sni)!=0)
    eids <- get.edge.ids(graph, sni)[nonzero]
    gr <- subgraph.edges(graph,eids,delete.vertices=T)
    return(sum(E(gr)$weight))
} # end function
# feature no. 122 average weighted degree
average.weighted.degree.graph <- function(graph,subgraph.node.index){
    neighborhood <- sort(unique(unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))))
    sni <- NULL
    for(j in 1:(length(neighborhood)-1)){
        for(k in (j+1):length(neighborhood)){
            sni <- c(sni,neighborhood[j],neighborhood[k])
        } # end for
    } # end for
    nonzero <- which(get.edge.ids(graph,sni)!=0)
    eids <- get.edge.ids(graph, sni)[nonzero]
    gr <- subgraph.edges(graph,eids,delete.vertices=T)
    return(mean(E(gr)$weight))
} # end function
# feature no. 123 standard deviation weighted degree
sd.weighted.degree.graph  <- function(graph,subgraph.node.index){
    neighborhood <- sort(unique(unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))))
    sni <- NULL
    for(j in 1:(length(neighborhood)-1)){
        for(k in (j+1):length(neighborhood)){
            sni <- c(sni,neighborhood[j],neighborhood[k])
        } # end for
    } # end for
    nonzero <- which(get.edge.ids(graph,sni)!=0)
    eids <- get.edge.ids(graph, sni)[nonzero]
    gr <- subgraph.edges(graph,eids,delete.vertices=T)
    return(sd(E(gr)$weight))
} # end function
# feature no. 124 skewness weighted degree
skewness.weighted.degree.graph <- function(graph,subgraph.node.index){
    neighborhood <- sort(unique(unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))))
    sni <- NULL
    for(j in 1:(length(neighborhood)-1)){
        for(k in (j+1):length(neighborhood)){
            sni <- c(sni,neighborhood[j],neighborhood[k])
        } # end for
    } # end for
    nonzero <- which(get.edge.ids(graph,sni)!=0)
    eids <- get.edge.ids(graph, sni)[nonzero]
    gr <- subgraph.edges(graph,eids,delete.vertices=T)
    return(skewness(E(gr)$weight))
} # end function
# feature no. 125 kurtosis weighted degree
kurtosis.weighted.degree.graph <- function(graph,subgraph.node.index){
    neighborhood <- sort(unique(unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))))
    sni <- NULL
    for(j in 1:(length(neighborhood)-1)){
        for(k in (j+1):length(neighborhood)){
            sni <- c(sni,neighborhood[j],neighborhood[k])
        } # end for
    } # end for
    nonzero <- which(get.edge.ids(graph,sni)!=0)
    eids <- get.edge.ids(graph, sni)[nonzero]
    gr <- subgraph.edges(graph,eids,delete.vertices=T)
    return(kurtosis(E(gr)$weight))
} # end function
#=============================================================================================
# head feature: absolute weighted degree graph
abs.weighted.degree.graph <- function(graph,subgraph.node.index){
    neighborhood <- sort(unique(unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))))
    if(length(neighborhood)==1){
        print("Node has no neighborhood. Thus, no edges incident on it. Returning 0")
        return(rep(0,5))
    }
    sni <- NULL
    for(j in 1:(length(neighborhood)-1)){
        for(k in (j+1):length(neighborhood)){
            sni <- c(sni,neighborhood[j],neighborhood[k])
        } # end for
    } # end for
    nonzero <- which(get.edge.ids(graph,sni)!=0)
    eids <- get.edge.ids(graph, sni)[nonzero]
    gr <- subgraph.edges(graph,eids,delete.vertices=T)
    E(gr)$weight <- abs(E(gr)$weight)
    return(c(sum(E(gr)$weight),mean(E(gr)$weight),sd(E(gr)$weight),
             skewness(E(gr)$weight),kurtosis(E(gr)$weight)))
} # end function
# feature no. 126 total absolute weighted degree
total.abs.weighted.degree.graph  <- function(graph,subgraph.node.index){
    neighborhood <- sort(unique(unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))))
    sni <- NULL
    for(j in 1:(length(neighborhood)-1)){
        for(k in (j+1):length(neighborhood)){
            sni <- c(sni,neighborhood[j],neighborhood[k])
        } # end for
    } # end for
    nonzero <- which(get.edge.ids(graph,sni)!=0)
    eids <- get.edge.ids(graph, sni)[nonzero]
    gr <- subgraph.edges(graph,eids,delete.vertices=T)
    return(sum(abs(E(gr)$weight)))
    
} # end function
# feature no. 127 average absolute weighted degree
average.abs.weighted.degree.graph <- function(graph,subgraph.node.index){
    neighborhood <- sort(unique(unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))))
    sni <- NULL
    for(j in 1:(length(neighborhood)-1)){
        for(k in (j+1):length(neighborhood)){
            sni <- c(sni,neighborhood[j],neighborhood[k])
        } # end for
    } # end for
    nonzero <- which(get.edge.ids(graph,sni)!=0)
    eids <- get.edge.ids(graph, sni)[nonzero]
    gr <- subgraph.edges(graph,eids,delete.vertices=T)
    return(mean(abs(E(gr)$weight)))
} # end function
# feature no. 128 standard deviation absolute weighted degree
sd.abs.weighted.degree.graph  <- function(graph,subgraph.node.index){
    neighborhood <- sort(unique(unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))))
    sni <- NULL
    for(j in 1:(length(neighborhood)-1)){
        for(k in (j+1):length(neighborhood)){
            sni <- c(sni,neighborhood[j],neighborhood[k])
        } # end for
    } # end for
    nonzero <- which(get.edge.ids(graph,sni)!=0)
    eids <- get.edge.ids(graph, sni)[nonzero]
    gr <- subgraph.edges(graph,eids,delete.vertices=T)
    return(sd(abs(E(gr)$weight)))
} # end function
# feature no. 129 skewness weighted degree
skewness.abs.weighted.degree.graph <- function(graph,subgraph.node.index){
    neighborhood <- sort(unique(unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))))
    sni <- NULL
    for(j in 1:(length(neighborhood)-1)){
        for(k in (j+1):length(neighborhood)){
            sni <- c(sni,neighborhood[j],neighborhood[k])
        } # end for
    } # end for
    nonzero <- which(get.edge.ids(graph,sni)!=0)
    eids <- get.edge.ids(graph, sni)[nonzero]
    gr <- subgraph.edges(graph,eids,delete.vertices=T)
    return(skewness(abs(E(gr)$weight)))
} # end function
# feature no. 130 kurtosis weighted degree
kurtosis.abs.weighted.degree.graph <- function(graph,subgraph.node.index){
    neighborhood <- sort(unique(unlist(neighborhood(graph,order=1,nodes=subgraph.node.index))))
    sni <- NULL
    for(j in 1:(length(neighborhood)-1)){
        for(k in (j+1):length(neighborhood)){
            sni <- c(sni,neighborhood[j],neighborhood[k])
        } # end for
    } # end for
    nonzero <- which(get.edge.ids(graph,sni)!=0)
    eids <- get.edge.ids(graph, sni)[nonzero]
    gr <- subgraph.edges(graph,eids,delete.vertices=T)
    return(kurtosis(abs(E(gr)$weight)))
} # end function
#=============================================================================================
geodesic.distances.graph <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)<2){
        print("WARNING! Single-noded graph. Geodesic distance cannot be computed. Returning NA")
        return(rep(NA,5))
    }
    # make sure to run a non-weighted graph
    t <- shortest.paths(graph,weights=NA)
    # only for subgraph nodes
    t <- t[subgraph.node.index,subgraph.node.index]
    t <- triangle(t)
    # exclude all 0 and infinity length connections
    x <- which(t==0)
    x <- c(x,which(t==Inf))
    if(length(x)==0){ # no zero or infinit connection found
        return(c(min(t),mean(t),sd(t),skewness(t),kurtosis(t)))
    } # end if
    else if(length(x)==(dim(t)[1]*dim(t)[2]))
    {
        return(rep(Inf,5))
    } # end else if
    else{
        return(c(min(t[-x]),mean(t[-x]),sd(t[-x]),skewness(t[-x]),kurtosis(t[-x])))
    }
} # end function
# feature no. 131 geodesic distance of the graph - nonweigthed
geodesic.distance.graph <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)<2){
        print("WARNING! Single-noded graph. Geodesic distance cannot be computed. Returning NA")
        return(NA)
    }
    # make sure to run a non-weighted graph
    t <- shortest.paths(graph,weights=NA)
    # only for subgraph nodes
    t <- t[subgraph.node.index,subgraph.node.index]
    t <- triangle(t)
    # exclude all 0 and infinity length connections
    x <- which(t==0)
    x <- c(x,which(t==Inf))
    if(length(x)==0){ # no zero or infinit connection found
        return(min(t))
    } # end if
    else if(length(x)==(dim(t)[1]*dim(t)[2]))
    {
        return(Inf)
    } # end else if
    else{
        return(min(t[-x]))
    }
} # end function
# feature no. 132 average nonweighted path length of graph
average.geodesic.distance.graph <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)<2){
        print("WARNING! Single-noded graph. Geodesic distance cannot be computed. Returning NA")
        return(NA)
    }
    # make sure to run a non-weigted graph
    t <- shortest.paths(graph,weights=NA)
    t <- t[subgraph.node.index,subgraph.node.index]
    t <- triangle(t)
    # exclude all 0 and infinity length connections
    x <- which(t==0)
    x <- c(x,which(t==Inf))
    if(length(x)==0){ # no zero or infinit connection found
        return(mean(t))
    } # end if
    else if(length(x)==(dim(t)[1]*dim(t)[2]))
    {
        return(Inf)
    } # end else if
    else{
        return(mean(t[-x]))
    }
} # end function
# feature no. 133 sd nonweighted path length of graph
sd.geodesic.distance.graph <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)<2){
        print("WARNING! Single-noded graph. Geodesic distance cannot be computed. Returning NA")
        return(NA)
    }
    # make sure to run a non-weighted graph
    t <- shortest.paths(graph,weights=NA)
    t <- t[subgraph.node.index,subgraph.node.index]
    t <- triangle(t)
    # exclude all 0 and infinity length connections
    x <- which(t==0)
    x <- c(x,which(t==Inf))
    if(length(x)==0){ # no zero or infinit connection found
        return(sd(t))
    } # end if
    else if(length(x)==(dim(t)[1]*dim(t)[2]))
    {
        return(Inf)
    } # end else if
    else{
        return(sd(t[-x]))
    }
} # end function
# feature no. 134 skewness nonweighte path length of graph
skewness.geodesic.distance.graph <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)<2){
        print("WARNING! Single-noded graph. Geodesic distance cannot be computed. Returning NA")
        return(NA)
    }
    # make sure to run a non-weighted graph
    t <- shortest.paths(graph,weights=NA)
    t <- t[subgraph.node.index,subgraph.node.index]
    t <- triangle(t)
    # exclude all 0 and infinity length connections
    x <- which(t==0)
    x <- c(x,which(t==Inf))
    if(length(x)==0){ # no zero or infinit connection found
        return(skewness(t))
    } # end if
    else if(length(x)==(dim(t)[1]*dim(t)[2]))
    {
        return(Inf)
    } # end else if
    else{
        return(skewness(t[-x]))
    }
} # end function
# feature no. 135 skewness nonweighte path length of graph
kurtosis.geodesic.distance.graph <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)<2){
        print("WARNING! Single-noded graph. Geodesic distance cannot be computed. Returning NA")
        return(NA)
    }
    # make sure to run a non-weighted graph
    t <- shortest.paths(graph,weights=NA)
    t <- t[subgraph.node.index,subgraph.node.index]
    t <- triangle(t)
    # exclude all 0 and infinity length connections
    x <- which(t==0)
    x <- c(x,which(t==Inf))
    if(length(x)==0){ # no zero or infinit connection found
        return(kurtosis(t))
    } # end if
    else if(length(x)==(dim(t)[1]*dim(t)[2]))
    {
        return(Inf)
    } # end else if
    else{
        return(kurtosis(t[-x]))
    }
} # end function

#=============================================================================================
# head feature weighted geodesic distances of graph
geodesic.distances.weighted.graph <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)<2){
        print("WARNING! Single-noded graph. Geodesic distance cannot be computed. Returning NA")
        return(rep(NA,5))
    }
    gr  <- graph
    E(gr)$weight  <-  abs(E(gr)$weight)
    t <- shortest.paths(gr,weights = NULL)
    t <- t[subgraph.node.index,subgraph.node.index]
    t <- triangle(t)
    x <- which(t==0)
    x <- c(x,which(t==Inf))
    if(length(x)==0){ # no zero or infinit connection found
        return(c(min(t),mean(t),sd(t),skewness(t),kurtosis(t)))
    } # end if
    else if(length(x)==(dim(t)[1]*dim(t)[2]))
    {
        return(rep(Inf,5))
    } # end else if
    else{
        return(c(min(t[-x]),mean(t[-x]),sd(t[-x]),skewness(t[-x]),kurtosis(t[-x])))
    }
} # end function
# feature no. 136 weighted geodesic distance of graph
geodesic.distance.weighted.graph <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)<2){
        print("WARNING! Single-noded graph. Geodesic distance cannot be computed. Returning NA")
        return(NA)
    }
    gr  <- graph
    E(gr)$weight  <-  abs(E(gr)$weight)
    t <- shortest.paths(gr,weights = NULL)
    t <- t[subgraph.node.index,subgraph.node.index]
    t <- triangle(t)
    x <- which(t==0)
    x <- c(x,which(t==Inf))
    if(length(x)==0){ # no zero or infinit connection found
        return(min(t))
    } # end if
    else if(length(x)==(dim(t)[1]*dim(t)[2]))
    {
        return(Inf)
    } # end else if
    else{
        return(min(t[-x]))
    }
} # end function
# feature no. 137 average weighted path length of graph
average.weighted.geodesic.graph <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)<2){
        print("WARNING! Single-noded graph. Geodesic distance cannot be computed. Returning NA")
        return(NA)
    }
    gr  <- graph
    E(gr)$weight  <- abs(E(gr)$weight)
    t <- shortest.paths(gr,weights=NULL)
    t <- t[subgraph.node.index,subgraph.node.index]
    
    # only upper triangle
    t <- triangle(t)
    x  <- which(t==Inf)
    # if only non-connections found
    if(length(x)==length(t)){
        return(Inf)
    }
    else if (length(x)==0){
        return(mean(t))
    } 
    else{
        return(mean(t[-x]))
    }
} # end function
# feature no. 138 sd weighted path length of graph
sd.weighted.geodesic.graph <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)<2){
        print("WARNING! Single-noded graph. Geodesic distance cannot be computed. Returning NA")
        return(NA)
    }
    gr  <- graph
    E(gr)$weight  <- abs(E(gr)$weight)
    t <- shortest.paths(gr,weights=NULL)
    t <- t[subgraph.node.index,subgraph.node.index]
    # only upper triangle
    t <- triangle(t)
    x  <- which(t==Inf)
    # if only non-connections found
    if(length(x)==length(t)){
        return(Inf)
    }
    else if (length(x)==0){
        return(sd(t))
    } 
    else{
        return(sd(t[-x]))
    }
} # end function
# feature no. 139 skewness weighted path length of graph
skewness.weighted.geodesic.graph <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)<2){
        print("WARNING! Single-noded graph. Geodesic distance cannot be computed. Returning NA")
        return(NA)
    }
    gr  <- graph
    E(gr)$weight  <- abs(E(gr)$weight)
    t <- shortest.paths(gr,weights=NULL)
    t <- t[subgraph.node.index,subgraph.node.index]
    # only upper triangle
    t <- triangle(t)
    x  <- which(t==Inf)
    # if only non-connections found
    if(length(x)==length(t)){
        return(Inf)
    }
    else if (length(x)==0){
        return(skewness(t))
    } 
    else{
        return(skewness(t[-x]))
    }
} # end function
# feature no. 140 sd weighted path length of graph
kurtosis.weighted.geodesic.graph <- function(graph,subgraph.node.index){
    if(length(subgraph.node.index)<2){
        print("WARNING! Single-noded graph. Geodesic distance cannot be computed. Returning NA")
        return(NA)
    }
    gr  <- graph
    E(gr)$weight  <- abs(E(gr)$weight)
    t <- shortest.paths(gr,weights=NULL)
    t <- t[subgraph.node.index,subgraph.node.index]
    # only upper triangle
    t <- triangle(t)
    x  <- which(t==Inf)
    # if only non-connections found
    if(length(x)==length(t)){
        return(Inf)
    }
    else if (length(x)==0){
        return(kurtosis(t))
    } 
    else{
        return(kurtosis(t[-x]))
    }
} # end function
#=============================================================================================
# feature no 141-144: kurtosis stress centrality of subgraph nodes
stress.centrality <- function(graph,subgraph.node.index){
    node.index <- c(1:vcount(graph))
    node.index <- node.index[-subgraph.node.index]
    # loop over all nodes and find geodesic distances between them
    stress.cent <- NULL  # variable storing the different stress centrality measures of the nodes in the subgraph
    for(l in subgraph.node.index){
        count.distances <- 0 # vector("list",0)
        for(i in 1:(length(node.index)-1)){
                distances <- get.all.shortest.paths(graph,node.index[i],node.index[-(1:i)],weights=NA,mode="all")$res
                # only if there is a connection between the two nodes
                if(length(distances)>0){
                    m <- lapply(distances,function(x) which(x==l))
                    count <- lapply(m,function(x) if(length(x)>0) return(1))
                    count.distances <- count.distances + sum(unlist(count))
                } # end if
        } # end for loop
        stress.cent <- c(stress.cent,count.distances)
    } # end for loop index l
    return(c(mean(stress.cent),sd(stress.cent),skewness(stress.cent),kurtosis(stress.cent)))
} # end function
# feature no 145-148: kurtosis stress centrality of subgraph nodes
weighted.stress.centrality <- function(graph,subgraph.node.index){
    options(warn = -1)
    node.index <- c(1:vcount(graph))
    node.index <- node.index[-subgraph.node.index]
    # set all weights to positive values
    E(graph)$weight  <- abs(E(graph)$weight)
    # loop over all nodes and find geodesic distances between them
    stress.cent <- NULL  # variable storing the different stress centrality measures of the nodes in the subgraph
    for(l in subgraph.node.index){
        count.distances <- 0 # vector("list",0)
        for(i in 1:(length(node.index)-1)){
            distances <- get.all.shortest.paths(graph,node.index[i],node.index[-(1:i)],weights=NA,mode="all")$res
            # only if there is a connection between the two nodes
            if(length(distances)>0){
                m <- lapply(distances,function(x) which(x==l))
                count <- lapply(m,function(x) if(length(x)>0) return(1))
                count.distances <- count.distances + sum(unlist(count))
            } # end if
        } # end for loop
        stress.cent <- c(stress.cent,count.distances)
    } # end for loop index l
    options(warn = 0)
    return(c(mean(stress.cent),sd(stress.cent),skewness(stress.cent),kurtosis(stress.cent)))
} # end function
## ------------------------------FUNCTIONS END-------------------------------------------------------------------##
