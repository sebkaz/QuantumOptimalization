# modularity dla dowolnego grafu 

import networkx as nx 
import networkx.algorithms.community as nx_comm


# modularity funkcja 
# graf i jego podzial  

def compute_modularity(graph= None, communities=None):
    G = nx.barbell_graph(3,0)
    return nx_comm.modularity(G, [{0,4,3}, {2,1,5}])

print(compute_modularity())

def cm(graph, comms, lam):
        return nx_comm.modularity(graph, comms) + lam*len(graph.nodes)

def compute_num(lam):
    
    graph = nx.barbell_graph(3,0)
    results = []
    nody = list(graph.nodes)
    ile = len(nody)
    while ile-1:
        podzial = nx_comm.louvain_communities(graph)
        nasza_wielkosc = cm(graph, podzial, lam)
        results.append(nasza_wielkosc)
        usun = -1*ile
        graph.remove_node(nody[usun])
        ile -= 1

    return results

