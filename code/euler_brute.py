import networkx as nx
import matplotlib.pyplot as plt
import itertools
from scipy.special import comb


def create_pulsar_graph(n1, m, n2):
    vertex_count = n1 + n1 + 1 + m + 1 + n2 + n2

    G = nx.Graph()

    G.add_nodes_from(range(vertex_count))

    # build the graph top to bottom
    a1 = n1 + n1
    a2 = n1 + n1 + 1 + m
    
    for i in range(n1):
        cip = i
        ci = cip + n1
        G.add_edge(cip, ci)
        #G.add_edge(ci, cip)
        
    for i in range(n1):
        ci = i + n1
        
        G.add_edge(ci, a1)
        #G.add_edge(a1, ci)
        
    for i in range(m):
        bi = n1 + n1 + 1 + i
        G.add_edge(a1, bi)
        #G.add_edge(bi, a1)

    for i in range(m):
        bi = n1 + n1 + 1 + i
        G.add_edge(bi, a2)
        #G.add_edge(a2, bi)

    for i in range(n2):
        di = n1 + n1 + 1 + m + 1 + i
        G.add_edge(a2, di)
        #G.add_edge(di, a2)
        
    for i in range(n2):
        di = n1 + n1 + 1 + m + 1 + i
        dip = di + n2
        G.add_edge(di, dip)
        #G.add_edge(dip, di)

    return G


def are_incident(edge1, edge2):
    return edge1[0] == edge2[0] \
        or edge1[0] == edge2[1] \
        or edge1[1] == edge2[0] \
        or edge1[1] == edge2[1]

def matchings(G, k):
    if k == 0:
        yield ()
        return

    def is_disjoint(edges):
        pairs = itertools.combinations(edges, 2)
        return not any((are_incident(edge1, edge2)) for edge1, edge2 in pairs)
    
    kcombs = itertools.combinations(G.edges(), k)
    
    for kcomb in kcombs:
        if is_disjoint(kcomb):
            yield kcomb


def count_kcubes(G, n, k):
    assert(k <= n)
    vertex_count = len(G.nodes)
    return len(list(matchings(G,k))) * comb(vertex_count - 2 * k, n - k)


def euler(G, n):
    s = 0
    for k in range(n+1):
        s += ((-1)**k) * count_kcubes(G, n, k)
    return s

n1 = 5
m = 3
n2 = 4
G = create_pulsar_graph(n1, m, n2)
#nx.draw(G, with_labels=True, font_weight='bold')
#plt.show()
print(len(list(matchings(G, 2))))
print(f"euler characteristic of UConf_3(pulsar graph({n1}, {m}, {n2}): {euler(G, 3)}")
print(f"Theta_m formula: {m * (m - 2) * (m - 7) / 6}")
#print(f"3-matchings {len(list(matchings(G,3)))}")
#print(f"2-matchings {len(list(matchings(G,2)))}")
#print(f"1-matchings {len(list(matchings(G,1)))}")

formula = 2 * n1 ** 3 + (m - 6) * n1 ** 2 + (m + 3) * n1


def create_complete_graph(n):
    G = nx.Graph()
    G.add_nodes_from(range(n))
    
    pairs = itertools.combinations(G.nodes, 2)
    for v1, v2 in pairs:
        G.add_edge(v1, v2)
    return G

#G = create_complete_graph(5)
#nx.draw(G, with_labels=True, font_weight='bold')
#plt.show()
def create_complete_bipartite_graph(m, n):
    G = nx.Graph()
    G.add_nodes_from(range(m + n))
    
    for v1 in range(m):
        for v2 in range(m, m+n, 1):
            G.add_edge(v1, v2)

    return G

print(f"Euler characteristic of UConf_2(K5): {euler(create_complete_graph(5), 2)}")
print(f"Euler characteristic of UConf_2(K3,3): {euler(create_complete_bipartite_graph(3, 3), 2)}")
print(f"Euler characteristic of UConf_3(K5): {euler(create_complete_graph(5), 3)}")
print(f"Euler characteristic of UConf_3(Theta_4): {euler(create_pulsar_graph(0, 4, 0), 3)}")
print(f"Euler characteristic of UConf_4(K3,3): {euler(create_complete_bipartite_graph(3, 3), 4)}")
