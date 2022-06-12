import networkx as nx
import gurobipy as gp
import itertools
import sys
import time

def find_srg(n, k, lmd, mu):
    """
    Tries to find a strongly regular graph with parameters (n, k, lmd, mu)
    :param n: number of nodes
    :param k: degree of each node
    :param lmd: number of vertices adjacent to a pair of adjacent vertices
    :param mu: number of vertices adjacent to a pair of non-adjacent vertices
    """
    m = gp.Model("model")
    BigM = n

    # variable for edge chosen
    x = m.addVars(
        itertools.combinations(range(n), 2),
        vtype=gp.GRB.BINARY,
        name="x"
    )
    # whether w has an edge to both u and v (helper variable)
    y = m.addVars(
        itertools.product(range(n), range(n), range(n)), # w, u, v
        vtype=gp.GRB.BINARY,
        name="y"
    )
    # degree is equal to k
    m.addConstrs(
        gp.quicksum(x[tuple(sorted([v, u]))] for u in range(n) if u != v) == k
        for v in range(n)
    )

    # y variables record the correct property
    #   require 1/2 (x_wl + x_wr) - 1/2 <= y_wuv <= 1/2 (x_wl + x_wr)
    m.addConstrs(
        -1 + x[tuple(sorted([w, u]))] + x[tuple(sorted([w, v]))] <= 2 * y[w, u, v]
        for w, u, v in itertools.product(range(n), range(n), range(n)) if w != u and w != v and u != v
    )
    m.addConstrs(
        2 * y[w, u, v] <= x[tuple(sorted([w, u]))] + x[tuple(sorted([w, v]))]
        for w, u, v in itertools.product(range(n), range(n), range(n)) if w != u and w != v and u != v
    )

    # for every pair of adjacent vertices, there are l vertices adjacent to both
    m.addConstrs(
        gp.quicksum(y[w, u, v] for w in range(n) if w != u and w != v) <= lmd + ((1 - x[u, v]) * BigM)
        for u, v in itertools.combinations(range(n), 2)
    )
    m.addConstrs(
        gp.quicksum(y[w, u, v] for w in range(n) if w != u and w != v) >= lmd - ((1 - x[u, v]) * BigM)
        for u, v in itertools.combinations(range(n), 2)
    )

    # for every pair of non-adjacent vertices, there are m vertices adjacent to both
    m.addConstrs(
        gp.quicksum(y[w, u, v] for w in range(n) if w != u and w != v) <= mu + (x[u, v] * BigM)
        for u, v in itertools.combinations(range(n), 2)
    )
    m.addConstrs(
        gp.quicksum(y[w, u, v] for w in range(n) if w != u and w != v) >= mu - (x[u, v] * BigM)
        for u, v in itertools.combinations(range(n), 2)
    )

    print("Start optimizing model")
    tstart = time.time()
    m.optimize()
    print("Finished optimizing in", time.time() - tstart, "seconds")

    # check if model has solution
    if m.getAttr("Status") != 2:
        print("Model has no solution\n")
        return

    # get edges in graph
    E = [e for e in itertools.combinations(range(n), 2) if x[e].X > 0.5]

    # check if graph strongly regular
    correct = True
    G = nx.empty_graph(n)
    G.add_edges_from(E)

    # check degree
    for v in range(n):
        correct &= G.degree(v) == k

    # check lambda and mu adjacency
    for u, v in itertools.combinations(range(n), 2):
        if x[(u, v)].X > 0.5: # adjacent
            cmp = lmd
        else: # non-adjacent
            cmp = mu

        correct &= cmp == len([w for w in range(n) if w != u and w != v and x[tuple(sorted([w, u]))].X > 0.5 and x[tuple(sorted([w, v]))].X > 0.5])

    # print solution
    print("Found", correct, "solution", E)
    print("Adjacency matrix:")
    print(nx.to_numpy_matrix(G))



def main():
    # Petersen graph (default)
    n = 10 # number of nodes
    k = 3 # degree
    lmd = 0 # lambda
    mu = 1 # mu

    # read input
    if len(sys.argv) < 2:
        print("Defaulting to n =", n, "k =", k, "lambda =", lmd, "mu =", mu)
    elif len(sys.argv) != 5:
        print("Usage: stronglyRegularGraphs.py <n k lambda mu>")
        return
    else:
        try:
            n = int(sys.argv[1])
            k = int(sys.argv[2])
            lmd = int(sys.argv[3])
            mu = int(sys.argv[4])
        except:
            print("Please specify integers n, k, lambda, and mu")
            print("Usage: stronglyRegularGraphs.py <n k lambda mu>")
            return

    find_srg(n, k, lmd, mu)



if __name__ == "__main__":
    main()