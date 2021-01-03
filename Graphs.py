from numpy import array, argmin, argmax, sqrt, sum, power, arctan2, pi as PI
# just a graph
class Graph:

    def __init__(self, V, edges):
        self._V = V
        self._adj = [[] for _ in range(self.V())]
        for i,j in edges:
            self._adj[i].append(j)
            self._adj[j].append(i)

    def V(self):
        return self._V

    def adj(self, i):
        return self._adj[i]

    def deg(self, i):
        return len(self.adj(i))

class CC:
    def __init__(self, G):
        self.G = G
        self.id = [-1 for v in range(self.G.V())]
        self.untagged = lambda v: self.id[v] == -1
        self.count = 0
        for s in range(self.G.V()):
            if self.untagged(s):
                self.dfs(s)
                self.count += 1

    def dfs(self, v):
        self.id[v] = self.count
        for w in self.G.adj(v):
            if self.untagged(w):
                self.dfs(w)

    def get_components(self):
        k = self.count
        components = [[] for i in range(k)]
        for v in range(self.G.V()):
            c = self.id[v]
            components[c].append(v)
        assert sum([len(component) for component in components]) == self.G.V()
        return components

# weighted graph
class WGraph(Graph):

    def __init__(self, V, edges):
        self._V = V
        self._adj = [[] for _ in range(self.V())]
        for i,j,w in edges:
            self._adj[i].append((j,w))
            self._adj[j].append((i,w))

# symmetric distance between x1 and x2
def euclidean(x1, x2):
    return sqrt(sum(power(array(x1)-array(x2),2), axis=0))

# angle from x1 to x2
def arctangent(x1, x2):
    return arctan2(x1[1]-x2[1], x1[0]-x2[0])

class EmbeddedGraph(Graph):

    # V : int
    # embedding : list of float D-tuples
    # edges : list of integer 2-tuples
    # dist : embedding^{2} -> reals
    def __init__(self, V, D, embedding, edges,
            dist=euclidean, angle=arctangent):
        # all nodes are embedded
        assert len(embedding) == V
        # consistent dimension
        assert [len(e) for e in embedding] == V * [D]
        super().__init__(V, edges)
        self._D = D
        self._embedding = embedding
        self._dist = dist
        self._angle = arctangent

    def D(self):
        return self._D

    def E(self):
        return sum([len(self._adj[v]) for v in range(self.V())]) // 2

    def pos(self, v):
        return self._embedding[v]

    def dist(self, v, w):
        return self._dist(self.pos(v), self.pos(w))

    def angle(self, v, w):
        return self._angle(self.pos(v), self.pos(w))


class EmbodiedGraph(EmbeddedGraph):

    def __init__(self, V, D, embedding, edges,
            mass=lambda v: 1,dist=euclidean,angle=arctangent):
        super().__init__(V, D, embedding, edges, dist=dist,angle=angle)
        self._mass = mass

        esum = [0.0 for d in range(self.D())]
        wsum = 0.0
        for v in range(self.V()):
            w = self.mass(v)
            wsum += w
            e = self.pos(v)
            for d in range(self.D()):
                esum[d] += e[d]
        self._pos_c = tuple([esum[d]/wsum for d in range(self.D())])
        self._c = argmin([dist(self.pos(v), self._pos_c) for v in range(self.V())])

    def mass(self, v):
        return self._mass(v)

    def pos_c(self):
        return self._pos_c

    def vtx_c(self):
        return self._c

    def dist(self, v, w=None):
        if w == None:
            return self._dist(self.pos(v), self.pos_c())
        else:
            return super().dist(v, w)

    def angle(self, v, w=None):
        if w == None:
            return self._angle(self.pos(v), self.pos_c())
        else:
            return super().angle(v, w)

# weighted, labeled graph
class WLGraph(WGraph):

    def __init__(self, edges, labels):
        super().__init__(len(labels), edges)
        self._labels = labels

    def get_label(self, i):
        return self._labels[i]

    def set_label(self, i, l):
        self._labels[i] = l


# constrained, labeled graph
class CLGraph(WLGraph):

    def __init__(self, edges, labels, constraint_fn):
        super().__init__(edges, labels)
        self.c = constraint_fn

    def adj(self, i, p):
        return [(j,w) for (j,w) in self._adj[i] if self.c(w,p)]

    # O(E + V)
    def apply_constraint(weighted_graph, constraint_fn):
        edges = []
        for i in range(weighted_graph.V()):
            # ensure i < j to avoid double counting edges
            adj = [(j,w) for (j,w) in weighted_graph.adj(i) if j > i]
            for j,w in adj:
                edges.append((i,j,w))
        labels = [i for i in range(weighted_graph.V())]
        return CLGraph(edges, labels, constraint_fn)
