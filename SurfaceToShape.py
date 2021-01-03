from Graphs import EmbeddedGraph, EmbodiedGraph
from numpy import array, argmin, argmax, sqrt, sum, power, arctan2, pi as PI
from numpy.fft import fft

# perimeter of a planar graph
class Surface2D:
    def __init__(self, embodied_graph, obj = None):
        self.G = embodied_graph
        self._label = [v for v in range(self.G.V())]
        # sort vertices by distance to the center of mass.
        vertices_by_distance = sorted([v for v in range(self.G.V())], key=lambda v: self.G.dist(v))
        # the vertex furthest from the center must lie on the surface,
        f = vertices_by_distance[-1]
        if obj == None:
            # use a simple maximization scheme to iteratively construct our perimeter cycle.
            obj = self.G.dist
        if obj == 'min-adj':
            obj = lambda v: 1 / self.G.deg(v)
        tail = f
        head = f
        self._label[head] = self._label[tail]
        cycle = []
        while self._hom_deg(tail) < 2 and len(cycle) < self.G.V():
            unvisited = self._het_adj(head)
            if len(unvisited) == 0:
                break
            cycle.append(head)
            head = unvisited[argmax([obj(v) for v in unvisited])]
            self._label[head] = self._label[tail]
        cycle.append(head)
        # the surface is a cycle ...
        print('cyclic', cycle[0] in self._hom_adj(cycle[-1]))
        n = len(cycle)
        #assert set(cycle) == set(vertices_by_distance[-n:])

        self._surface = cycle
        self._n = n

    # heterogeneous adjacent
    def _het_adj(self, v):
        return [w for w in self.G.adj(v) if self._label[v] != self._label[w]]

    # heterogeneous degree
    def _het_deg(self, v):
        return len(self._het_adj(v))

    # homogeneous adjacent
    def _hom_adj(self, v):
        return [w for w in self.G.adj(v) if self._label[v] == self._label[w]]

    # homogeneous degree
    def _hom_deg(self, v):
        return len(self._hom_adj(v))

    def perimeter(self):
        return self._n

    def vertices(self):
        return self._surface

    def points(self):
        return [self.G.pos(v) for v in self.vertices()]

    def distance_sequence(self):
        return [self.G.dist(v) for v in self.vertices()]

    def characteristic(self):
        seq = self.distance_sequence()
        return fft(seq) / self.perimeter()

def main():
    import sys, re
    tupler = lambda edge: re.sub(r'[()]', '', edge).split(',')
    to_ints = lambda x: (int(x[0]), int(x[1]))
    parse = lambda edge: to_ints(tupler(edge))
    G = None
    ext = sys.argv[1][-4:]
    if ext == '.txt':
        with open(sys.argv[1], 'r') as graph_str:
            data = [line.replace('\n', '') for line in graph_str]
            V = int(data[0])
            D = int(data[1])
            E = int(data[2+V])
            embedding = [parse(embed) for embed in data[2:2+V]]
            edges = [parse(edge) for edge in data[3+V:3+V+E]]
            G = EmbodiedGraph(V, D, embedding, edges)
            assert G.E() == E
    elif ext == '.png':
        from skimage.io import imread
        from skimage.filters import threshold_otsu
        blobsource = imread(sys.argv[1], as_gray=True)
        D = len(blobsource.shape)
        assert D == 2
        w, h = blobsource.shape
        print(w,h)
        is_node = blobsource > threshold_otsu(blobsource)
        nbhd_rule = 4
        if nbhd_rule == 4:
            adj = lambda x, y: [(x-1,y),(x,y-1),(x+1,y),(x,y+1)]
        elif nbhd_rule == 8:
            adj = lambda x, y: [(x-1,y),(x-1,y-1),(x,y-1),(x+1,y-1),(x+1,y),(x+1,y+1),(x,y+1),(x-1,y+1)]
        inv_emb = dict()
        V = 0
        for x in range(w):
            for y in range(h):
                if is_node[x,y] == 1:
                    inv_emb[(x,y)] = V
                    V += 1
        edges = []
        embedding = list(inv_emb.keys())
        for v in range(V):
            x1,y1 = embedding[v]
            for x2,y2 in adj(x1,y1):
                if is_node[x2,y2]:
                    w = inv_emb[x2,y2]
                    if v < w:   edges.append((v,w))
        G = EmbodiedGraph(V, D, embedding, edges)
    if G:
        center = G.c()
        print('center')
        print(center, G.pos(center), G.dist(center), G.angle(center))
        surface = Surface2D(G)
        print('surface')
        print(surface.perimeter(), surface.points(), surface.distance_sequence(), surface.characteristic())

if __name__ == '__main__':
    main()
