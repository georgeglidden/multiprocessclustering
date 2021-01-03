from MPC import Graph
from numpy import array, argmin, argmax, sqrt, sum, power, arctan2, pi as PI
from numpy.fft import fft

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
        tail = f
        head = f
        self._label[head] = self._label[tail]
        cycle = []
        while self._hom_deg(tail) < 2 and len(cycle) < self.G.V():
            unvisited = self._het_adj(head)
            cycle.append(head)
            head = unvisited[argmax([obj(v) for v in unvisited])]
            self._label[head] = self._label[tail]
        cycle.append(head)
        # the surface is a cycle ...
        assert cycle[0] in self._hom_adj(cycle[-1])
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

class Surface2D_momentum(Surface2D):
    def __init__(self, embodied_graph, sustain = 40):
        self.G = embodied_graph
        self._label = [v for v in range(self.G.V())]
        vertices_by_distance = sorted([v for v in range(self.G.V())], key=lambda v: self.G.dist(v))
        f = vertices_by_distance[-1]
        tail = f
        head = f
        self._label[head] = self._label[tail]
        cycle = []
        while self._hom_deg(tail) < 2 and len(self._het_adj(head)) > 0:
            cycle.append(head)
            trail = cycle[-sustain:]
            t = len(trail)
            delta_angle = lambda v, w: self.G.angle(w) - self.G.angle(v)
            adjusted_delta_angle = lambda v, w: min(delta_angle(v, w), 2*PI - delta_angle(v,w))
            if t > 0:
                sum_delta_angle = 0.0
                for v in range(t-1):
                    w = v+1
                    sum_delta_angle += adjusted_delta_angle(v,w)
                avg_delta_angle = sum_delta_angle / t
                # momentum objective
                #obj = lambda v: self.G.dist(v) - (sustain * abs(avg_delta_angle - delta_angle(head, v)))
            obj = lambda v: self.G.dist(v) - (abs(avg_delta_angle - adjusted_delta_angle(head, v)))
            unvisited = self._het_adj(head)
            head = unvisited[argmax([obj(v) for v in unvisited])]
            self._label[head] = self._label[tail]
            if t > 0:
                print(f'{len(cycle)} {self.G.dist(head)}\t{abs(avg_delta_angle - adjusted_delta_angle(cycle[-1], head))}\t{self.G.dist(head) - (abs(avg_delta_angle - adjusted_delta_angle(cycle[-1], head)))}')

        cycle.append(head)
        # the surface is a cycle
        assert cycle[0] in self._hom_adj(cycle[-1])
        n = len(cycle)

        self._surface = cycle
        self._n = n

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
