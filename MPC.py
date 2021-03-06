from Graphs import Graph, CC, CLGraph

# union all nodes and relabel accordingly
def default_sum(process_a, process_b):
    # sanity check
    assert process_a._graph == process_b._graph
    G = process_a._graph
    # minimize the number of mislabled nodes
    if process_b.n() < process_a.n():
        process_i = process_b
        process_j = process_a
    else:
        process_i = process_a
        process_j = process_b
    print('i nodes', process_i.nodes)
    print('j nodes', process_j.nodes)
    # choose an ID and rename mislabeled nodes
    sum_id = process_i.process_id
    for v in process_j.nodes:
        G.set_label(v, sum_id)
    print('relabeled:', G._labels)
    # collect nodes
    sum_nodes = process_i.nodes.extend(process_j.nodes)
    sum_active = process_i.active.extend(process_j.active)
    print('Enodes', sum_nodes)
    print('Eactive', sum_active)
    # return the summative process
    return BFSProcess(sum_id, G, nodes=sum_nodes, active=sum_active)

class BFSProcess:

    def __init__(self, seed, constrained_graph_with_labels,
            sum_op = default_sum, nodes = None, active = None):
        self.process_id = -seed
        if nodes == None:
            self.nodes = []
        else:
            self.nodes = nodes
        self._graph = constrained_graph_with_labels
        if active == None:
            self.active = [seed]
            self._graph.set_label(seed, self.process_id)
        else:
            self.active = active
        self.sum = sum_op
        print(f'created process {self.process_id}\n\tnodes: {self.nodes}\n\tactive: {self.active}')

    def __add__(self, other):
        return self.sum(self, other)

    def n(self):
        return len(self.nodes)

    # perform a BFS
    def tick(self):
        signal = []
        new_nodes = []
        print(20*'/')
        print(f'pretick on process {self.process_id}\n\tnodes: {self.nodes}\n\tactive: {self.active}')
        for i in self.active:
            for j,w in self._graph.adj(i, self):
                label = self._graph.get_label(j)
                # ignore nodes already in the Process
                if label == self.process_id:
                    continue
                # signal to merge the intersecting processes
                elif label != j:
                    signal.append((self.process_id, label))
                    print(f'process {self.process_id} intersected process {label}')
                    continue
                # otherwise, subsume the node
                else:
                    new_nodes.append(j)
                    self._graph.set_label(j, self.process_id)
                    print(f'process {self.process_id} subsuming node {label}')
        self.nodes.extend(self.active)
        self.active = new_nodes
        print(f'posttick on process {self.process_id}\n\tnodes: {self.nodes}\n\tactive: {self.active}\n\t new: {new_nodes}')
        print(20*'\\')
        return signal

    def is_complete(self):
        if len(self.active) == 0:
            return True
        else:
            return False

    def summary(self):
        return (self.nodes, self.active, self.is_complete())

def MultiProcessCluster(seeds, N, edges, constraint_fn,
        sum_op = default_sum):
    # initialize the constrained graph with labels
    clgraph = CLGraph(edges, [n for n in range(N)], constraint_fn)
    # initialize the processes
    processes = [BFSProcess(seed, clgraph, sum_op = sum_op) for seed in seeds]
    print('init processes', [P.summary() for P in processes])
    # iterate processes to completion
    while (False in [process.is_complete() for process in processes]):
        print(50*'-')
        # tick processes and accumulate intersection signals
        signals = []
        for P in processes:
            signals.extend(P.tick())
        print('mid processes', [P.summary() for P in processes])
        print('labels', clgraph._labels)
        print('signals', signals)
        # group processes by intersection
        merge_graph = Graph(len(processes), signals)
        intersections = CC(merge_graph).get_components()
        print('intersections', intersections)
        # merge intersecting processes
        new_processes = []
        for component in intersections:
            sum_process = processes[component[0]]
            for p in component[1:]:
                sum_process += processes[p]
            new_processes.append(sum_process)
        processes = new_processes
        print('new processes', [P.summary() for P in processes])
    # extract clusters from the final process nodes
    print('all processes complete')
    clusters = [P.nodes for P in processes]
    return clusters

def main():
    V = 10
    #G = Graph(V, [(2,3),(3,1),(9,1),(4,9),(5,7),(7,8)])
    #cc = CC(G)
    #print(cc.get_components())
#
    #WG = WGraph(V, [(2,3,-1),(3,1,1),(9,1,-1),(4,9,1),(5,7,-1),(7,8,1)])
    #print(WG.adj(7))
#
    #CLG = CLGraph.apply_constraint(WG, lambda w,p: w > 0)
    #print(CLG.adj(7, 1))
    #print(CLG.adj(7, -1))
    #process = BFSProcess(7, CLG)
    #print(CLG._labels)
    #print(process.active)
    #while not process.is_complete():
    #    process.tick()
    #    print(CLG._labels)
    #    print(process.active)

    edges = [(2,3,-1),(3,1,1),(9,1,-1),(4,9,1),(5,7,-1),(7,8,1)]
    seeds = [1, 5, 6]
    constraint_fn = lambda w,p: True
    clusters = MultiProcessCluster(seeds, V, edges, constraint_fn)
    print(clusters)


if __name__ == '__main__':
    main()
