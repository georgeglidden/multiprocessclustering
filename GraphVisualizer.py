from SurfaceToShape import EmbodiedGraph, Surface2D
import sys, re
from random import randint
tupler = lambda edge: re.sub(r'[()]', '', edge).split(',')
to_ints = lambda x: (int(x[0]), int(x[1]))
parse = lambda edge: to_ints(tupler(edge))
try:
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
        nbhd_rule = int(input('4-adj or 8-adj? (enter an integer): '))
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
except FileNotFoundError:
    print(f'the file {sys.argv[1]} could be located.')

#except Exception as ex:
#    print(f'failed to load graph from {sys.argv[1]}.')
#    print(ex)
else:
    print(
        f'loaded an embedded graph from {sys.argv[1]}.\n'
        f'V: {G.V()}\tE: {G.E()}\tD: {G.D()}')

    surface = Surface2D(G, obj='min-adj')

    print('\ninitializing ux...')
    import pygame
    pygame.init()
    w,h = int(sys.argv[2]), int(sys.argv[3])
    screen = pygame.display.set_mode((w,h))

    pad = float(sys.argv[4])
    x_min, x_max = min(embedding, key=lambda e:e[0]), max(embedding, key=lambda e:e[0])
    G_w = x_max[0] - x_min[0]
    x_offset = w * pad
    x_scale = (w - (2*x_offset)) / G_w
    y_min, y_max = min(embedding, key=lambda e:e[1]), max(embedding, key=lambda e:e[1])
    G_h = y_max[1] - y_min[1]
    y_offset = h * pad
    y_scale = (h - (2*y_offset)) / G_h
    scale_embedding = lambda x, y : (x_offset + x_scale * (x - x_min[0]), y_offset + y_scale * (y - y_min[1]))
    points = [scale_embedding(x,y) for (x,y) in embedding]
    center_point = scale_embedding(G.pos_c()[0], G.pos_c()[1])
    print(
        f'\nrescaled graph embedding to fit screen {w,h} with {100*pad}% padding.\n'
        f'G_w: {G_w} G_h: {G_h} x_s: {x_scale} y_s: {y_scale} y_o: {y_offset} x_o: {x_offset} \n'
        f'original: {embedding[:min(G.V(), 20)]}\n'
        f'rescaled: {points[:min(G.V(), 20)]}')

    dot_size, line_size = int(sys.argv[5]), int(sys.argv[6])
    vtx_color =  [pygame.Color(100,50,randint(0,50)) for v in range(G.V())]
    for i in range(surface.perimeter()):
        v = surface.vertices()[i]
        vtx_color[v] = pygame.Color(0, 100+int(155 * i/surface.perimeter()), 100+int(155 * (surface.perimeter() - i)/surface.perimeter()))
    center_color = pygame.Color(255, 205, 30)
    edge_color = pygame.Color(30,30,30)
    bg_color = pygame.Color(255,255,255)
    while True:
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                print('stopping ux...')
                pygame.quit()
                print('goodbye')
                sys.exit()
        screen.fill(bg_color)
        pygame.draw.circle(screen, center_color, center_point, 2*dot_size)
        for v in range(G.V()):
            x1,y1 = points[v]
            for w in G.adj(v):
                if v > w:
                    x2,y2 = points[w]
                    pygame.draw.line(screen, edge_color, (x1,y1), (x2,y2), line_size)
        for v in range(G.V()):
            x,y = points[v]
            color = vtx_color[v]
            pygame.draw.circle(screen, color, (x,y), dot_size)
        pygame.display.update()
