from SurfaceToShape import EmbeddedGraph
import sys, re
from random import randint
tupler = lambda edge: re.sub(r'[()]', '', edge).split(',')
to_ints = lambda x: (int(x[0]), int(x[1]))
parse = lambda edge: to_ints(tupler(edge))
try:
    with open(sys.argv[1], 'r') as graph_str:
        data = [line.replace('\n', '') for line in graph_str]
        V = int(data[0])
        D = int(data[1])
        E = int(data[2+V])
        embedding = [parse(embed) for embed in data[2:2+V]]
        edges = [parse(edge) for edge in data[3+V:3+V+E]]
        G = EmbeddedGraph(V, D, embedding, edges)
        assert G.E() == E
except FileNotFoundError:
    print(f'the file {sys.argv[1]} could be located.')
except:
    print(f'failed to load graph from {sys.argv[1]}.')
else:
    print(
        f'loaded an embedded graph from {sys.argv[1]}.\n'
        f'V: {G.V()}\tE: {G.E()}\tD: {G.D()}')

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
points = [(x_offset + x_scale * x, y_offset + y_scale * y) for (x,y) in embedding]
print(
    f'\nrescaled graph embedding to fit screen {w,h} with {100*pad}% padding.\n'
    f'original: {embedding}\n'
    f'rescaled: {points}')

dot_size, line_size = int(sys.argv[5]), int(sys.argv[6])
vtx_color =  [pygame.Color(randint(0,200),randint(0,200),randint(0,200)) for v in range(G.V())]
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
