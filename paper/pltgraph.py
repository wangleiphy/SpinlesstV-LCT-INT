#import matplotlib 
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from xml.etree import ElementTree
import sys

def parse(fn):
    root = ElementTree.parse(fn).getroot()

    vertices = {}
    for vertex in root.findall('VERTEX'):
        vid = vertex.get('id')
        vpos = tuple([float(x) for x in vertex.find('COORDINATE').text.split()])
        vertices[vid] = dict( {'pos':vpos,'type':vertex.get('type') }  )
    
    edges = {}
    for edge in root.findall('EDGE'):
        edges[edge.get('id')] = dict([(k,edge.get(k)) for k in ['source', 'target', 'type', 'vector']])

   
    #print 'vertices', vertices
    #print 'edges', edges
    return (vertices,edges)

def showgraph(graph):
    vertices = graph[0]
    edges = graph[1]
    
    #bond type cmap 
    bond_map = {'0':'k','1':'k','2':'b','3':'r','4':'c','5':'m','6':'y','7':'k'} 
    #site type cmap 
    site_map = {'0':'k','1':'k','2':'b','3':'r','4':'c','5':'m','6':'y','7':'k'} 

    symbol_map = {'0':'k', '1':'w'}

    for edge in edges.values():
        s = edge['source']
        t = edge['target']
        c = edge['type']

        p0 = vertices[s]['pos']
        p1 = vertices[t]['pos']
       
        #draw bond type 
        plt.plot([p0[0], p1[0]], [p0[1], p1[1]], c=bond_map[c],lw=1,zorder=1)
    
        #u = p1[0] - p0[0] 
        #v = p1[1] - p0[1] 
        #if u>=0 and v>=0:
        #    print p0, p1
        #    plt.quiver(p0[0], p0[1], u, v, color=bond_map[c], zorder=1, scale_units='xy', angles='xy', scale=1)
        

    x = [v['pos'][0] for v in vertices.values()]
    y = [v['pos'][1] for v in vertices.values()]

    c = [site_map[v['type']] for v in vertices.values()]
    fc = [symbol_map[v['type']] for v in vertices.values()]

    plt.scatter(x, y, s=120, zorder=10, facecolors=fc, edgecolors=c, linewidth=2)

    #for k, v in vertices.items():
    #    print k, v['pos']
    #    plt.annotate('%s' %k, v['pos'])

if __name__ == '__main__':
    import matplotlib.pyplot as plt 
    import sys 
    plt.figure(figsize=(6,4))
    plt.box('off')

    graph = parse(sys.argv[1])
    
    showgraph(graph)
    #plt.show()
    plt.savefig('lattice.pdf', transparent=True)
    sys.exit(1)

    #read lattice 
    Row = []
    Col = []
    Val = []
    edges = graph[1]
    for edge in edges.values():
        Row.append(int(edge['source'])-1)
        Col.append(int(edge['target'])-1)
        Val.append(int(edge['type']))

    import scipy.sparse as sps 
    Nsite = max(max(Row), max(Col)) + 1
    Kmat = sps.csr_matrix((Val, (Row, Col)), shape=(Nsite, Nsite))
    
    Ham = Kmat.todense()
    Ham += Ham.transpose()
    import matplotlib.pyplot as plt 
    plt.spy(Ham , markersize=10 ,marker='.')
    plt.show()
