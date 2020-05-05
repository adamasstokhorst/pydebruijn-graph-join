import sympy
import networkx as nx
import itertools
import pydebruijn
from pydebruijn.helpers import powerset


def parse_anf(s, degree=None):
    """
    Converts a list of tuples into an ANF.
    
    [(0, 2), (1,), ()] -> 1 + x_1 + x_0*x_2
    """
    # Implicitly derive the (minimum) degree of the ANF, if not supplied
    if degree is None:
        degree = max(map(max, s)) + 1
    syms = sympy.symbols('x_:{}'.format(degree), integer=True)
    
    # Form the ANF and return it with the degree
    anf = sympy.Integer(0)
    for term in s:
        anf += reduce(lambda a, b: a*b, map(syms.__getitem__, term), sympy.Integer(1))
    return anf, degree


def state_graph(anf, degree):
    """
    Constructs the (directed) state graph of a given ANF.
    """
    # Convert ANF to a string, for title
    title = str(anf).replace('*', '').replace('_', '')
    # terms = str(anf).replace('*', '').split('+')
    # grouped = list(zip(*[iter(terms)]*7)) + ([terms[-(len(terms)%7):]] if len(terms)%7 != 0 else [])
    # title = '${}$'.format('$\n+$'.join(map(lambda a: '+'.join(a), grouped)))
    
    # Construct the digraph
    syms = sympy.symbols('x_:{}'.format(degree), integer=True)
    graph = nx.DiGraph()
    # This map call creates padded binary strings (e.g. 00010) depending on the needed length
    graph.add_nodes_from(map(lambda n: '{:>0{width}}'.format(bin(n)[2:], width=degree), xrange(2**degree)))
    for state in map(lambda n: '{:>0{width}}'.format(bin(n)[2:], width=degree), xrange(2**degree)):
        args = zip(syms, map(int, state))
        next_bit = anf.subs(args) % 2
        next_state = state[1:] + str(next_bit)
        graph.add_edge(state, next_state)
    # Set some properties for drawing the graph
    graph.graph['edge'] = {'splines': 'curved'}
    graph.graph['graph'] = {'label': title, 'scale': '3'}
    return graph


def components(graph):
    """
    Determines the number of components in a simple graph, via traversal.
    
    Returns a list of elements (chosen arbitrarily) belonging to each component.
    """
    graph = nx.Graph(graph)
    visited = []
    seeds = []
    # Iterate through all the nodes, skipping the ones that has been
    # visited and traversing nodes that haven't been visited.
    for e in graph.nodes():
        if e in visited:
            continue
        seeds.append(e)
        stack = [e]
        # Keep adding a node's neighbors to the stack to traverse,
        # dropping it if it has been visited before.
        while stack:
            cur_node = stack.pop()
            visited.append(cur_node)
            neighbors = list(graph[cur_node])
            stack.extend([x for x in neighbors if x not in visited])
    return seeds


def cycles(graph, seeds):
    """
    Returns a list of eventual cycles in a digraph.
    
    This function is intended for state graphs only.
    """
    out = []
    for e in seeds:
        visited = []
        while e not in visited:
            visited.append(e)
            e = list(graph[e])[0]  # Every vertex only has one out-edge
        pos = visited.index(e)
        cycle = visited[pos:]

        try:
            # This only works if nodes are binary strings
            state = list(cycle[0])
            for s in cycle[1:]:
                if len(state) >= len(cycle):
                    break
                state += s[-1]
            out.append(''.join(state[:len(cycle)]))
        except TypeError:
            out.append(len(cycle))
    return out


def strip_self_loops(graph):
    """
    Returns a copy of a graph without self-loops
    """
    graph_copy = graph.copy()
    for u, v in graph.edges():
        if u == v:
            graph_copy.remove_edge(u, v)
    return graph_copy


def simplify(graph):
    """
    Simplifies parallel edges into a single edge, labelling it with
    its multiplicity.
    """
    import collections
    
    graph_copy = graph.__class__()
    graph_copy.add_nodes_from(graph.nodes)
    for u, v, _ in graph.edges:
        if not graph_copy.has_edge(u, v):
            graph_copy.add_edge(u, v, label=0)
        # This probably doesn't work if `graph` isn't a multigraph
        graph_copy.edges[u, v, 0]['label'] += 1
    for u, v, _ in graph_copy.edges:
        if graph_copy.edges[u, v, 0]['label'] == 1:
            graph_copy.edges[u, v, 0]['label'] = ''
    return graph_copy


def eval_anf(expression, degree, state):
    """
    Evaluates a given ANF at the specified values.
    """
    args = zip(sympy.symbols('x_:{}'.format(degree), integer=True), state)
    return expression.subs(args) % 2


if __name__ == '__main__':
    p_anf = [(1,)]
    order = 7
    
    func, order = parse_anf(p_anf, order)
    print 'input: ANF {} of order {}'.format(func, order)
    
    G = state_graph(func, order)
    cyc = cycles(G, components(G))
    e_cyc = [a * (-(-order / len(a)) + 1) for a in cyc]
    # for easy look-up on original length of cycles
    cyc_lens = {b: len(a) for a, b in itertools.izip(cyc, e_cyc)}
    print '{} cycles: {}'.format(len(cyc), cyc)
    
    # tree adjacency graph
    H = nx.MultiDiGraph()
    H.add_nodes_from(e_cyc)
    
    # populate tree adjacency graph
    for extended_cycle in e_cyc:
        for i in range(cyc_lens[extended_cycle]):
            state = extended_cycle[i:i+order]
            conj_state = state[:-1] + str(1 - int(state[-1]))
            
            # check if conjugate has in-edge
            if G.in_edges(nbunch=[conj_state]):
                continue
            
            # or if it feeds back to current cycle
            cs = map(int, conj_state)
            while True:
                cs = cs[1:] + [eval_anf(func, order, cs)]
                cs_str = ''.join(map(str, cs))
                if cs_str in extended_cycle:
                    break
                elif any([cs_str in a for a in e_cyc]):
                    H.add_edge(extended_cycle, [a for a in e_cyc if cs_str in a][0], state=state)
                    break

    H2 = simplify(nx.relabel_nodes(H, {a: a[:cyc_lens[a]] for a in e_cyc}, copy=True))
    H2.graph['edge'] = {'splines': 'curved'}
    H2.graph['graph'] = {'scale': '3'}
    nx.drawing.nx_agraph.to_agraph(H2).draw('preference_adjacency_graph.png', prog='dot')
    
    import pprint
    H3 = nx.relabel_nodes(H, {a: a[:cyc_lens[a]] for a in e_cyc}, copy=True)
    pcp_table = {}
    for v in H3.nodes:
        for e in H3.out_edges(v, data=True):
            if (e[0], e[1]) not in pcp_table:
                pcp_table[(e[0], e[1])] = []
            pcp_table[(e[0], e[1])].append(e[2]['state'])
    print '\nPreference companion pairs:'
    pprint.pprint(pcp_table)
    
    import sys
    sys.exit(0)
    
    # find spanning trees
    for tree in pydebruijn.helpers.spanning_trees(nx.Graph(H)):
        # each edge can be in two directions
        # each direction may have multiple directed edge
        # then choose an initial state...
        for edge_flips in itertools.product([False, True], repeat=len(tree)):
            directed_edges = [e if not f else e[::-1] for e, f in zip(tree, edge_flips)]
            multi_edges = [H.get_edge_data(*e) for e in directed_edges]
            node_source = set([e[0] for e in directed_edges])
            node_sinks = set([e[1] for e in directed_edges])
            node_root = node_sinks - node_source
            if len(node_root) != 1:
                continue
            c = node_root.pop()
            print 'spanning tree: {}'.format(directed_edges)
            for roots in itertools.product(*[[v['state'] for v in d.itervalues()] for d in multi_edges]):
                for i in range(cyc_lens[c]):
                    print 'adjacent states: {}, initial state: {}'.format(roots, c[i:i+order])
                    init_state = map(int, c[i:i+order])
                    visited_roots = {a: False for a in roots}

                    # this part is borrowed from GPO routine
                    syms = sympy.symbols('x_:{}'.format(order), integer=True)
                    table = {'0' * (order - 1): 1, '1' * (order - 1): 1}
                    state = init_state[:]
                    seq = ''.join(map(str, state))

                    while True:
                        args = zip(syms, state)
                        next_bit = int(func.subs(args)) % 2
                        last_bit = state[0]
                        state = state[1:] + [next_bit]
                        if ''.join(map(str, state)) not in roots or visited_roots[''.join(map(str, state))]:
                            if ''.join(map(str, state)) not in seq:
                                state[-1] = 1 - state[-1]
                        else:
                            visited_roots[''.join(map(str, state))] = True
                        # print '-> {}'.format(''.join(map(str, state))),
                        seq += str(state[-1])
                        table[''.join(map(str, state[:-1]))] = state[-1] ^ last_bit
                        if len(table) == 2 ** (order - 1):
                            break
                    fsr = pydebruijn.FeedbackShiftRegister(pydebruijn.helpers.anf_from_truth_table(table, order),
                                                           order, [0] * order)
                    # print
                    print '    dbseq: {} ({})'.format(''.join(map(str, fsr.sequence())), fsr.anf)
                    print
