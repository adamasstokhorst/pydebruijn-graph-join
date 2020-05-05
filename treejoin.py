from __future__ import print_function
from helpers import *

# //////////////////////////////// #
# //                            // #
# //  MODIFY THIS SECTION ONLY  // #
# //                            // #
# //////////////////////////////// #

# Provide input ANF as list of tuples, where each tuple
# describes a single term, and the list is a summation
# of all the terms -- e.g. [(), (1,), (0, 2)] -> 1 + x_1 + x_0*x_2
input_anf = [(1,), ()]
input_order = 4

# //////////////////////////////// #
# //////////////////////////////// #

def join_trees_from_anf(p_anf, order, verbose=False):
    """
    Yields `FeedbackShiftRegister`s when imported.
    Refer to `pydebruijn` module for more details.
    """
    to_console = lambda x: x
    if verbose:
        to_console = print
    
    func, order = parse_anf(p_anf, order)
    to_console('Input: ANF {} of order {}'.format(func, order))
    
    # State graph
    G = state_graph(func, order)
    cyc = cycles(G, components(G))
    e_cyc = [a * (-(-order / len(a)) + 1) for a in cyc]
    # For easy look-up on original length of cycles
    cyc_lens = {b: len(a) for a, b in itertools.izip(cyc, e_cyc)}
    to_console('{} cycles: {}'.format(len(cyc), cyc))
    
    # Tree adjacency graph (a.k.a. preference adjacency graph)
    H = nx.MultiDiGraph()
    H.add_nodes_from(e_cyc)
    
    # Populate tree adjacency graph
    for extended_cycle in e_cyc:
        for i in range(cyc_lens[extended_cycle]):
            state = extended_cycle[i:i+order]
            conj_state = state[:-1] + str(1 - int(state[-1]))
            
            # Check if conjugate has in-edge...
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

    # This bit saves an image of the preferred adjacency graph to a file -- requires pygraphviz
    # H2 = simplify(nx.relabel_nodes(H, {a: a[:cyc_lens[a]] for a in e_cyc}, copy=True))
    # H2.graph['edge'] = {'splines': 'curved'}
    # H2.graph['graph'] = {'scale': '3'}
    # nx.drawing.nx_agraph.to_agraph(H2).draw('preference_adjacency_graph.png', prog='dot')
    
    # Find spanning trees
    for tree in pydebruijn.helpers.spanning_trees(nx.Graph(H)):
        # Each edge can be in two directions
        # Each direction may have multiple directed edge
        # Then choose an initial state
        for edge_flips in itertools.product([False, True], repeat=len(tree)):
            directed_edges = [e if not f else e[::-1] for e, f in zip(tree, edge_flips)]
            multi_edges = [H.get_edge_data(*e) for e in directed_edges]
            node_source = set([e[0] for e in directed_edges])
            node_sinks = set([e[1] for e in directed_edges])
            node_root = node_sinks - node_source
            # There needs to be exactly one node which all other nodes eventually feed into
            if len(node_root) != 1:
                continue
            c = node_root.pop()
            to_console('Spanning tree: {}'.format(directed_edges))
            for roots in itertools.product(*[[v['state'] for v in d.itervalues()] for d in multi_edges]):
                for i in range(cyc_lens[c]):
                    to_console('Adjacent states: {}, initial state: {}'.format(roots, c[i:i+order]))
                    init_state = map(int, c[i:i+order])
                    visited_roots = {a: False for a in roots}

                    # This part is borrowed from GPO routine
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
                            state[-1] = 1 - state[-1]
                            if ''.join(map(str, state)) in seq:
                                state[-1] = 1 - state[-1]
                        else:
                            visited_roots[''.join(map(str, state))] = True
                        seq += str(state[-1])
                        table[''.join(map(str, state[:-1]))] = state[-1] ^ last_bit
                        if len(table) == 2 ** (order - 1):
                            break
                    # After the trees are joined, use pydebruijn module to construct a FSR
                    fsr = pydebruijn.FeedbackShiftRegister(pydebruijn.helpers.anf_from_truth_table(table, order),
                                                           order, [0] * order)
                    print('    dBSeq: {} ({})'.format(''.join(map(str, fsr.sequence())), fsr.anf))
                    print()
                    
                    yield fsr


if __name__ == '__main__':
    for _ in join_trees_from_anf(input_anf, input_order, verbose=True):
        continue
