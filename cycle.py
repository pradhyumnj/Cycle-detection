import networkx as nx
import json
import matplotlib.pyplot as plt

# Use of any function from networkx.algorithms module is strictly not allowed.
# Matplotlib can be used for visualization purposes

# Add your functions here if needed

global colormap
global reaction_ids

def construct_metabolic_graph(name_of_json_file):
    """
    Given the reaction dict, return the metabolic directed bi-partite networkx graph
    """
    
    global colormap # To store the color of nodes
    global reaction_ids
    G = nx.DiGraph()
    obj = open("e_coli_core.json")
    data = json.load(obj)
    reactions = data.get("reactions")
    reaction_ids = [reaction["id"] for reaction in reactions]
    for reaction in reactions:
        reaction_id = reaction['id']
        metabolites = reaction['metabolites']
        reactants = [x for x,y in metabolites.items() if y > 0]
        products = [x for x,y in metabolites.items() if y < 0]
        edges_to_reaction = ((x, reaction_id) for x in reactants) 
        edgest_from_reaction = ((reaction_id, x) for x in products)
        G.add_edges_from(edges_to_reaction)
        G.add_edges_from(edgest_from_reaction)
    colormap = ['red' if node in reaction_ids else "blue" for node in G]
    return G

def Find_Cycles_In_Metabolic_Graph(G):
    """
    Given a bi-partite networkx graph, return the list of list of metabolites involved in the cycle i.
    """
    a = []
    for x in G.nodes():
        state = {node : len(list(G.predecessors(node))) for node in G}
        paths = []
        src = x
        def DFS(node,p):
            if p[-1] == src:
                paths.append(p)
            if state[node] == 0:
                return
            state[node] -= 1
            for n in G.neighbors(node):
                DFS(n,p + [n])
        DFS(src,[src])
        paths[1:]
        a += paths
    return a

# For the purpose of parsing, look under "metabolites" sub-section of "reaction" section. The stoichiometric coefficients will help you in determining if a metabolite is a reactant or a product.
G = construct_metabolic_graph("e_coli_core.json")

# Add code to visualize G below
%matplotlib inline
fig, ax = plt.subplots(figsize = (30,30))
nx.draw(G, ax = ax, node_color = colormap, with_labels = True)
# For better visualization run the following line and open "graph.pdf"
# plt.savefig("graph.pdf")

cycles = Find_Cycles_In_Metabolic_Graph(G)
print(G.number_of_edges())
print(len(cycles))