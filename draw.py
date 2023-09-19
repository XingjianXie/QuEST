import json
import pydot

# Load the JSON data
with open('result.json') as f:
    data = json.load(f)

# Create a directed graph
G = pydot.Dot(graph_type='digraph')

# Set the rank direction to left-to-right ('LR')
G.set_rankdir('LR')

l = set()

# Recursive function to add nodes and edges
def add_dependencies(item, dependencies):
    keep = item['type'] == 'internal'
    for dep in dependencies:
        keep_sub = dep['type'] == 'internal'
        dep_name = dep['name']
        dep_type = dep['type']
        sub_node = None
        sub_edge = None
        if not any(n.get_name() == dep_name for n in G.get_nodes()):
            color = 'red' if dep_type == 'internal' else 'black'
            if dep_type == 'internal':
                l.add(dep_name)
            sub_node = pydot.Node(dep_name, type=dep_type, color=color)
        if not G.get_edge(item['name'], dep_name):  # Only add the edge if it does not exist
            sub_edge = pydot.Edge(item['name'], dep_name)
        if 'dependencies' in dep:
            keep_sub = add_dependencies(dep, dep['dependencies']) or keep_sub
        if keep_sub:
            if sub_node:
                G.add_node(sub_node)
            if sub_edge:
                G.add_edge(sub_edge)
        keep = keep or keep_sub
    return keep

# Iterate through the data and add nodes and edges
for item in data:
    if item['type'] == 'unitary':
        G.add_node(pydot.Node(item['name'], type='unitary', color='blue'))
        if 'dependencies' in item:
            add_dependencies(item, item['dependencies'])

# Add Legend
legend = pydot.Cluster(graph_name="legend", label="Legend", fontsize="20")
legend.add_node(pydot.Node("Unitary", fontsize="10", color="blue"))
legend.add_node(pydot.Node("Common", fontsize="10", color="black"))
legend.add_node(pydot.Node("Internal", fontsize="10", color="red"))
G.add_subgraph(legend)

# Save the graph
G.write_png('result.png')

print(json.dumps(list(l), indent=2))
