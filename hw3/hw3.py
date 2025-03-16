import re
from collections import defaultdict
import matplotlib.pyplot as plt # type: ignore

def hamming_distance(p, q):  # BA1G
    mismatch = 0
    for i in range(0, len(p)):
        if p[i] != q[i]:
            mismatch += 1

    return mismatch

def composition(k, text):  # BA3A kmer composition of a string
    kmers = []
    for i in range(0, len(text) - k + 1):
        k_mer = text[i:i + k]
        kmers.append(k_mer)

    return kmers


def reconstruct(k_mers):  # BA3B Reconstruct a String from its Genome Path
    new_string = k_mers[0]
    for kmer in k_mers[1:]:
        new_string += kmer[-1]

    return new_string


def formatting(nodes):  # printing nodes and edges
    for key in nodes:
        if len(nodes[key]) != 0:
            print(key, "->", ",".join(nodes[key]))


def overlap(kmers):  # BA3C overlap
    nodes = {}
    k = len(kmers[0])
    for kmer in kmers:
        nodes[kmer] = []
        for k_mer in kmers:
            if kmer != k_mer:
                if kmer[1:] == k_mer[:k - 1]:
                    nodes[kmer].append(k_mer)

    formatting(nodes)
    # return nodes


def de_bruijn_graph_from_string(k, text):  # BA3D
    nodes = {}

    for i in range(len(text) - k + 1):
        k_mer = text[i:i + k]
        prefix = k_mer[:-1]
        suffix = k_mer[1:]

        if prefix not in nodes:
            nodes[prefix] = []

        if suffix not in nodes[prefix]:
            nodes[prefix].append(suffix)

    # formatting(nodes)
    return nodes


def de_bruijn_graph_from_kmers(kmers):
    k = len(kmers[0]) - 1
    nodes = {}
    for kmer in kmers:
        for i in range(0, len(kmer) - k + 1):
            k_mer = kmer[i:i + k]
            nodes[k_mer] = []

    for pattern in kmers:
        prefix = pattern[0:len(pattern) - 1]
        suffix = pattern[1:len(pattern)]

        if prefix in nodes:
            nodes[prefix].append(suffix)

    # formatting(nodes)
    return nodes

def eulerian_cycle(graph):  # traverses every edge once and returns back
    start_node = list(graph.keys())[0]
    for node in graph:
        if len(graph.get(node)) != 0:  # non-zero degree
            start_node = node
            break

    s = []  # stack
    l = []  # to store cycles

    s.append(start_node)
    while len(s) != 0:
        u = s[-1]  # top of stack
        if len(graph.get(u)) != 0:  # unused edges going out
            w = graph[u].pop()  # remove edge (u,w) from the graph
            s.append(w)  # push w onto s

        else:
            u = s.pop()
            l.append(u)

    l = l[::-1]
    # return "->".join(l)
    return l


def is_unbalanced(graph):
    edges_dict = defaultdict(lambda: [0, 0])  # 0-th idx outgoing 1-th idx incoming
    unbalanced_nodes = []
    for key, value in graph.items():
        edges_dict[key][0] = len(graph[key])
        for node in value:
            edges_dict[node][1] += 1

    for k in edges_dict:
        if edges_dict[k][0] != edges_dict[k][1]:
            unbalanced_nodes.append(k)

    if len(unbalanced_nodes) == 0:
        return False
    else:
        return True


def eulerian_path(graph):
    # print(graph)
    edges_dict = defaultdict(lambda: [0, 0])  # 0-th idx outgoing 1-th idx incoming
    unbalanced_nodes = []
    for key, value in graph.items():
        edges_dict[key][0] = len(graph[key])
        for node in value:
            edges_dict[node][1] += 1

    v1, v2 = None, None
    for n in edges_dict:
        if edges_dict[n][0] > edges_dict[n][1]:  # out degrees greater
            v1 = n
        elif edges_dict[n][1] > edges_dict[n][0]:  # in degrees greater
            v2 = n

    if v2 not in graph:
        graph[v2] = []
    graph[v2].append(v1)

    cycle = eulerian_cycle(graph)

    i = 0
    for i in range(0, len(cycle) - 2):
        if cycle[i] == v2 and cycle[i + 1] == v1:
            break

    cycle = cycle[i + 1:-1] + cycle[0:i + 1]
    return cycle
    # return "->".join(cycle)


def string_reconstruction(patterns):
    graph = de_bruijn_graph_from_kmers(patterns)
    path = eulerian_path(graph)
    return reconstruct(path)


def de_bruijn_from_pairs(pairs):
    graph = {}
    for pair in pairs:  # iterate through each pair and split
        new_pair = pair.split("|")
        left, right = new_pair[0], new_pair[1]
        prefix = left[:-1] + "-" + right[:-1]
        suffix = left[1:] + "-" + right[1:]
        if prefix not in graph:
            graph[prefix] = []

        graph[prefix].append(suffix)

    return graph


def reconstruct_from_pairs(graph, k, d):
    cycle = eulerian_path(graph)
    left_parts = []
    right_parts = []
    for node in cycle:
        temp = node.split("-")
        left_parts.append(temp[0])
        right_parts.append(temp[1])

    string_1 = reconstruct(left_parts)
    string_2 = reconstruct(right_parts)
    return string_1 + string_2[-(k+d):]


def k_d_mer_composition(genome, k, d):
    read_pairs = []
    for i in range(0, len(genome)-((2*k)+d)+1):
        kmer_1 = genome[i:i+k]
        j = i + k + d
        kmer_2 = genome[j:j+k]
        read_pair = kmer_1 + "|" + kmer_2
        read_pairs.append(read_pair)
    
    return read_pairs



def find_best_k_d_reconstruction(genome, max_k, max_d):
    best_k = None
    best_d = None
    best_score = float('inf') 
    best_d_for_k = {} 

   
    for k in range(3, max_k + 1): 
        best_d_for_current_k = None
        best_score_for_current_k = float('inf')  
        for d in range(3, max_d + 1): 
            print(f"Evaluating k = {k}, d = {d}...")
            
            
            read_pairs = k_d_mer_composition(genome, k, d)
            
           
            graph = de_bruijn_from_pairs(read_pairs)
            
            
            reconstructed_genome = reconstruct_from_pairs(graph, k, d)
            
            
            score = hamming_distance(genome, reconstructed_genome)
            
            
            if score < best_score_for_current_k:
                best_score_for_current_k = score
                best_d_for_current_k = d
                
       
        best_d_for_k[k] = best_d_for_current_k
        
        
        if best_score_for_current_k < best_score:
            best_score = best_score_for_current_k
            best_k = k
            best_d = best_d_for_current_k

    return best_k, best_d, best_score, best_d_for_k



if __name__ == "__main__":
    # Load the genome
    with open("genome.txt", "r") as file:
        genome = file.read()
    genome = "".join(genome.split())  # Cleaning string

    # Reconstruct using De Bruijn Graph
    kmers = composition(4, genome) 
    gen1 = string_reconstruction(kmers)  # Genome reconstruction using De_Bruijn Graph
   
    
    # Reconstruct using paired De Bruijn Graph
    read_pairs = k_d_mer_composition(genome, 4, 2)  # (k,d)-mer composition
    graph = de_bruijn_from_pairs(read_pairs)  # Paired de_bruijn graph
    final_genome = reconstruct_from_pairs(graph, 4, 2)
    
    
    # Compare using Hamming distance
    print("Hamming Distance between Actual Genome and De Bruijn Graph Reconstruction:")
    hamming_dist_1 = hamming_distance(genome, gen1)
    print(f"Hamming Distance: {hamming_dist_1}\n")
    
    print("Hamming Distance between Actual Genome and Paired De Bruijn Graph Reconstruction:")
    hamming_dist_2 = hamming_distance(genome, final_genome)
    print(f"Hamming Distance: {hamming_dist_2}")


    # Find the best k and d with the minimum Hamming distance
    best_k, best_d, best_score, best_d_for_k = find_best_k_d_reconstruction(genome, max_k=10, max_d=50)

    print(f"Best k = {best_k}, Best d = {best_d} with minimum Hamming distance of {best_score}")

    # Print the best d value needed for each k
    print("Best d values for each k:")
    for k, d in best_d_for_k.items():
        print(f"k = {k}: Best d = {d}")

    # **Plot the graph**
    k_values = list(best_d_for_k.keys())
    d_values = list(best_d_for_k.values())

    plt.figure(figsize=(10, 5))
    plt.plot(k_values, d_values, marker='o', linestyle='-', color='b', label="Best d for each k")

    # Labels and title
    plt.xlabel("k values")
    plt.ylabel("Best d values")
    plt.title("Best d values for each k in genome reconstruction")
    plt.xticks(k_values, rotation=45)  # Rotate x-ticks for better visibility
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend()

    # Show plot
    plt.show()
