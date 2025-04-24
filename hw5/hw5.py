import time
import random
from collections import defaultdict

class TrieNode:
    def __init__(self):
        self.children = {}
        self.is_word = False
        self.edge = [0, 0]

class Trie:
    def __init__(self):
        self.root = TrieNode()

    def add(self, patterns):
        num_nodes = 0
        for pattern in patterns:
            prev = None
            curr_node = self.root
            for char in pattern:
                if char not in curr_node.children:
                    num_nodes += 1
                    curr_node.children[char] = TrieNode()
                    temp = curr_node.children[char]
                    if prev is None:
                        temp.edge[1] = num_nodes
                    elif prev is not None:
                        temp.edge[0] = curr_node.edge[1]
                        temp.edge[1] = num_nodes

                prev = curr_node
                curr_node = curr_node.children[char]
            curr_node.is_word = True

def dfs(node, file_path='output.txt'):  # BA9A had to change because recursive dfs has a recursive depth limit.
    with open(file_path, 'w') as f:
        stack = [node]

        while stack:
            node = stack.pop()
            for key, child in node.children.items():
                print(f"{child.edge[0]}->{child.edge[1]}:{key}")
                f.write(f"{child.edge[0]}->{child.edge[1]}:{key}\n")
                stack.append(child)


def generate_kmers(k, bases):
    if k == 1:
        return bases

    small_kmers = generate_kmers(k - 1, bases)
    k_mers = []
    for kmer in small_kmers:
        for b in bases:
            k_mers.append(kmer + b)

    return k_mers


def generate_random_dna(length):
    return ''.join(random.choices('ATCG', k=length))
        
def suffix_array(text, file_path='suffix_array.txt'): # BA9G
    # seqs = dict((i, text[i:]) for i in range(len(text)))
    # return sorted(seqs.keys(), key=lambda x: seqs[x])

    suffixes = {i: text[i:] for i in range(len(text))}
    sorted_indices = sorted(suffixes.keys(), key=lambda x: suffixes[x])

    with open(file_path, 'w') as f:
        for index in sorted_indices:
            f.write(f"{index}\t{suffixes[index]}\n")

    print(f"Suffix array written to {file_path}")
    return sorted_indices


def evaluate_compression(original_seq, rle_seq):
    original_size = len(original_seq)
    compressed_size = len(rle_seq)
    compression_ratio = compressed_size / original_size

    print(f"Original size: {original_size} characters")
    print(f"Compressed size: {compressed_size} characters")
    print(f"Compression ratio: {compression_ratio:.2f}")

    if compressed_size < original_size:
        print("✅ RLE compression was effective.")
    else:
        print("❌ RLE compression was not effective.")


def create_bwt(text, output_file="bwt_output.txt"):
    cyclic_permutations = []
    for i in range(0, len(text)):
        cyclic_permutations.append(text[i:] + text[:i])
    cyclic_permutations = sorted(cyclic_permutations)
    bwt = ''.join(st[-1] for st in cyclic_permutations)

    with open(output_file, "w") as file:
        file.write(bwt)
    return bwt


def reconstruct_text(bwt_string):
    last_col = list(bwt_string)
    first_col = sorted(last_col)

    text = "$"
    idx = 0
    next_symbol = last_col[idx]
    while next_symbol != "$":
        text = next_symbol + text
        symbol_count = 0
        for sym in last_col[:idx+1]:
            if sym == next_symbol:
                symbol_count += 1

        count = 0
        for i in range(1, len(first_col)):
            if first_col[i] == next_symbol:
                count += 1
                if count == symbol_count:
                    idx = i
                    break

        next_symbol = last_col[idx]

    return text


def preprocess_bwt(bwt):
    first_col = sorted(bwt)
    first_occurrence = {}
    for i, char in enumerate(first_col):
        if char not in first_occurrence:
            first_occurrence[char] = i

    count_matrix = []
    tally = defaultdict(int)
    for char in bwt:
        tally[char] += 1
        count_matrix.append(tally.copy())

    def annotate(column):
        counts = {}
        result = []
        for char in column:
            counts[char] = counts.get(char, 0) + 1
            result.append((char, counts[char]))
        return result

    last_annotated = annotate(bwt)
    first_annotated = annotate(first_col)

    last_to_first = {}
    used = defaultdict(list)
    for idx, item in enumerate(first_annotated):
        used[item].append(idx)
    for idx, item in enumerate(last_annotated):
        last_to_first[idx] = used[item].pop(0)

    return first_occurrence, count_matrix, last_to_first


def count_symbol(count_matrix, symbol, row):
    if row < 0:
        return 0
    return count_matrix[row].get(symbol, 0)


def bw_matching(bwt, pattern, first_occurrence, count_matrix, last_to_first):
    top = 0
    bottom = len(bwt) - 1
    while top <= bottom:
        if pattern:
            symbol = pattern[-1]
            pattern = pattern[:-1]
            top_count = count_symbol(count_matrix, symbol, top - 1)
            bottom_count = count_symbol(count_matrix, symbol, bottom)
            if bottom_count - top_count > 0:
                top = first_occurrence[symbol] + top_count
                bottom = first_occurrence[symbol] + bottom_count - 1
            else:
                return 0
        else:
            return bottom - top + 1
    return 0


def match_all(bwt_string, patterns):
    first_occurrence, count_matrix, last_to_first = preprocess_bwt(bwt_string)
    results = {}
    for pat in patterns:
        count = bw_matching(bwt_string, pat, first_occurrence, count_matrix, last_to_first)
        results[pat] = count
    return results


def run_length_encoding(sequence):
    compressed = ""
    i = 0

    while i < len(sequence):
        count = 1
        j = i + 1
        while j < len(sequence) and sequence[i] == sequence[j]:
            count += 1
            j += 1

        compressed += str(count)
        compressed += sequence[i]
        i = j  

    return compressed



if __name__ == "__main__":
    sequence = ""

    with open('sequence.txt', 'r') as file:
        for line in file:
            sequence += line.strip()


    suffixes = []
    for i in range(len(sequence)):
        suffix = sequence[i:]
        suffixes.append(suffix)


    start_time1 = time.time()
    my_trie = Trie()
    my_trie.add(suffixes)
    end_time1 = time.time()
    elapsed_time = end_time1 - start_time1
    print(f"Suffix trie created in {elapsed_time:.6f} seconds") 
    # dfs(my_trie.root) Takes too long


    start_time2 = time.time()
    suf_array = suffix_array(sequence)
    end_time2 = time.time()
    elapsed_time = end_time2 - start_time2
    print(f"Suffix Array created in {elapsed_time:.6f} seconds")


    start_time3 = time.time()
    bwt_string = create_bwt(sequence)
    end_time3 = time.time()
    elapsed_time = end_time3 - start_time3
    print(f"BWT Array created in {elapsed_time:.6f} seconds")


    rle = run_length_encoding(sequence)
    print("rle: ",len(rle), " ", "sequence: ",len(sequence))

    random_dna = generate_random_dna(len(sequence))
    random_dna_rle = run_length_encoding(random_dna)
    print(evaluate_compression(random_dna, random_dna_rle))


    kmers_occurences = {}
    kmers = generate_kmers(5, ['A', 'T', 'G', 'C'])
    start_time4 = time.time()
    for kmer in kmers:
        count = 0
        
        for i in range(0, len(sequence) - len(kmer) + 1):
            if sequence[i:i + len(kmer)] == kmer:
                count += 1
        kmers_occurences[kmer] = count
    end_time4 = time.time()
    elapsed_time = end_time4 - start_time4
    print(f"Linear Search in {elapsed_time:.6f} seconds")
    print(kmers_occurences, "\n")
    

    start_time5 = time.time()
    occurences = match_all(bwt_string, kmers)
    end_time5 = time.time()
    elapsed_time = end_time5 - start_time5
    print(f"BWT Search in {elapsed_time:.6f} seconds")
    print(print(occurences), "\n")