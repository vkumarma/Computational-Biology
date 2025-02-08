import os
from collections import defaultdict

def hamming_distance(p, q):  # BA1G
    mismatch = 0
    for i in range(0, len(p)):
        if p[i] != q[i]:
            mismatch += 1

    return mismatch


def approx_find(pattern, text, d):  # BA1H
    indices = []
    s = ""
    for i in range(0, len(text) - len(pattern) + 1):
        if hamming_distance(text[i:i + len(pattern)], pattern) <= d:
            indices.append(i)

    for i in range(len(indices)):
        s += str(indices[i])
        if i != len(indices) - 1:
            s += " "
    return s, indices


def generate_kmers(k, bases):
    if k == 1:
        return bases

    small_kmers = generate_kmers(k - 1, bases)
    k_mers = []
    for kmer in small_kmers:
        for b in bases:
            k_mers.append(kmer + b)
    
    
    return k_mers


def reverse_complement(pattern):
    mapping = {'G':'C', 'C':'G', 'A':'T', 'T':'A'}
    new_string = ""
    for i in range(0, len(pattern)):
        new_string += mapping[pattern[i]]

    new_string = new_string[::-1]
    return new_string


def frequent_words_with_mismatches_and_reverse_complement(text, k, d, kmers):  # BA1J
    visited = set()
    k_mers = kmers
    k_mers_count = {}
    res = []
    for kmer in k_mers:
        if kmer not in visited:
            r_c_kmer = reverse_complement(kmer)
            visited.add(r_c_kmer)
            freq1 = len(approx_find(kmer, text, d)[1])
            freq2 = len(approx_find(r_c_kmer, text, d)[1])
            k_mers_count[kmer] = freq1 + freq2


    maximum = max(k_mers_count.values())
    for key in k_mers_count:
        if k_mers_count[key] == maximum:
            res.append((key, k_mers_count[key]))
            res.append((reverse_complement(key), k_mers_count[key]))

   
    return res


def skew(genome):
    arr = [0] * (len(genome) + 1)
    g_s, c_s = 0, 0
    for i in range(0, len(genome)):
        if genome[i] == 'A' or genome[i] == 'T':
            arr[i + 1] = arr[i]
        else:
            if genome[i] == 'G':
                g_s += 1

            elif genome[i] == 'C':
                c_s += 1

            arr[i + 1] = g_s - c_s

    minimum = min(arr)
    indices = [i for i in range(0, len(arr)) if arr[i] == minimum]
    return indices


def get_batches(sequence, index, batch_size=500):
    left_batch = sequence[max(0, index - batch_size):index]
    right_batch = sequence[index + 1:index + 1 + batch_size]
    
    return left_batch, right_batch

def create_batches_for_indexes(sequence, indices):
    return [batch for index in indices for batch in get_batches(sequence, index)]


def potential_DnaA_boxes(sequences, k_mers):
    dnaA_boxes = set()
    for i in range(len(sequences)):
        patterns = frequent_words_with_mismatches_and_reverse_complement(sequences[i], 9, 2, k_mers)
        for pattern in patterns:
            dnaA_boxes.add(pattern)
        
        print("Potential DnaA Boxes: ", dnaA_boxes)

    return dnaA_boxes


if __name__ == "__main__":
    with open("genome.txt", "r") as file:
        genome = file.read()

    genome = "".join(genome.split())  # cleaning string
    indices = skew(genome)
    potential_oris = create_batches_for_indexes(genome,indices)
    generated_kmers = generate_kmers(9,['A', 'T', 'G', 'C'])  # all 9-mers 
    print("Final Potential DnaA Boxes: ", potential_DnaA_boxes(potential_oris, generated_kmers))