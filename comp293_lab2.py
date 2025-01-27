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


print(len(generate_kmers(4, ['A', 'T', 'G', 'C'])))


def frequent_words_with_mismatches(text, k, d):  # BA1I
    bases = ['A', 'T', 'G', 'C']
    s = ""
    k_mers = generate_kmers(k, bases)
    k_mers_count = {}
    for kmer in k_mers:
        freq = len(approx_find(kmer, text, d)[1])
        k_mers_count[kmer] = freq

    maximum = max(k_mers_count.values())
    res = [key for key in k_mers_count if k_mers_count[key] == maximum]

    for i in range(len(res)):
        s += res[i]
        if i != len(res) - 1:
            s += " "

    return s

# print(frequent_words_with_mismatches("CTTAGATAGTCGCAAAAAAAAACTTGCCCGGCTTGCCCGGCAAAAAACACGCGGCTTTGCCCGGCTTGCCCGGCTCGCAAAAAAAAACAAAAAACTCGCAAATTGCCCGGCTTGCCCGGCCTTAGATAGTTGCCCGGCTCGCAAACTTAGATAGTTGCCCGGCACGCGGCTCTTAGATAGAAAAAACCTTAGATAGTTGCCCGGCACGCGGCTACGCGGCTAAAAAACACGCGGCTAAAAAACACGCGGCTTTGCCCGGCCTTAGATAGAAAAAACACGCGGCTCTTAGATAGACGCGGCTACGCGGCTCTTAGATAGCTTAGATAGCTTAGATAGACGCGGCTACGCGGCTCTTAGATAGTCGCAAATTGCCCGGCTCGCAAATTGCCCGGCTCGCAAATCGCAAATCGCAAAACGCGGCTTCGCAAATCGCAAATTGCCCGGCAAAAAACTCGCAAACTTAGATAGAAAAAACAAAAAACACGCGGCTTTGCCCGGCCTTAGATAGAAAAAACTTGCCCGGCTTGCCCGGCAAAAAACCTTAGATAGACGCGGCTAAAAAACACGCGGCTTCGCAAAACGCGGCTAAAAAACACGCGGCTAAAAAACTTGCCCGGCCTTAGATAGAAAAAACTCGCAAATTGCCCGGCACGCGGCTTCGCAAAAAAAAACTCGCAAAAAAAAACTTGCCCGGCAAAAAACACGCGGCTCTTAGATAGAAAAAACAAAAAACTTGCCCGGCAAAAAACCTTAGATAGTCGCAAAAAAAAACTTGCCCGGCTTGCCCGGCACGCGGCTCTTAGATAGAAAAAAC", 5, 3))

def reverse_complement(pattern):
    mapping = {'G':'C', 'C':'G', 'A':'T', 'T':'A'}
    new_string = ""
    for i in range(0, len(pattern)):
        new_string += mapping[pattern[i]]

    new_string = new_string[::-1]
    return new_string

def frequent_words_with_mismatches_and_reverse_complement(text, k, d):  # BA1J
    bases = ['A', 'T', 'G', 'C']
    s = ""
    visited = set()
    k_mers = generate_kmers(k, bases)
    k_mers_count = {}
    for kmer in k_mers:
        if kmer not in visited:
            r_c_kmer = reverse_complement(kmer)

            freq1 = len(approx_find(kmer, text, d)[1])
            freq2 = len(approx_find(r_c_kmer, text, d)[1])
            k_mers_count[kmer] = freq1 + freq2

    maximum = max(k_mers_count.values())
    res = [key for key in k_mers_count if k_mers_count[key] == maximum]

    for i in range(len(res)):
        s += res[i]
        if i != len(res) - 1:
            s += " "

    return s

print(frequent_words_with_mismatches_and_reverse_complement("ATGTTGTTGAGCTGGGATAGAGGTGTTTCTCTTTTTCTCTTGAGCTGGGATGTTGTTCGACTGCCGACTGCTTTCTCTTATGTTGTTGAGCTGGGCGACTGCATGTTGTTATAGAGGTGTTTCTCTTGAGCTGGGCGACTGCATAGAGGTGTTTCTCTTGAGCTGGGATGTTGTTGAGCTGGGATGTTGTTTTTCTCTTCGACTGCATGTTGTTTTTCTCTTGAGCTGGGGAGCTGGGCGACTGCATGTTGTTTTTCTCTTCGACTGCATAGAGGTGGAGCTGGGGAGCTGGGATGTTGTTGAGCTGGGCGACTGCTTTCTCTTATGTTGTTCGACTGCGAGCTGGGTTTCTCTTTTTCTCTTATAGAGGTGCGACTGCATAGAGGTGTTTCTCTTATAGAGGTGATAGAGGTGGAGCTGGGCGACTGCATGTTGTTATGTTGTTGAGCTGGGGAGCTGGGCGACTGCTTTCTCTTTTTCTCTTATGTTGTTCGACTGCGAGCTGGGGAGCTGGGGAGCTGGGGAGCTGGGTTTCTCTTATGTTGTTGAGCTGGGTTTCTCTTATGTTGTTGAGCTGGGGAGCTGGGATGTTGTTATAGAGGTGTTTCTCTTATAGAGGTGTTTCTCTTATAGAGGTGATAGAGGTGTTTCTCTTTTTCTCTTCGACTGCCGACTGCCGACTGCATGTTGTTCGACTGCGAGCTGGGTTTCTCTTTTTCTCTTATAGAGGTGATGTTGTTGAGCTGGGGAGCTGGGCGACTGCATGTTGTTATGTTGTTATAGAGGTGGAGCTGGGATGTTGTT", 7, 3))