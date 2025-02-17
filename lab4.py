def hamming_distance(p, q):  # BA1G
    mismatch = 0
    for i in range(0, len(p)):
        if p[i] != q[i]:
            mismatch += 1

    return mismatch


def profile_most_probable_kmer(text, k, profile_matrix):  # BA2C
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    k_mer_probability = {}

    for i in range(len(text) - k + 1):
        k_mer = text[i:i + k]
        p = 1

        for j in range(len(k_mer)):
            p *= profile_matrix[mapping[k_mer[j]]][j]

        k_mer_probability[k_mer] = p

    maximum = max(k_mer_probability.values())
    for key in k_mer_probability:
        if k_mer_probability[key] == maximum:
            return key


def create_profile_matrix(motifs):  # can start with initial matrix with 1.0 and normalize using (num_motifs+4) will be same as pseudocounts
    rows = 4  # base nucleotides
    cols = len(motifs[0])
    matrix = [[0.0] * cols for _ in range(rows)]
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    num_motifs = float(len(motifs))  # for normalization

    for motif in motifs:
        for i in range(len(motif)):
            matrix[mapping[motif[i]]][i] += 1.0

    for i in range(rows):
        for j in range(cols):
            matrix[i][j] = round((matrix[i][j] / num_motifs), 2)

    return matrix


def get_count_matrix(motifs):
    rows = 4  # base nucleotides
    cols = len(motifs[0])
    matrix = [[0.0] * cols for _ in range(rows)]
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    num_motifs = float(len(motifs))  # for normalization

    for motif in motifs:
        for i in range(len(motif)):
            matrix[mapping[motif[i]]][i] += 1.0

    return matrix


def create_profile_matrix_with_pseudocounts(motifs):
    counts = get_count_matrix(motifs)
    for i in range(len(counts)):
        for j in range(len(counts[0])):
            counts[i][j] += 1

    for i in range(len(counts)):
        for j in range(len(counts[0])):
            counts[i][j] = round(counts[i][j] / (len(motifs)+4), 2)

    return counts

def score(motifs):
    consensus = ""
    rows = len(motifs)
    cols = len(motifs[0])
    score = 0

    for j in range(0, cols):
        freq = {}
        for i in range(0, rows):
            freq[motifs[i][j]] = freq.get(motifs[i][j], 0) + 1

        consensus += max(freq, key=freq.get)

    for motif in motifs:
        score += hamming_distance(motif, consensus)

    return score


def greedy_motif_search(dna, k, t):  # BA2D
    best_motifs = []
    for seq in dna:
        best_motifs.append(seq[:k])

    seq_1 = dna[0]
    for i in range(len(seq_1) - k + 1):
        motifs = [seq_1[i:i + k]]
        profile = create_profile_matrix(motifs)
        for seq in dna[1:]:
            new_kmer = profile_most_probable_kmer(seq, k, profile)
            motifs.append(new_kmer)
            profile = create_profile_matrix(motifs)

        if score(motifs) < score(best_motifs):
            best_motifs = motifs

    s = ""
    for i in range(len(best_motifs)):
        s += best_motifs[i]
        if i != len(best_motifs) - 1:
            s += " "

    return s


def greedy_motif_search_with_pseudocounts(dna, k, t):  # BA2E
    best_motifs = []
    for seq in dna:
        best_motifs.append(seq[:k])

    seq_1 = dna[0]
    for i in range(len(seq_1) - k + 1):
        motifs = [seq_1[i:i + k]]
        profile = create_profile_matrix_with_pseudocounts(motifs)
        for seq in dna[1:]:
            new_kmer = profile_most_probable_kmer(seq, k, profile)
            motifs.append(new_kmer)
            profile = create_profile_matrix_with_pseudocounts(motifs)

        if score(motifs) < score(best_motifs):
            best_motifs = motifs

    s = ""
    for i in range(len(best_motifs)):
        s += best_motifs[i]
        if i != len(best_motifs) - 1:
            s += " "

    return s

with open("dna.txt", "r") as file:
    data = [line.strip() for line in file]  # Removes trailing newlines
    print(data)
    print("\n")
print(greedy_motif_search_with_pseudocounts(data, 12, 25))

