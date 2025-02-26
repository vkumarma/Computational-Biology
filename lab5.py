import random
import string


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


def get_count_matrix(motifs):
    rows = 4  # base nucleotides
    cols = len(motifs[0])
    matrix = [[0.0] * cols for _ in range(rows)]
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    for motif in motifs:
        for i in range(len(motif)):
            matrix[mapping[motif[i]]][i] += 1.0

    return matrix


def create_profile_matrix_with_pseudocounts(motifs):
    counts = get_count_matrix(motifs)
    for i in range(0, len(counts)):
        for j in range(len(counts[0])):
            counts[i][j] += 1

    for i in range(0, len(counts)):
        for j in range(len(counts[0])):
            counts[i][j] = round(counts[i][j] / (len(motifs) + 4), 2)  # for normalization

    return counts


def score(motifs):
    if not all(len(motif) == len(motifs[0]) for motif in motifs):
        print(motifs)
        raise ValueError("All motifs must have the same length!")

    consensus = ""
    rows = len(motifs)
    cols = len(motifs[0])
    score = 0

    # for i in range(len(motifs)):
    #     motifs[i] = list(motifs[i])

    for j in range(0, cols):
        freq = {}
        for i in range(0, rows):
            freq[motifs[i][j]] = freq.get(motifs[i][j], 0) + 1

        consensus += max(freq, key=freq.get)

    for motif in motifs:
        score += hamming_distance(motif, consensus)

    return score


def profile_randomly_generated_kmer(text, k, profile_matrix):  # probability of getting picket is proportional to the probability of the kmer higher probability higher likelihood of getting picked
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    k_mer_probability = {}
    big_bucket = []
    for i in range(len(text) - k + 1):
        k_mer = text[i:i + k]
        p = 1

        for j in range(len(k_mer)):
            p *= profile_matrix[mapping[k_mer[j]]][j]

        k_mer_probability[k_mer] = p

    probs_sum = sum(k_mer_probability.values())  # sum all probabilities
    for key in k_mer_probability:
        k_mer_probability[key] /= probs_sum  # normalize them so that they add upto 1

    for kmer in k_mer_probability:
        n = k_mer_probability[kmer] * 10000
        for index in range(0, int(n)):
            big_bucket.append(kmer)

    return random.choice(big_bucket)


def randomized_motif_search(dna, k, t):
    best_motifs = []
    for seq in dna:
        random_num = random.randint(0, len(seq) - k)
        best_motifs.append(seq[random_num:random_num + k])

    while True:
        new_motifs = []
        profile = create_profile_matrix_with_pseudocounts(best_motifs)
        for i in range(len(dna)):
            new_motifs.append(profile_most_probable_kmer(dna[i], k, profile))
        if score(new_motifs) < score(best_motifs):
            best_motifs = new_motifs
        else:
            return best_motifs


def gibbs_sampler(dna, k, t, N):
    best_motifs = []
    for seq in dna:
        random_num = random.randint(0, len(seq) - k)
        best_motifs.append(seq[random_num:random_num + k])

    for j in range(N):
        i = random.randint(0, t - 1)
        motifs = [best_motifs[index] for index in range(0, len(best_motifs)) if index != i]
        profile = create_profile_matrix_with_pseudocounts(motifs)
        new_kmer = profile_randomly_generated_kmer(dna[i], k, profile)
        motifs = best_motifs[:]
        motifs[i] = new_kmer
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    return best_motifs


# dna = []
# with open("dna.txt", "r") as file:
#     for line in file:
#         words = line.strip().split()  # Removes '\n' and splits by space
#         dna.extend(words)  # Adds all words to the list
#
# best_result = randomized_motif_search(dna, 15, 20)
# for m in range(1000):
#     result = randomized_motif_search(dna, 15, 20)
#     if score(result) < score(best_result):
#         best_result = result
#
# for i in best_result:
#     print(i)


dna = []
with open("dna.txt", "r") as file:
    for line in file:
        words = line.strip().split()  # Removes '\n' and splits by space
        dna.extend(words)  # Adds all words to the list

best_result = gibbs_sampler(dna, 15, 20, 2000)
for m in range(20):
    result = gibbs_sampler(dna, 15, 20, 2000)
    if score(result) < score(best_result):
        best_result = result

for i in best_result:
    print(i)
