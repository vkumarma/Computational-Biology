import random
random.seed(4269)

class Chromosome:
    def __init__(self, encoding, fitness=None):
        self.encoding = encoding  #kmer string
        self.fitness = fitness

    def __repr__(self):
        return f"Chromosome(encoding='{self.encoding}', fitness={self.fitness})"


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


def distance_between_pattern_and_string(pattern, dna):  # BA2H
    k = len(pattern)
    distance = 0
    for gene in dna:
        ham_dis = float('inf')
        for i in range(0, len(gene)-k+1):
            if ham_dis > hamming_distance(pattern, gene[i:i+k]):
                ham_dis = hamming_distance(pattern, gene[i:i+k])

        distance = distance + ham_dis

    return distance



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

def consensus_motif(motifs):
    consensus = ""
    rows = len(motifs)
    cols = len(motifs[0])

    for j in range(0, cols):
        freq = {}
        for i in range(0, rows):
            freq[motifs[i][j]] = freq.get(motifs[i][j], 0) + 1

        consensus += max(freq, key=freq.get)

    return consensus

def profile_randomly_generated_kmer(text, k, profile_matrix):  
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

    return best_motifs


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


def generate_random_kmer(k):
    bases = ['A', 'T', 'G', 'C']
    return ''.join(random.choice(bases) for _ in range(k))


def create_initial_population(population_size, k):
    population = []
    for _ in range(population_size):
        kmer = generate_random_kmer(k)
        chromosome = Chromosome(encoding=kmer)
        population.append(chromosome)
    return population


def roulette_wheel_selection(population, num_parents=2, scaling_factor=1000):
    bucket = []
    for chromosome in population:
        count = max(1, int(chromosome.fitness * scaling_factor))
        bucket.extend([chromosome] * count)

    selected_parents = [random.choice(bucket) for _ in range(num_parents)]
    return selected_parents


def crossover(parent1, parent2, k): # k point crossover
    point = random.randint(0, k)  
    child_encoding = parent1.encoding[:point] + parent2.encoding[point:]
    return Chromosome(encoding=child_encoding)


def mutate(chromosome, mutation_rate=0.05):
    
    if random.random() < mutation_rate:
        bases = ['A', 'T', 'G', 'C']
        encoding = list(chromosome.encoding)
        pos = random.randint(0, len(encoding) - 1)

        original_base = encoding[pos]
        new_base = random.choice([b for b in bases if b != original_base])

        encoding[pos] = new_base
        chromosome.encoding = ''.join(encoding)
    
    return chromosome


def genetic_algorithm(dna, k, N, pop_size, MUTATION_RATE):
    population = create_initial_population(pop_size, k)
    gen = 0

    while gen <= N:
        print(f"Generation {gen}")

        for chromosome in population:
            dist = distance_between_pattern_and_string(chromosome.encoding, dna)
            chromosome.fitness = 1 / (dist + 1e-6)

        if gen % 20 == 0:
            best_so_far = max(population, key=lambda chromo: chromo.fitness)
            print(f"Best motif at Generation {gen}: {best_so_far.encoding} (Fitness: {best_so_far.fitness})")

        parent_pool = [roulette_wheel_selection(population) for _ in range(pop_size // 2)]
        # we can directly copy some parent pool chromosomes to the next generation and create rest using crossover and mutation but in this algorithm all go through atleast crossover to get next generation

        new_population = []
        for parents in parent_pool:
            parent1, parent2 = parents
            child1 = crossover(parent1, parent2, k)
            child2 = crossover(parent2, parent1, k)
            mutate(child1, MUTATION_RATE)
            mutate(child2, MUTATION_RATE)
            new_population.extend([child1, child2])

        if len(new_population) > pop_size:
            new_population = new_population[:pop_size]

        population = new_population
        gen += 1


    # final generation results or population of chromosomes
    for chromosome in population:
        dist = distance_between_pattern_and_string(chromosome.encoding, dna)
        chromosome.fitness = 1 / (dist + 1e-6)

    final_kmers = [chromosome.encoding for chromosome in population]
    consensus = consensus_motif(final_kmers)
    consensus_fitness = 1 / (distance_between_pattern_and_string(consensus, dna) + 1e-6)

    return consensus, consensus_fitness, population



# def generate_kmer_with_limited_iupac(k):

#     allowed_chars = ['A', 'T', 'G', 'C', 'M', 'R', 'W', 'S', 'Y', 'K']
#     return ''.join(random.choice(allowed_chars) for _ in range(k))



# def weighted_selection(motifs, scores, k):  # roulette wheel selection
#     total = sum(scores)
#     if total == 0:
#         return random.choices(motifs, k=k)
    
#     probabilities = [score / total for score in scores]
#     selected = random.choices(motifs, weights=probabilities, k=k)
#     return selected



# def crossover_2(parent1, parent2):
#     k = len(parent1)
#     point = random.randint(1, k - 1)
#     child = parent1[:point] + parent2[point:]
#     return child


# def total_fitness_score(dna, d, kmer):  # fitness score for each candidate motif
#     total_score = 0
#     k = len(kmer)
    
#     for gene in dna:
#         max_match_score = 0.0  
        
#         for i in range(len(gene) - k + 1): # sequence score
#             substring = gene[i:i + k]
#             ham_dist = hamming_distance(substring, kmer)
            
#             if ham_dist == 0:
#                 match_score = 1.0  
#             elif ham_dist <= d:
#                 match_score = (k - ham_dist) / k  
#             else:
#                 match_score = 0.0
            
#             max_match_score = max(max_match_score, match_score)
#         total_score += max_match_score

#     return total_score




# def fmga_implementation(dna, k, N, pop_size):
#     population = [generate_kmer_with_limited_iupac(k) for _ in range(pop_size)]
#     index = 0
#     M = 100



if __name__ == "__main__":
    
    with open("genome.txt", "r") as file:
        genome = file.read()
    genome = "".join(genome.split())  # cleaning string

    with open("DosR.txt", "r") as file:
        data = [line.strip() for line in file]  # Removes trailing newlines
    
    
    # Parameters which can be set using config file
    KMER = 15  
    NUM_TIMES = 1000 # number of generations for ga
    INNER_LOOP, OUTER_LOOP = 2000, 50 
    MUTATION_RATE = 0.05 # for ga
    POPULATION_SIZE = 100 # for ga
    
    genetic_motifs = genetic_algorithm(data, KMER, NUM_TIMES,POPULATION_SIZE, MUTATION_RATE)
    print("Best Motif: ", genetic_motifs[0], "Fitness: ", genetic_motifs[1])

    best_motifs = genetic_motifs[2]
    for em in best_motifs:
        print("Motif: ", em.encoding, "Fitness: ", em.fitness)

    genetic_score = score([chromo.encoding for chromo in best_motifs])
    print("Genetic Algorithm Motifs Score:", genetic_score, "\n")

       
    greedy_motifs = greedy_motif_search_with_pseudocounts(data, KMER, 10)
    greedy_score = score(greedy_motifs)
    print("Greedy Search Motifs: ", greedy_motifs, "score: ", greedy_score, "\n")
    

    random_best_result = randomized_motif_search(data, KMER, 10)
    random_best_score = score(random_best_result)
    for m in range(NUM_TIMES):
        result = randomized_motif_search(data, KMER, 10)
        result_score = score(result)
        if result_score < random_best_score:
            random_best_result = result
            random_best_score = result_score

    print("Randomized Motifs: ", random_best_result, "score: ", random_best_score, "\n")

   
    gibbs_best_result = gibbs_sampler(data, KMER, 10, INNER_LOOP)
    gibbs_best_score = score(gibbs_best_result)
    for m in range(OUTER_LOOP):
        result = gibbs_sampler(data, KMER, 10, INNER_LOOP)
        result_score = score(result)
        if result_score < gibbs_best_score:
            gibbs_best_result = result
            gibbs_best_score = result_score

    print("Gibbs Sampling Motifs: ", gibbs_best_result, "score: ", gibbs_best_score, "\n")

   
   
if greedy_score < random_best_score and greedy_score < gibbs_best_score:
    final_result = greedy_motifs
    final_score = greedy_score

elif random_best_score < greedy_score and random_best_score < gibbs_best_score:
    final_result = random_best_result
    final_score = random_best_score

elif gibbs_best_score < greedy_score and gibbs_best_score < random_best_score:
    final_result = gibbs_best_result
    final_score = gibbs_best_score

else:  # if scores are same
    if greedy_score == random_best_score and greedy_score != gibbs_best_score:
        final_result = greedy_motifs + random_best_result
    elif greedy_score == gibbs_best_score and greedy_score != random_best_score:
        final_result = greedy_motifs + gibbs_best_result
        
    elif random_best_score == gibbs_best_score and greedy_score != random_best_score:
        final_result = random_best_result + gibbs_best_result
        
    else:
        final_result = greedy_motifs + random_best_result + gibbs_best_result

    consensus1 = consensus_motif(final_result)  # regular scoring function to get consensus motif
    print("Other Algorithms","Fitness: ", (1 / distance_between_pattern_and_string(consensus1, data)))
    
    freq = {}  # consensus motif minimizing distance between all collection of motifs and kmer patterns in the dna sequences
    for motif in final_result:
        freq[motif] = distance_between_pattern_and_string(motif, data)
    minimum = min(freq.values())
    for motif in freq:
        if freq[motif] == minimum:
            consensus2 = motif

    
    print("Consensus String1: ", f"'{consensus1}'","indices: ", approx_find(consensus1, genome, 3)[1])
    print("Consensus String2: ", f"'{consensus2}'","indices: ", approx_find(consensus2, genome, 3)[1])