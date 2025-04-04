import random
char_map = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19}

def generate_random_dna_sequence(length=20):
    return ''.join(random.choices('ACGT', k=length))

def global_alignment(v, w, score, indel_penalty):  # BA5E
    m = len(v)  # len of v, col
    n = len(w)  # len of w, row
    s = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    backtrack = [[0 for _ in range(m + 1)] for _ in range(n + 1)]

    for i in range(1, n+1):
        s[i][0] = -i * indel_penalty
    for j in range(1, m+1):
        s[0][j] = -j * indel_penalty


    for i in range(1, n+1): # 1 -> n
        for j in range(1, m+1):  # 1 -> m
            match_score = score[char_map[w[i-1]]][char_map[v[j-1]]]

            diag = s[i-1][j-1] + match_score
            down = s[i-1][j] - indel_penalty
            right = s[i][j-1] - indel_penalty

            s[i][j] = max(diag, down, right)

            if s[i][j] == diag:
                backtrack[i][j] = "diag"
            if s[i][j] == down:
                backtrack[i][j] = "down"
            if s[i][j] == right:
                backtrack[i][j] = "right"


    # Backtracking
    top = ""
    bottom = ""
    i = len(w)
    j = len(v)
    while i > 0 and j > 0:
        if backtrack[i][j] == "down":
            top = "-" + top
            bottom = w[i-1] + bottom
            i -= 1

        if backtrack[i][j] == "right":
            top = v[j-1] + top
            bottom = "-" + bottom
            j -= 1

        if backtrack[i][j] == "diag":
            top = v[j-1] + top
            bottom = w[i-1] + bottom
            i -= 1
            j -= 1

    while i > 0:
        top = "-" + top
        bottom = w[i-1] + bottom
        i -= 1

    while j > 0:
        top = v[j-1] + top
        bottom = "-" + bottom
        j -= 1

    return s[n][m], top, bottom
    

    


def local_alignment(v, w, score, indel_penalty):  # BA5F
    m = len(v)  # len of v, col
    n = len(w)  # len of w, row
    s = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    backtrack = [[None for _ in range(m + 1)] for _ in range(n + 1)]

    for i in range(1, n + 1):
        for j in range(1, m + 1):

            diag = s[i - 1][j - 1] + score[char_map[w[i - 1]]][char_map[v[j - 1]]]
            down = s[i - 1][j] - indel_penalty
            right = s[i][j - 1] - indel_penalty

            s[i][j] = max(diag, down, right, 0)

            if s[i][j] == 0:
                backtrack[i][j] = "end"
            elif s[i][j] == diag:
                backtrack[i][j] = "diag"
            elif s[i][j] == down:
                backtrack[i][j] = "down"
            elif s[i][j] == right:
                backtrack[i][j] = "right"

    max_val = 0
    max_i = 0
    max_j = 0
    for i in range(0, n + 1):
        for j in range(0, m + 1):
            if s[i][j] > max_val:
                max_val = s[i][j]
                max_i = i
                max_j = j


    #Backtracking
    aligned_v = ""
    aligned_w = ""
    i, j = max_i, max_j
    while i > 0 and j > 0:
        if backtrack[i][j] == "diag":
            aligned_v = v[j - 1] + aligned_v
            aligned_w = w[i - 1] + aligned_w
            i -= 1
            j -= 1
        if backtrack[i][j] == "down":
            aligned_v = "-" + aligned_v
            aligned_w = w[i - 1] + aligned_w
            i -= 1
        if backtrack[i][j] == "right":
            aligned_v = v[j - 1] + aligned_v
            aligned_w = "-" + aligned_w
            j -= 1
        if backtrack[i][j] == "end":
            break

    return max_val, aligned_v, aligned_w



def fitting_alignment(v, w, score, indel_penalty):  # BA5H
    m = len(v)  # len of v, col
    n = len(w)  # len of w, row
    s = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    backtrack = [[None for _ in range(m + 1)] for _ in range(n + 1)]

    for i in range(0, n + 1):  # first column has gap penalty
        s[i][0] = -i * indel_penalty

    for i in range(1, n + 1):  # iterating matrix
        for j in range(1, m + 1):

            # Get the score for the match/mismatch
            match_score = score[char_map[w[i - 1]]][char_map[v[j - 1]]]

            diag = s[i - 1][j - 1] + match_score  
            down = s[i - 1][j] - indel_penalty 
            right = s[i][j - 1] - indel_penalty  

            s[i][j] = max(diag, down, right) 

            # Keep track of the direction from which the score came
            if s[i][j] == diag:
                backtrack[i][j] = "diag"
            elif s[i][j] == down:
                backtrack[i][j] = "down"
            elif s[i][j] == right:
                backtrack[i][j] = "right"

    # Finding the maximum score in the last row, as this is fitting alignment
    max_val = s[n][0]
    max_j = len(v)

    for j in range(0, m + 1):  # iterating last row columns
        if s[-1][j] > max_val:
            max_val = s[-1][j]
            max_j = j

    # Backtracking to retrieve the optimal alignment
    aligned_v = ""
    aligned_w = ""
    i, j = n, max_j
    while i > 0:
        if backtrack[i][j] == "diag":
            aligned_v = v[j - 1] + aligned_v
            aligned_w = w[i - 1] + aligned_w
            i -= 1
            j -= 1
        elif backtrack[i][j] == "down":
            aligned_v = "-" + aligned_v
            aligned_w = w[i - 1] + aligned_w
            i -= 1
        elif backtrack[i][j] == "right":
            aligned_v = v[j - 1] + aligned_v
            aligned_w = "-" + aligned_w
            j -= 1

    return max_val, aligned_v, aligned_w




def overlap_alignment(v, w, indel, mismatch, match):
    n, m = len(v), len(w)

    # Initialize score and backtrack matrices with all zeros
    score_mat = [[0] * (m + 1) for _ in range(n + 1)]
    backtrack_mat = [[None] * (m + 1) for _ in range(n + 1)]

    # Fill the matrices
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_score = match if v[i - 1] == w[j - 1] else mismatch
            scores = [
                score_mat[i - 1][j - 1] + match_score,  
                score_mat[i - 1][j] + indel,  
                score_mat[i][j - 1] + indel  
            ]
            score_mat[i][j] = max(scores)
            backtrack_mat[i][j] = ['diag', 'down', 'right'][scores.index(score_mat[i][j])]

    
    max_score, max_j = max((score_mat[n][j], j) for j in range(m + 1))
    i, j = n, max_j

    
    aligned_v, aligned_w = "", ""
    while i > 0 and j > 0:
        if backtrack_mat[i][j] == 'diag':
            aligned_v = v[i - 1] + aligned_v
            aligned_w = w[j - 1] + aligned_w
            i, j = i - 1, j - 1
        elif backtrack_mat[i][j] == 'down':
            aligned_v = v[i - 1] + aligned_v
            aligned_w = "-" + aligned_w
            i -= 1
        elif backtrack_mat[i][j] == 'right':
            aligned_v = "-" + aligned_v
            aligned_w = w[j - 1] + aligned_w
            j -= 1

    return max_score, aligned_v, aligned_w



def alignment_information(seq1, seq2):  # gives sequences alignment information (matches, indels, mismatches)
    matches = 0
    indels = 0
    mismatches = 0
    
    for i in range(len(seq1)):
        s1 = seq1[i]
        s2 = seq2[i]
        
        if s1 == s2:  
            matches += 1
        elif '-' in [s1, s2]:  
            indels += 1
        else:  
            mismatches += 1
    
    return matches, indels, mismatches


def find_best_alignment_for_each_sequence(sequences, alignment_func_type, score_matrix, indel_penalty, human_sequence=None):
    # Map alignment function names to their respective functions
    alignment_functions = {
        'global': global_alignment,
        'local': local_alignment,
        'fitting': fitting_alignment
    }

    if alignment_func_type not in alignment_functions:
        raise ValueError(f"Invalid alignment function type: {alignment_func_type}. Expected 'global', 'local', or 'fitting'.")

    
    alignment_func = alignment_functions[alignment_func_type]
    
    best_alignments = {}

    
    if alignment_func_type == 'fitting' and human_sequence is None:
        raise ValueError("For fitting alignment, the human sequence must be provided.")

    # Compare every sequence with every other sequence (including itself)
    if alignment_func_type == 'fitting':
        for org, seq in sequences.items():
            score, alignment_v, alignment_w = alignment_func(human_sequence, seq, score_matrix, indel_penalty)
            matches, indels, mismatches = alignment_information(alignment_v, alignment_w)
            
            best_alignments[(human_sequence, org)] = {
                "alignment": (alignment_v, alignment_w),
                "score": score,
                "matches": matches,
                "indels": indels,
                "mismatches": mismatches
            }
    
    else:
       
        for org1, seq1 in sequences.items():
            max_score = float('-inf')
            best_alignment = None
            best_pair = None

            for org2, seq2 in sequences.items():
                if org1 == org2:
                    continue  

                score, alignment_v, alignment_w = alignment_func(seq1, seq2, score_matrix, indel_penalty)
                
                if score > max_score:
                    max_score = score
                    best_alignment = (alignment_v, alignment_w)
                    best_pair = (org1, org2)
                    
            if best_pair:
                matches, indels, mismatches = alignment_information(*best_alignment)

                best_alignments[best_pair] = {
                    "alignment": best_alignment,
                    "score": max_score,
                    "matches": matches,
                    "indels": indels,
                    "mismatches": mismatches
                }

    return best_alignments

 
def apply_overlap_alignment(seq1, seq2, indel_penalty, mismatch_penalty, match_score):
    # Alignment of suffix of seq1 to prefix of seq2
    suffix_seq1_prefix_seq2 = overlap_alignment(seq1[-10:], seq2[:10], indel_penalty, mismatch_penalty, match_score)

    # Alignment of suffix of seq2 to prefix of seq1
    suffix_seq2_prefix_seq1 = overlap_alignment(seq2[-10:], seq1[:10], indel_penalty, mismatch_penalty, match_score)

    return suffix_seq1_prefix_seq2, suffix_seq2_prefix_seq1


if __name__ == "__main__":
    score = []
    with open('matrix.txt', 'r') as file:
        for line in file:
            values = line.split()
            row = list(map(int, values[1:]))
            score.append(row)

    sequences = {}
    with open("sequences.txt", "r") as file:
        current_key = None
        current_sequence = ""

        for line in file:
            line = line.strip()
            if line.startswith(">"):  
                if current_key: 
                    sequences[current_key] = current_sequence
                current_key = line[1:] 
                current_sequence = ""  
            else:
                current_sequence += line  

        if current_key: 
            sequences[current_key] = current_sequence

   
    human_sequence = "MGDVEKGKKIFIMKCSQCHTVEKGGKHKTGPNLHGLFGRKTGQAPGYSYTAANKNKGIIWGEDTLMEYLENPKKYIPGTKMIFVGIKKKEERADLIAYLKKATNE"  # conserved

    indel_penalty = 5
    alignment_type = "global" # "global",  "local"
    best_alignments = find_best_alignment_for_each_sequence(sequences, alignment_type, score, indel_penalty, human_sequence)

    
    for pair, result in best_alignments.items():
        org1, org2 = pair
        alignment_v, alignment_w = result['alignment']
        print(f"Best alignment for pair: {org1} vs {org2}")
        print(f"Alignment Type: {alignment_type.capitalize()}")
        print(f"Score: {result['score']}")
        print(f"Matches: {result['matches']}")
        print(f"Mismatches: {result['mismatches']}")
        print(f"Indels: {result['indels']}")
        print(f"Length of alignment: {len(alignment_v)}")
        print("Alignment:")
        print(alignment_v)
        print(alignment_w)
        print("-" * 120)



        # Generate random sequences
    seq1 = generate_random_dna_sequence(20)
    print(seq1)
    seq2 = generate_random_dna_sequence(20)
    print(seq2)

    # Test with different parameter sets
    params = [
        (-5, -1, 1),   # Indel = -5, Mismatch = -1, Match = 1
        (-1, -5, 1),   # Indel = -1, Mismatch = -5, Match = 1
        (-1, -1, 5)    # Indel = -1, Mismatch = -1, Match = 5
    ]

    # Apply the function for each parameter set
    for indel, mismatch, match in params:
        print(f"Parameters - Indel: {indel}, Mismatch: {mismatch}, Match: {match}")
        align1, align2 = apply_overlap_alignment(seq1, seq2, indel, mismatch, match)
        print("Alignment 1 (Suffix seq1 to Prefix seq2):")
        print(f"Score: {align1[0]}, Alignment: {align1[1]} | {align1[2]}")
        print("Alignment 2 (Suffix seq2 to Prefix seq1):")
        print(f"Score: {align2[0]}, Alignment: {align2[1]} | {align2[2]}")
        print('-' * 50)


    



