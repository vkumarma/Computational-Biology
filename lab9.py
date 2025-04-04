char_map = {'A': 0, 'C': 1, 'D': 2, 'E': 3, 'F': 4, 'G': 5, 'H': 6, 'I': 7, 'K': 8, 'L': 9, 'M': 10, 'N': 11, 'P': 12, 'Q': 13, 'R': 14, 'S': 15, 'T': 16, 'V': 17, 'W': 18, 'Y': 19}



file_path = 'matrix.txt'
score = []
with open(file_path, 'r') as file:
    for line in file:
        values = line.split()
        row = list(map(int, values[1:]))
        score.append(row)

# for row in score:
#     print(row)
#
# print(score[char_map['V']][char_map['W']])




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

    return s[n][m], backtrack

def global_output_alignment(backtrack, v, w):
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

    return top, bottom


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



def fitting_alignment(v, w, indel_penalty):  # BA5H
    m = len(v)  # len of v, col
    n = len(w)  # len of w, row
    s = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    backtrack = [[None for _ in range(m + 1)] for _ in range(n + 1)]

    for i in range(0, n + 1):  # first column has gap penalty
        s[i][0] = -i * indel_penalty

    for i in range(1, n + 1):  # iterating matrix
        for j in range(1, m + 1):

            diag = s[i - 1][j - 1]
            down = s[i - 1][j]
            right = s[i][j - 1]

            if  w[i - 1] == v[j - 1]: # if match
                diag = s[i - 1][j - 1] + indel_penalty
                down = s[i - 1][j] - indel_penalty
                right = s[i][j - 1] - indel_penalty

            elif w[i - 1] != v[j - 1]: # if mis-match
                diag = s[i - 1][j - 1] - indel_penalty
                down = s[i - 1][j] - indel_penalty
                right = s[i][j - 1] - indel_penalty

            s[i][j] = max(diag, down, right)
            if s[i][j] == diag:
                backtrack[i][j] = "diag"
            elif s[i][j] == down:
                backtrack[i][j] = "down"
            elif s[i][j] == right:
                backtrack[i][j] = "right"

    max_val = s[n][0]
    max_j = len(v)

    for j in range(0, m+1):  # iterating last row columns
        if s[-1][j] > max_val:
            max_val = s[-1][j]
            max_j = j

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
                score_mat[i - 1][j - 1] + match_score,  # Diagonal
                score_mat[i - 1][j] + indel,  # Up (insert gap in w)
                score_mat[i][j - 1] + indel  # Left (insert gap in v)
            ]
            score_mat[i][j] = max(scores)
            backtrack_mat[i][j] = ['d', 'u', 'l'][scores.index(score_mat[i][j])]

    # Find the position with the maximum score in the last row or column
    max_score, max_j = max((score_mat[n][j], j) for j in range(m + 1))
    i, j = n, max_j

    # Traceback to construct the alignment
    aligned_v, aligned_w = "", ""
    while i > 0 and j > 0:
        if backtrack_mat[i][j] == 'd':
            aligned_v = v[i - 1] + aligned_v
            aligned_w = w[j - 1] + aligned_w
            i, j = i - 1, j - 1
        elif backtrack_mat[i][j] == 'u':
            aligned_v = v[i - 1] + aligned_v
            aligned_w = "-" + aligned_w
            i -= 1
        elif backtrack_mat[i][j] == 'l':
            aligned_v = "-" + aligned_v
            aligned_w = w[j - 1] + aligned_w
            j -= 1

    return max_score, aligned_v, aligned_w


v = "GACTGGTTTAAACCAAAGATACGATGGTTCTCCAAGCCATCCGGACTATAAAGGTGCCGGTGGAACTATCGGCTTGGTTTTCGGAGGTGTGACGAGCAAGGTCATCGTTCGGGGGCGCGCGGCTGCGCAGATGCCAACGGTATAGTACACATTATCATGTCTCTTCGCATCTTATTCACCACGGTCTCAGTAGAACTTCTAGGGATATACAAATCACCGTTTCGGTGTGTTCTTCCACAGGGGTTATTCACGTTGGTACATGCCCGGTTGTTCCACTGAATGGCGAGGAGTTGAAGTATACTCATATAGATAGGGAGAGGACCGCCTATGGACCTGTACCCGACCTGGCGTTGTGCCAGACTTATCTCTAACAAATGTGAAAGCAAGCACCCACAGCCATACGGACAACAGGCGTGAATGCCAACCAATCCCTCTCCCGCGTGCTACAGGCAGAAAAGTCATGAAATCGACTACCACCACCGAGACGAACTAGCCCGAAGAAATTGCGCGCCGTTGGAGGGTGTCTAACGATGACGTAGAGCAAATGCTAGCGATACGCGGAGTCAGATCGACCCGGGCGTCAAAGGCCAATCGCAAGTTTTCTCATAGAATTCCATCCCGGTTAGTCATACGGCGCGTAGAGCTCGCACACGTATAGCTGCAGTCCGTATGGCCCAGCTTTGTTCCCCATGATGGGTTTCCCGTGTATTCGGCAAGAGTGATTGTACTGCTGATGGTTGTATGGGTTGGTCGCCTCTGAGTGATAGCATTAGTTTCGTTTGTACGCGTAACTAGGACTGGCACGTGCGACAGCCATGGGTCACCCGCTGGTACCTTCAATCGTACGCGTGTCTTCG"
w = "ATCGCAGCCTTCTCATGGGTCACAGCCCCGGTTAGTCGTCGCGCCCGCAGTGCTCGCACACCGATTAACCTGCAGCTCGTGGCCCGACTTTGTCCTCAGGGTGGGGTGTTCCTGGTATTTAGGCATAAGTGGTTACTCGCTGATGTTGTAAGGGTGCCGTGGCTCTGACCCGCTACATTAACATAGCTTGGACGGTAACTAAGGCAGGCCCTTCAGCAGCCACGGGTAAATCCGGTGAGTCCCTTTAAGTATAGCTTTATTGCCTACAGCGAACTGAAGGACCCTCCTGACCGAAGGTATATATGAAAACCGTACTGTGCACCAATTATACACTGCATCGGTTCTACCGTTTTTTGACTGACAGGGGCATTGGGTCTGCTATACGCGAGGAGCGAGAGCGCAGGGCCAATCCCTAGCTATTTGGACGAAGCCACACCCTACTGCAGCAAGGCGTAGGGCCCCCGACCAAACCACGTTTAATATTACACTTTTAAGACGAAGGGCAAATATTTCTCGGTTAAGACCACCTGACGTTAATCTCTATGCGCCTGAACCATCCACTAACCCACGCACTTGTTAAGGTTCTCTATGTCCCAACCGCGGTAACTCCTCTCGCAAACACTCGGTCCTGCGTGTGTATCGACCACTACTAGAGCAGACCGTTAACGCGTCCTCCACTATATCAGTGCTCGCCCCGATTAGGAGTCCGTACGAACTCCCAAGTTAGCTGTAACCTCGTACTGAGGTAAATCCTGTGAAATTTCGCTGCACAAATTGCGGACGCCTAGTCCGGATTGCGCTCCCGACTTATTAATGACGTGGTAATATCTTACTATGTTGATTTTGCTCGTACACGT"

s, str1, str2 = overlap_alignment(v, w, -5,-1, 1)
print(s)
print(str1)
print(str2)

