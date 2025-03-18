def manhattan_tourist_problem(n, m, down, right):
    s = [[0] * (m + 1) for _ in range(n + 1)]
    s[0][0] = 0

    for i in range(1, n+1):
        s[i][0] = down[i-1][0] + s[i-1][0]

    for j in range(1, m+1):
        s[0][j] = right[0][j-1] + s[0][j-1]

    for i in range(1, n+1):
        for j in range(1, m+1):
            s[i][j] = min(s[i-1][j] + down[i-1][j], s[i][j-1] + right[i][j-1])

    return s[n][m]


def format_input(data):
    data_parts = data.split("-")

    first_part = data_parts[0].strip().splitlines()
    second_part = data_parts[1].strip().splitlines()

    n, m = [int(x) for x in first_part[0].split()]
    down_matrix = []
    for i in range(1, n + 1):
        down_matrix.append([int(x) for x in first_part[i].split()])

    right_matrix = []
    for line in second_part:
        right_matrix.append([int(x) for x in line.split()])

    return n, m, down_matrix, right_matrix


def longestCommonSubsequence(text1, text2):
    string_a = text1
    string_b = text2
    matrix = [[None for _ in range(len(string_a) + 1)] for _ in range(len(string_b) + 1)]


    for i in range(len(string_a) + 1):
        matrix[0][i] = 0
    for j in range(len(string_b) + 1):
        matrix[j][0] = 0


    for i in range(1, len(string_b) + 1):
        for j in range(1, len(string_a) + 1):
            if string_b[i - 1] == string_a[j - 1]:
                matrix[i][j] = matrix[i - 1][j - 1] + 1
            else:
                matrix[i][j] = max(matrix[i][j - 1], matrix[i - 1][j])


    # backtracking
    i, j = len(string_b), len(string_a)
    lcs = []
    while i > 0 and j > 0:
        if string_b[i - 1] == string_a[j - 1]:
            lcs.append(string_b[i - 1])
            i -= 1
            j -= 1
        elif matrix[i - 1][j] > matrix[i][j - 1]:
            i -= 1
        else:
            j -= 1


    return ''.join(reversed(lcs))


input_data = """11 4
4 2 3 2 4
0 3 1 2 3
4 0 3 1 0
4 4 0 1 3
1 4 0 1 3
1 0 3 0 0
4 1 3 0 0
4 4 4 4 1
0 3 2 4 4
2 0 2 3 4
0 2 4 1 3
-
4 0 0 0
0 3 4 1
3 2 0 3
1 3 2 0
0 2 2 0
4 3 0 2
3 3 0 4
1 1 0 1
0 4 2 1
2 0 3 1
2 4 2 2
0 4 4 2
"""


n, m, down_matrix, right_matrix = format_input(input_data)
print(manhattan_tourist_problem(n, m, down_matrix, right_matrix))


text1 = "GTGGCTGGCACATCGTAGTACCACAACGTTGACTGAGTTATGTGCGGAGTGCCAGTTTGAGGGCAGTGGCCGACTGGAAAGTATAGCTTGCGATGCAGTGGTGTGCACGGAAGGCTTTTCCGGCGCCCAAGTTTCGTTTTGCAACAAAGGCTCAACCCGAAAATCCACTGCCATTGCTAGCTAGCTACCAAGGCAGGTAACCAGTGTGTTTGGTGATCTGCGGCGCGTGGGTAGGCCTTGGAGTGGCAATTTAGCGATCTGCCCTCTGTGGATTCCGAGCGCTGCCCTGAATGTCCACCAAAACGCCGTTACAACTGAGAGGTTTTAGCCTACATATGAACTAGCCAATGGGGGCGTTGATATGTACGTATCTGCGAATACGCTTGGTGTCTCCTTCATTTCCGATTTTTTTAAAGTGACCATTCTCCCATGAGTAGCGTAGCGCTATTGGACTAAGGACATAAACAGTGCGTTTTCCGTGGTAATTCGGCAACTTACGTTGTTCTTTTATTATGGATGAGAAGTAGGCCCCTTCCGCTAGATTGGGCAACGAGAGTCTAATCATATTCATAATGGAGCACTTCACTGCAAACCTAGAAGGTGTGAATTTGACGAGAACGTTTGTAGGGTTAGTGGATTAAACCCATCCTATACCATTAGAATGTGTTGGGCCACAGCAGTGTGATGTCCGTGGTCGGGAGCAAAATTGTTTAATAGCCCCTTCGAATTGACTCAAATTTTTGCGTTTCTAATTCCCTAAAAAATGCTACGAGACTGCGTACATGTGAAATGCCCGGGCGTTGTCAC"
text2 = "GTGGGCAGACTCTTCCGCCCTTCCTAGTTGCGAACGGCTTAATTGACACCAAGAATTGAAGTCTCGACTTCACCCAAACTCGACTAGTTGGTTTTCCACGTCGGGCAGATGCGTGACTGCAGTTACCACAACAAAGCGCCGGCTGCCTAAGTCGGCTCACAAGCATCAGCTTATCCGACCCTCGGTAGACAAGATTCAGGTAAGGCGACTCGACTACACTGGAGAACTAAGAATGAGGGAGTATTTATGGCACGAGACACTCTCTCTCAATTCAATATCTTGCCGGAGGGCTGTTAGGTACAGATCCTGGTCAGGTTACACCAGTATCGCGTCGGAACCGTGTATTAGTTTAAAACCTCGCCACAGTGAATTCAAATGTGTTTAAAGATGCGTAGGTACTTCAGCCGGAATATTGTTCCAGTCTCACCGAGGTCCCGCGTGAGTTGGTTATCACGCTGGCCAAGCAGGCGACAGCGTCTCCGCATAGAGAGGACAACATCCCTGATTGCCATAACGTCGAGGTCCCTGTGGTAACCAGACATACTGTAACTGGTACGGCTGGATCTGAAGCTCCCGCCAACAATACAAGACTTCTGGCGTGGGTTGAGCAGTGTCTAGAACACTTTCGACCTGAGTTATATTCGCAGGAATTATCCGGTTATTGAGTAAAACAGATACTGAGGATTGGGTTTTGGGGTTATGAAGTGTCGGTAATACTAGCAGTGGACCGCTGTCTGCGTTCTCGGCAAAGCTTCACGTACCACACCATGTGTCGGAGCTTGGTAGACGGACATACCGACCGCGTCATTGGCGGATTCCATTACTCCCTAGCCCTATGCGTTATTGTAGTGGCACTTCGGAGACCCCCCGGCACTACTATGGTATCTGGTGTAA"
print(longestCommonSubsequence(text1, text2))