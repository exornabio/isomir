COMPARISON = {"A": {"A"},
              "C": {"C"},
              "G": {"G"},
              "T": {"T"},
              "R": {"A", "G"},
              "Y": {"C", "T"},
              "S": {"C", "G"},
              "W": {"A", "T"},
              "K": {"G", "T"},
              "M": {"A", "C"},
              "B": {"C", "G", "T"},
              "D": {"A", "G", "T"},
              "H": {"A", "C", "T"},
              "V": {"A", "C", "G"},
              "N": {"A", "C", "G", "T"}}

AMBIGUOUS_LETTERS = {"N", "R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"}


class Isoform:
    def __init__(self, mirna_id: str, read_id: str, read_seq: str, read_num: int, dist: int):
        self.mirna_id = mirna_id
        self.read_id = read_id
        self.read_seq = read_seq
        self.read_num = read_num
        self.dist = dist


class Mirna:
    def __init__(self, id: str, motif: str, consensus: str):
        self.id = id
        self.motif = motif
        self.consensus = consensus


def edit_dist(str1: str, str2: str):
    """A Dynamic Programming to find minimum number operations to convert str1 to str2
    """
    m = len(str1)
    n = len(str2)
    # Create a table to store results of subproblems
    dp = [[0 for x in range(n + 1)] for x in range(m + 1)]

    # Fill d[][] in bottom up manner
    for i in range(m + 1):
        for j in range(n + 1):

            # If first string is empty, only option is to
            # insert all characters of second string
            if i == 0:
                dp[i][j] = j  # Min. operations = j

            # If second string is empty, only option is to
            # remove all characters of second string
            elif j == 0:
                dp[i][j] = i  # Min. operations = i

            # If last characters are same, ignore last char
            # and recur for remaining string
            elif str1[i - 1] == str2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1]

            # If last character are different, consider all
            # possibilities and find minimum
            else:
                dp[i][j] = 1 + min(dp[i][j - 1],  # Insert
                                   dp[i - 1][j],  # Remove
                                   dp[i - 1][j - 1])  # Replace

    return dp[m][n]


def compare_mirna_seq_nt(mirna_nt: str, seq_nt: str, ) -> bool:
    """Compare the mirna nt and read nt using COMPARISON
    
    Parameters
    ----------
    mirna_nt: str, mirna sequence nt
    seq_nt: str, read sequence nt

    Returns
    -------
    true mean same, false is diffrent
    """
    codes = COMPARISON[seq_nt]
    return mirna_nt in codes


def compare_mirna_seq(mirna: str, seq: str) -> bool:
    """Compare the mirna and read sequence with same length

    Parameters
    ----------
    mirna: str, mirna sequence
    seq: str, read sequence

    Returns
    -------
    true is same, unwise false
    """
    if len(seq) != len(mirna):
        print("length of mirna and seq is not equal!")
        exit(1)

    for i in range(len(seq)):
        if not compare_mirna_seq_nt(mirna[i], seq[i]):
            return False
    return True


def find_motif_in_seq(motif: str, seq: str, is_amb: bool) -> int:
    """Find the positon of motif on seq

    Parameters
    ----------
    motif: str
    seq: str

    Returns
    ----------
    position of motif on seq, -1 not found
    """
    if is_amb:
        motif_len = len(motif)
        for i in range(len(seq) - motif_len + 1):
            if compare_mirna_seq(motif, seq[i:motif_len + i]):
                return i
        return -1
    else:
        return seq.find(motif)


def seq_contain_n(seq: str) -> bool:
    """Check if read sequence contain N

    Returns
    ----------
    
    Returns
    true or false
    """
    return seq.find("N") != -1
