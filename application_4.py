"""
A set of functions + script for comparing two genetic sequences represented as strings,
as well as implementing a spelling correction based on the edit distance concept.
"""
import math
import random
import urllib.request

import matplotlib.pyplot as plt
import project_4 as helper_module


def read_scoring_matrix(filename):
    """
    Read a scoring matrix from the file named filename.  

    Argument:
    filename -- name of file containing a scoring matrix

    Returns:
    A dictionary of dictionaries mapping X and Y characters to scores
    """
    scoring_dict = {}
    scoring_file = urllib.request.urlopen(filename)
    ykeys = scoring_file.readline().decode('utf-8')
    ykeychars = ykeys.split()
    lines = scoring_file.readlines()
    for line in lines:
        line_dec = line.decode('utf-8')
        vals = line_dec.split()
        xkey = vals.pop(0)
        scoring_dict[xkey] = {}
        for ykey, val in zip(ykeychars, vals):
            scoring_dict[xkey][ykey] = int(val)
    return scoring_dict


def read_protein(filename):
    """
    Read a protein sequence from the file named filename.

    Arguments:
    filename -- name of file containing a protein sequence

    Returns:
    A string representing the protein
    """
    protein_file = urllib.request.urlopen(filename)
    protein_seq = str(protein_file.read(), "utf-8")
    protein_seq = protein_seq.rstrip()
    return protein_seq


def read_words(filename):
    """
    Load word list from the file named filename.

    Returns a list of strings.
    """
    # load assets
    word_file = urllib.request.urlopen(filename)
    
    # read in files as string
    words = str(word_file.read(), "utf-8")
    
    # template lines and solution lines list of line string
    word_list = words.split('\n')
    print ("Loaded a dictionary with", len(word_list), "words")
    return word_list


def generate_null_distribution(seq_x, seq_y, scoring_matrix, num_trials):
    """
    Computes the null distribution of scores between two sequences
    by shuffling one of them (num_trials) times and calculating the score
    of local alignment for each trial.

    Returns a dictinary with keys corresponding to score values and values
    to the number of occurencies of the particular score.
    """
    scoring_distribution = {}
    for _ in range(num_trials):
        rand_y = list(seq_y)
        random.shuffle(rand_y)
        align_m = helper_module.compute_alignment_matrix(seq_x, rand_y, scoring_matrix, False)
        score = helper_module.compute_local_alignment(seq_x, rand_y, scoring_matrix, align_m)[0]
        if score in scoring_distribution:
            scoring_distribution[score] += 1
        else:
            scoring_distribution[score] = 1
    
    return scoring_distribution


def plot_distribution(dist, num_trials):
    """
    Plots a bar chart for the normalized
    null distribution with (num_trials) 
    using two protein sequences.
    """
    # sorted by key, return a list of tuples
    lists = sorted(dist.items())

    # unpack a list of pairs into two tuples
    x_vals, y_vals = zip(*lists)

    # create normalized score values
    y_vals_norm = []
    for y_val in y_vals:
        y_vals_norm.append(y_val / num_trials)

    # plot the graph
    plt.bar(x_vals, y_vals_norm)

    # mark 'x' and 'y' axis accordingly
    plt.xlabel('Score')
    plt.ylabel('Fraction of total trials')        
    
    # assign a title
    plt.title('Normalized null distribution')

    plt.legend(loc = 'upper right', title = f'Number of trials:\n{NUM_TRIALS}')
    
    # exhibit the plot
    plt.show()


def check_spelling(checked_word, dist, word_list):
    """
    iterates through (word_list) and returns the set of all words 
    that are within edit distance (dist) of the string (checked_word).
    """
    DIAG_SCORE = 2 
    OFF_DIAG_SCORE = 1 
    DASH_SCORE = 0
    ALPHABET = set(["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", \
                    "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z"])
    scoring_matrix = helper_module.build_scoring_matrix(ALPHABET, DIAG_SCORE, OFF_DIAG_SCORE, DASH_SCORE)
    a_like = set()
    for word in word_list:
        alig_m = helper_module.compute_alignment_matrix(checked_word, word, scoring_matrix, True)
        glob_alig = helper_module.compute_global_alignment(checked_word, word, scoring_matrix, alig_m)
        ed_dist = len(checked_word) + len(word) - glob_alig[0]
        if ed_dist == dist:
            a_like.add(word)
    
    return a_like


if __name__ == "__main__":
    """
    Runs a script for solving problems from the Application #4
    """
    # URLs for data files
    PAM50_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_PAM50.txt"
    HUMAN_EYELESS_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_HumanEyelessProtein.txt"
    FRUITFLY_EYELESS_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_FruitflyEyelessProtein.txt"
    CONSENSUS_PAX_URL = "http://storage.googleapis.com/codeskulptor-alg/alg_ConsensusPAXDomain.txt"
    WORD_LIST_URL = "http://storage.googleapis.com/codeskulptor-assets/assets_scrabble_words3.txt"

    # loading the necessary external data sets
    scoring_matrix = read_scoring_matrix(PAM50_URL)
    human_protein = read_protein(HUMAN_EYELESS_URL)
    fly_protein = read_protein(FRUITFLY_EYELESS_URL)
    consensus = read_protein(CONSENSUS_PAX_URL)
    word_list = read_words(WORD_LIST_URL)

    # Question 1
    allignment_matrix = helper_module.compute_alignment_matrix(human_protein, fly_protein, scoring_matrix, False)
    score, human_local, fly_local = helper_module.compute_local_alignment(human_protein, fly_protein, scoring_matrix, allignment_matrix)
    print("Human local alignment:", human_local)
    print("Fly local alignment:", fly_local)
    print(len(human_local) == len(fly_local))
    print("Human - fly local alignment score:", score)

    # Question 2
    human_local = human_local.replace("-", "")
    fly_local = fly_local.replace("-", "")

    allignment_matrix_human = helper_module.compute_alignment_matrix(human_local, consensus, scoring_matrix, True)
    human_consensus = helper_module.compute_global_alignment(human_local, consensus, scoring_matrix, allignment_matrix_human)
    allignment_matrix_fly = helper_module.compute_alignment_matrix(fly_local, consensus, scoring_matrix, True)
    fly_consensus = helper_module.compute_global_alignment(fly_local, consensus, scoring_matrix, allignment_matrix_fly)

    counter = 0
    for idx in range(len(human_consensus[1])):
        if human_consensus[1][idx] == human_consensus[2][idx]:
            counter += 1
    human_perc = counter / len(human_consensus[1]) * 100
    print (f"Human vs consensus: {human_perc}%")

    counter = 0
    for idx in range(len(fly_consensus[1])):
        if fly_consensus[1][idx] == fly_consensus[2][idx]:
            counter += 1
    fly_perc = counter / len(fly_consensus[1]) * 100
    print (f"Fly vs consensus: {fly_perc}%")

    # Question 3
    amino_acids = list("ACBEDGFIHKMLNQPSRTWVYXZ")
    rand_seq = ""
    for _ in range(len(consensus)):
        rand_gene = random.choice(amino_acids)
        rand_seq += rand_gene
    
    counter = 0
    for idx in range(len(consensus)):
        if rand_seq[idx] == consensus[idx]:
            counter += 1
    rand_perc = counter / len(consensus) * 100
    print (f"Random vs consensus: {rand_perc}%")

    # Question 4
    NUM_TRIALS = 1000
    null_dist = generate_null_distribution(human_protein, fly_protein, scoring_matrix, NUM_TRIALS)
    plot_distribution(null_dist, NUM_TRIALS)

    # Question 5
    mean = sum([key * val for key, val in null_dist.items()]) / NUM_TRIALS
    st_dev = math.sqrt(sum([((key - mean) ** 2) * val for key, val in null_dist.items()]) / NUM_TRIALS)
    z_score = (score - mean) / st_dev
    print("Mean of the distribution:", mean)
    print("Standard deviation of the distribution:", st_dev)
    print("Z-score of a human-fly local alignment:", z_score)

    # Question 8
    humble_1 = check_spelling("humble", 1, word_list)
    firefly_2 = check_spelling("firefly", 2, word_list)
    print("Humble - 1 editing distance:", humble_1)
    print("Firefly - 2 editing distances:", firefly_2)

