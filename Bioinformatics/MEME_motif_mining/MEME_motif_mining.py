import numpy as np
from sklearn.preprocessing import normalize
alphabet = "ATGC"
sequences = [
    "AGG", "CGT", "AGC"
]
motif_len = 2
pseudocount_constant = 1

def to_base(seq):
    return ''.join(["ATGC "[x] for x in seq])

def likelihood(x, f, is_motif):
    if is_motif == True:
        prob = [f[pos+1, c] for (pos, c) in enumerate(x)]
    else:
        prob = [f[0, c] for c in x]
    return np.prod(prob)

def full_log_likelihood(X, Z, lamda, f):
    log_likelihood = 0
    for i in range(len(X)):
        for j in (0, 1):
            log_likelihood += Z[i, j] * (np.log(likelihood(X[i], f, j)) + np.log(lamda[j]))
    return log_likelihood

kmer_split = sum([[seq[i:i+motif_len] for i in range(len(seq)-motif_len+1)] for seq in sequences], [])
alphabet_dict = {'A':0, 'T':1, 'G':2, 'C':3}
sequences = [list(map(lambda c : alphabet_dict[c], seq)) for seq in kmer_split]
num_alphabet = len(alphabet)
num_sequences = len(sequences)

# Initialize theta = (lambda, f)
lamda_old = np.float32((0.5, 0.5))
f_old = np.float32(
    [ 
        [.25, .25, .25, .25],
        [1/6, 1/6, 1/2, 1/6],
        [1/6, 1/6, 1/6, 1/2]
    ]
)

for iteration in range(10):
    # E Step : Compute expected value of latent variable
    z = np.zeros((num_sequences, 2))
    for i in range(num_sequences):
        z[i, 0] = likelihood(sequences[i], f_old, False) * lamda_old[0] # Background
        z[i, 1] = likelihood(sequences[i], f_old, True)  * lamda_old[1] # Motif
    z = normalize(z, norm='l1').round(2)

    # M step : Using z, compute MLE within model parameters lambda and f
    lamda_new = normalize(z.sum(axis=0, keepdims=True), norm='l1').reshape(-1, 1)
    c = np.zeros_like(f_old)
    for i in range(num_sequences):
        for j in range(motif_len):
            c[0, sequences[i][j]] += z[i, 0]
    for i in range(num_sequences):
        for j in range(motif_len):
            c[j+1, sequences[i][j]] += z[i, 1]
    f_new = np.zeros_like(f_old)
    for i in range(motif_len+1):
        for j in range(num_alphabet):
            f_new[i, j] = (c[i, j] + pseudocount_constant) / (c[i].sum() + pseudocount_constant * num_alphabet)