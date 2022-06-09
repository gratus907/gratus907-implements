import numpy as np
from sklearn.preprocessing import normalize

sequences = [
    "AATCCACGTCGT",
    "ATCATGCAGCTA",
    "TCTCGATTCCGA",
    "CTAACCGATGAC",
    "AGATCCTGCACT",
    "CAGAATCCTGTC"
]
motif_len = 4

def sequence_kmer_split(seq): 
    starting_point = min([seq.index(i) for i in range(4)])
    return [seq[i:i+motif_len] for i in range(starting_point, starting_point + sequence_length - motif_len+1)], starting_point

f = np.float32(
    [ 
        [.25, .25, .25, .25],
        [ 0.4,  0.3 , 0.1 , 0.2 ],
        [ 0.3,  0.5 , 0.1 , 0.1 ],
        [ 0.1,  0.2 , 0.4 , 0.3 ],
        [ 0.2,  0.1 , 0.6 , 0.1 ]
    ]
)

def likelihood(x, f, is_motif):
    if is_motif == True:
        prob = [f[pos+1, c] for (pos, c) in enumerate(x)]
    else:
        prob = [f[0, c] for c in x]
    return np.prod(prob)

def to_base(seq):
    return ''.join(["ATCG "[x] for x in seq])

def highlight_motif(base_sequence, motif_position):
    return base_sequence[:motif_position] + '\033[4m\033[91m' + base_sequence[motif_position:motif_position+motif_len] + '\033[0m' + base_sequence[motif_position+motif_len:]

PADDING = 4
sequence_length = len(sequences[0])
motif_position = sequence_length
alphabet_dict = {'X': PADDING, 'A':0, 'T':1, 'C':2, 'G':3}
sequences = [list(map(lambda c : alphabet_dict[c], 'X'*sequence_length+s)) for s in sequences]
pseudocount_constant = 0.1
ALPHABETS = [0, 1, 2, 3]

for iteration in range(50):
    current_idx = iteration % len(sequences)
    kmers, starting_point = sequence_kmer_split(sequences[current_idx])
    # 0 = non-motif, 1 = motif

    # Compute param f
    cnt = np.zeros_like(f)
    for idx, seq in enumerate(sequences):
        if idx == current_idx: continue
        for (i, s) in enumerate(seq):
            if s == PADDING: continue
            if motif_position <= i and i < motif_position + motif_len:
                cnt[i-motif_position+1, s] += 1
            else: cnt[0, s] += 1
    for i in range(motif_len):
        for j in ALPHABETS:
            f[i+1, j] = (cnt[i+1, j] + pseudocount_constant) / (len(sequences) - 1 + len(ALPHABETS)*pseudocount_constant)
    f[0, :] = (cnt[0, :] + pseudocount_constant) / ((len(sequences)-1)*(sequence_length - motif_len) + len(ALPHABETS)*pseudocount_constant)
    
    l = np.zeros((len(kmers), 2))
    for (i, s) in enumerate(kmers):
        for j in (0, 1):
            l[i, j] = likelihood(s, f, j)
    position_score = np.squeeze(normalize(np.expand_dims(l[:, 1] / l[:, 0], axis=0), norm='l1'))
    reposition = np.random.choice(np.array(list(range(len(position_score)))), p = position_score)
    sequences[current_idx] = [PADDING]*(sequence_length - reposition) + sequences[current_idx][starting_point:starting_point+sequence_length] + [PADDING]*(reposition)
    motif_segment = sequences[current_idx][motif_position:motif_position+motif_len]

    print(highlight_motif(to_base(sequences[current_idx]), motif_position), to_base(motif_segment))
    if iteration % 6 == 5 : print()