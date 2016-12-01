"""Dataset generation module."""
import random
import numpy as np


def random_string(length, alphabet_list):
    """Generate a random string."""
    rand_str = ''.join(random.choice(alphabet_list) for i in range(length))
    return rand_str


def perturb(seed, alphabet_list, p=0.5):
    """Randomize a string."""
    seq = ''
    for c in seed:
        if random.random() < p:
            c = random.choice(alphabet_list)
        seq += c
    return seq


def inflate_normalize(pwms=None, exp=2):
    """Inflate and normalize PWM to utilize noise level."""
    num_motives = pwms.shape[0]

    for j in range(num_motives):    # inflate
        pwms[j] = pwms[j] ** exp

    for j in range(num_motives):    # normalize
        pwms[j] = pwms[j] / pwms[j].sum(axis=0)

    return pwms


def get_pwms(alphabet='ACGT', num=2, length=6, exp=2):
    """Generate PWMs for every motif."""
    letters = len(alphabet)
    pwms = []
    for i in range(num):
        i_pwm = np.random.random_sample((letters, length))
        i_pwm = i_pwm / i_pwm.sum(axis=0)    # normalize
        pwms.append(i_pwm)
    pwms = np.array(pwms)
    pwms = inflate_normalize(pwms=pwms, exp=exp)

    return pwms


def motif_from_pwm(alphabet_list, pwm):
    """Create motif string from the PWM."""
    seq = ""
    length = pwm.shape[1]

    for i in range(length):
        alphabet_dist = pwm[:, i]
        c = np.random.choice(a=alphabet_list, p=alphabet_dist)
        seq += c
    return seq


def make_artificial_dataset(alphabet='ACGT', motif_length=6, sequence_length=100,
                            n_sequences=1000, n_motives=2, p=0.2, random_state=1):
    """Generate artificial dataset.

    Returns: motives - list of motives used in sequences
             seq - dataset as list of sequences
             binary_seq - a sequence of 0's & 1's which can be used for computing ROC score.
    """
    random.seed(random_state)
    alphabet_list = [c for c in alphabet]

    pwms = get_pwms(alphabet=alphabet, num=n_motives, length=motif_length, exp=2 - p)

    sequence_length = sequence_length / n_motives
    flanking_length = (sequence_length - motif_length) / 2
    n_seq_per_motif = n_sequences

    counter = 0
    seqs = []
    motives = []
    for i in range(n_seq_per_motif):
        total_seq = ''
        motif = []
        for j in range(n_motives):
            left_flanking = random_string(flanking_length, alphabet_list)
            right_flanking = random_string(flanking_length, alphabet_list)
            noisy_motif = motif_from_pwm(alphabet_list, pwms[j])
            seq = left_flanking + noisy_motif + right_flanking
            total_seq += seq
            motif.append(noisy_motif)
        seqs.append(('ID%d' % counter, total_seq))
        motives.append(motif)
        counter += 1
    binary_skeleton = '0' * flanking_length + \
        '1' * motif_length + '0' * flanking_length
    binary_seq = binary_skeleton * n_motives
    return motives, seqs, binary_seq
