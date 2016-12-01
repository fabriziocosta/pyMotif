"""Script for finding optimal SMoD parameters for a given noise level."""
import sys


# In[1]:
from smod_wrapper import SMoDWrapper
from sklearn.cluster import KMeans
import numpy as np
from sklearn.metrics import roc_auc_score
import datetime
import time

# In[2]:

from eden.util import configure_logging
import logging


# In[3]:
from dataset_generator import make_artificial_dataset

noise_level = float(sys.argv[1])

# In[4]:


def score_seqs(seqs, n_motives, tool):
    """Score every sequences in the given test data set."""
    scores = []
    if tool is None:
        return scores

    for j in range(len(seqs)):
        seq_scr = []
        iters = tool.nmotifs
        for k in range(iters):
            scr = tool.score(motif_num=k + 1, seq=seqs[j][1])
            seq_scr.append(scr)

        # taking average over all motives for a sequence
        if len(seq_scr) > 1:
            x = np.array(seq_scr[0])
            for l in range(1, iters):
                x = np.vstack((x, seq_scr[l]))
            seq_scr = list(np.mean(x, axis=0))
            scores.append(seq_scr)
        elif len(seq_scr) == 1:
            scores.append(np.array(seq_scr[0]))
        else:
            raise ValueError("no sequence score")
    return scores


# In[5]:

def get_dataset(sequence_length=200,
                n_sequences=200,
                motif_length=10,
                n_motives=2,
                p=0.2,
                random_state=1):
    """Generate, preprocess and return the dataset."""
    motives, pos_seqs, binary_seq = make_artificial_dataset(alphabet='ACGT',
                                                            sequence_length=sequence_length,
                                                            n_sequences=n_sequences,
                                                            motif_length=motif_length,
                                                            n_motives=n_motives,
                                                            p=p,
                                                            random_state=random_state)

    from eden.modifier.seq import seq_to_seq, shuffle_modifier
    neg_seqs = seq_to_seq(
        pos_seqs, modifier=shuffle_modifier, times=2, order=2)
    neg_seqs = list(neg_seqs)

    block_size = n_sequences / 8

    pos_size = len(pos_seqs)
    train_pos_seqs = pos_seqs[:pos_size / 2]
    test_pos_seqs = pos_seqs[pos_size / 2:]

    neg_size = len(neg_seqs)
    train_neg_seqs = neg_seqs[:neg_size / 2]
    # test_neg_seqs = neg_seqs[neg_size / 2:]

    true_score = [float(int(i)) for i in binary_seq]
    return (block_size, train_pos_seqs, train_neg_seqs, test_pos_seqs, n_motives, true_score)


# In[6]:
def test_on_datasets(n_sets=5, param_setting=None, p=0.2, max_roc=0.5, std_roc=0.01):
    """Test the parameter setting on a given number of datasets."""
    dataset_score = []
    seeds = [i * 2000 for i in range(1, n_sets + 1)]
    for k in range(n_sets):
        print "n_sets: %d" % (k + 1)
        # Generate data set
        seed = seeds[k]
        data = get_dataset(sequence_length=300, n_sequences=1000, motif_length=10, n_motives=4, p=p, random_state=seed)
        block_size = data[0]
        train_pos_seqs = data[1]
        train_neg_seqs = data[2]
        test_pos_seqs = data[3]
        n_motives = data[4]
        true_score = data[5]

        smod = SMoDWrapper(alphabet='dna',
                           scoring_criteria='pwm',

                           complexity=5,
                           n_clusters=10,
                           min_subarray_size=8,
                           max_subarray_size=12,
                           clusterer=KMeans(),
                           pos_block_size=block_size,
                           neg_block_size=block_size,
                           sample_size=300,
                           p_value=param_setting['p_value'],
                           similarity_th=param_setting['similarity_th'],
                           min_score=param_setting['min_score'],
                           min_freq=param_setting['min_freq'],
                           min_cluster_size=param_setting['min_cluster_size'],
                           regex_th=param_setting['regex_th'],
                           freq_th=param_setting['freq_th'],
                           std_th=param_setting['std_th'])

        smod.fit(train_pos_seqs, train_neg_seqs)

        try:
            scores = score_seqs(seqs=test_pos_seqs,
                                n_motives=n_motives,
                                tool=smod)
        except:
            continue

        mean_score = np.mean(scores, axis=0)
        roc_score = roc_auc_score(true_score, mean_score)

        """
        # if a parameter setting performs poorly, don't test on other datasets
        z-score = (x - mu)/sigma
        if ((roc_score - max_roc)/std_roc) > 2:"""
        if roc_score < 0.6:
            # print "discarding parameter setting..."
            break

        dataset_score.append(roc_score)
    return dataset_score


# In[ ]:

def check_validity(key, value, noise):
    """Check if the generated value for a particular parameter is fit for use."""
    if key == 'min_score':    # atleast greater than (motif_length)/2
        if value >= 5:
            return True, int(round(value))
    elif key == 'min_cluster_size':
        if value >= 3:
            return True, int(round(value))
    elif key == 'min_freq':    # atmost (1 - noise_level)
        if value > 0 and value <= (1 - noise):
            return True, value
    elif key == 'p_value':
        if value <= 1.0 and value >= 0.0:
            return True, value
    elif key == 'similarity_th':
        if value <= 1.0 and value >= 0.8:
            return True, value
    elif key == 'regex_th':
        if value > 0 and value <= 0.3:
            return True, value
    elif key == 'freq_th':
        if value <= 1.0 and value > 0:
            return True, value
    elif key == 'std_th':
        if value <= 1.0 and value > 0:
            return True, value
    else:
        raise ValueError('Invalid key: ', key)
    return False, value


def random_setting(parameters=None, best_config=None, noise=None):
    """Generate a random parameter setting."""
    parameter_setting = {}
    max_iter = 1000
    # use best_configuration of last run as initial setting
    if not parameters['min_score']:
        for key in parameters.keys():
            parameters[key].append(best_config[key])
            parameter_setting[key] = best_config[key]
    else:
        for key in parameters.keys():
            success = False
            n_iter = 0
            mu = np.mean(parameters[key])
            sigma = np.std(parameters[key])
            if sigma == 0:
                sigma = 0.01
            while not success:
                # if max_iterations exceeded, return mean as value
                if n_iter == max_iter:
                    value = mu
                    if key in ['min_score', 'min_cluster_size']:
                        value = int(round(value))
                    # parameter_setting[key] = value
                    break
                value = np.random.normal(mu, 2 * sigma)
                n_iter += 1
                success, value = check_validity(key, value, noise)
            parameter_setting[key] = value
    return parameter_setting


# In[ ]:

# %%time

logger = logging.getLogger()
configure_logging(logger, verbosity=1)

filename = "Result_at_" + str(noise_level) + ".txt"

best_config = {'min_score': 6,  # atleast motif_length/2
               'min_freq': 0.1,  # can not be more than (1- noise level)
               'min_cluster_size': 3,  # atleast 3
               'p_value': 0.1,  # atleast 0.1
               'similarity_th': 0.8,  # 0.8
               'regex_th': 0.3,  # max 0.3
               'freq_th': 0.05,  # 0.05
               'std_th': 0.2}  # 0.2


results_dic = {}

REPS = 120    # different settings to be tried

# for i in param:
parameters = {'min_freq': [],
              'min_cluster_size': [],
              'p_value': [],
              'similarity_th': [],
              'min_score': [],
              'regex_th': [],
              'freq_th': [],
              'std_th': []}
max_roc = 0.5
std_roc = 0.01
# parameters = generate_dist(parameters, best_config)

with open(filename, "w") as f:
    f.write(datetime.datetime.fromtimestamp(time.time()).strftime('%H:%M:%S') + " Starting experiment...\n")

    for j in range(REPS):
        print
        print "REPS: %d" % (j + 1)
        # i)    # Randomize Parameter setting
        param_setting = random_setting(parameters, best_config, noise_level)
        N_SETS = 5    # Different data sets
        dataset_score = test_on_datasets(n_sets=N_SETS,
                                         param_setting=param_setting,
                                         p=noise_level,
                                         max_roc=max_roc,
                                         std_roc=std_roc)
        if not dataset_score:
            continue
        mean_roc = np.mean(dataset_score)
        std = np.std(dataset_score)

        if mean_roc > max_roc:
            max_roc = mean_roc
            std_roc = std
            current_time = datetime.datetime.fromtimestamp(time.time()).strftime('%H:%M:%S')
            f.write(current_time + " Better Configuration found at perturbation prob = " + str(noise_level) + '\n')
            f.write("ROC: " + str(mean_roc) + '\n')
            f.write("Parameter Configuration: " + str(param_setting) + '\n\n')

            best_config = param_setting
            param_setting["ROC"] = mean_roc
            results_dic[noise_level] = param_setting

            if mean_roc > 0.97:
                break
    current_time = datetime.datetime.fromtimestamp(time.time()).strftime('%H:%M:%S')
    f.write(current_time + " Experiment finished.\n")
