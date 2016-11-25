"""Wrapper for motif discovery tool DREME."""
import os

from time import time

from subprocess import PIPE, Popen

from utilities import MotifWrapper

from Bio.Alphabet import IUPAC
from Bio import SeqIO

import logging
logger = logging.getLogger(__name__)


class Dreme(MotifWrapper):
    """
    Wrapper for DREME 4.11.0.

    Usage: dreme [options] -p <primary sequence file> [-n <control sequence file>]
    """

    def __init__(self,
                 output_dir='dreme_out',
                 e_threshold=None,    # motifs with E-value less than this threshold will be recorded

                 # alphabet
                 alphabet='dna',    # ['dna','rna','protein']
                 gap_in_alphabet=True,
                 scoring_criteria='pwm',
                 threshold=None,
                 k=1,    # return top k-scores for hmm

                 # general
                 norc=True,
                 ngen=100,
                 seed=1,

                 # stopping conditions
                 e_value=0.05,
                 m_motifs=None,
                 t_seconds=None,

                 # set core motif width
                 mink=3,
                 maxk=7,
                 k_width=None,

                 verbosity=2,

                 # parameters for (parent class) MotifWrapper)
                 pseudocounts=0,
                 muscle_obj=None,
                 weblogo_obj=None
                 ):
        """Initialize a Dreme Object."""
        self.output_dir = output_dir
        self.e_threshold = e_threshold
        self.alphabet = alphabet
        self.gap_in_alphabet = gap_in_alphabet
        self.scoring_criteria = scoring_criteria
        if threshold is None:
            if scoring_criteria == 'pwm':
                self.threshold = 1.0e-9
            else:
                self.threshold = 0.8
        else:
            self.threshold = threshold
        self.k = k

        self.norc = norc
        self.ngen = ngen
        self.seed = seed
        self.e_value = e_value
        self.m_motifs = m_motifs
        self.t_seconds = t_seconds
        self.mink = mink
        self.maxk = maxk
        self.k_width = k_width
        self.verbosity = verbosity

        self.muscle_obj = muscle_obj
        self.weblogo_obj = weblogo_obj

        # no. of seqs in input file, to be set by fit()
        self.n_seqs = 0
        # to store the names of sequences as given in input file to fit()
        self.seq_names = list()
        # over-rides same attribute of MotifWrapper class
        self.pseudocounts = pseudocounts
        # number of motives found
        self.nmotifs = None
        # list-of-strings representation of motifs
        self.motives_list = list()
        # aligned list-of-strings of motifs
        self.aligned_motives_list = list()
        # list of sequence logos created with WebLogo
        self.logos = list()
        # list of PWMS created by DREME
        self.pwms_list = list()
        # list of motif widths
        self.widths = list()
        # list of motif consensus sequences
        self.consensus = list()
        # record of motif data
        self.record = None
        # list of sequence logos created with Weblogo
        self.logos = list()

    def _make_param_string(self):
        # creates a string of parameters
        params = '-oc ' + self.output_dir

        if self.alphabet == 'dna':
            params += ' -dna'
        elif self.alphabet == 'rna':
            params += ' -rna'
        else:
            params += ' -protein'

        if self.norc is True:
            params += ' -norc'

        if self.ngen != 100:
            params += ' -g ' + str(self.ngen)

        if self.seed != 1:
            params += ' -s ' + str(self.seed)

        if self.e_value != 0.05:
            params += ' -e ' + str(self.e_value)

        if self.m_motifs is not None:
            params += ' -m ' + str(self.m_motifs)

        if self.t_seconds is not None:
            params += ' -t ' + str(self.t_seconds)

        if self.mink != 3:
            params += ' -mink ' + str(self.mink)

        if self.maxk != 7:
            params += ' -maxk ' + str(self.maxk)

        if self.k_width is not None:
            params += ' -k ' + str(self.k_width)

        if self.verbosity != 2:
            params += ' -verbosity ' + str(self.verbosity)

        return params

    def _command_exec(self, primary_file, control_file, params):
        cmd = 'dreme ' + params + ' -p ' + str(primary_file)

        if control_file is not None:
            cmd += cmd + ' -n ' + str(control_file)

        io = Popen(cmd.split(" "), stdout=PIPE, stderr=PIPE)
        (stderr, stdout) = io.communicate()

    def fit(self, fasta_file='', control_file=None):
        """Save the output of DREME and parse it."""
        start = time()
        if not fasta_file:
            return NameError('Input fasta file not specified')
        cmd_params = self._make_param_string()
        self._command_exec(fasta_file, control_file, cmd_params)
        end = time()
        logger.debug('DREME finished in %d s' % (end - start))

        start = time()
        filename = os.path.join(self.output_dir, 'dreme.txt')
        # record is useful information retrieved from parsing output
        record = Record(threshold=self.e_threshold)
        with open(filename) as handle:
            record._get_data(handle)

        self.record = record
        self.consensus = record.consensus_seqs[:]
        self.widths = record.widths[:]
        self.nmotifs = len(record.consensus_seqs)
        headers, seqs = self._parse_fasta(fasta_file)
        self.nseqs = len(headers)
        self.seq_names = headers[:]
        self.motives_list = self._get_motives_list(self.nmotifs, headers, seqs)
        super(Dreme, self).fit(motives=self.motives_list)
        end = time()
        logger.debug('Processing DREME output finished in %d s' % (end - start))

    def _parse_fasta(self, filename):
        headers = []
        seqs = []
        for seq_record in SeqIO.parse(filename, "fasta"):
            headers.append(seq_record.id)
            seqs.append(str(seq_record.seq))
        return headers, seqs

    def _get_motif_seqs(self, scores, headers, seqs, width):
        motif_occs = []
        for scr in scores:
            i, j, k = scr
            header_x = headers[k]
            seq_x = seqs[k]
            motif = seq_x[j:j + width]
            motif_occs.append((header_x, motif))
        return motif_occs

    def _get_scores_list(self, seqs, width, alpha, pm):
        scores = []    # best segment score for each sequence
        for k, seq in enumerate(seqs):
            seq_scores = []
            for i in xrange(len(seq) - width):
                segment_scr = 0
                for j in xrange(width):
                    try:
                        alpha_id = alpha.index(seq[i + j])
                        segment_scr += pm[j][alpha_id]

                    except ValueError:
                        # ignore any alphabet not in the alphabet list
                        continue
                seq_scores.append((segment_scr, i))
            max_segment = max(seq_scores)
            data = (max_segment[0], max_segment[1], k)    # (score, start_index, seq_id)
            scores.append(data)
        return scores

    def _modify_chars(self, motives_list):
        mod_motif_list = []
        for motif_i in motives_list:
            heads, seq = [list(x) for x in zip(*motif_i)]
            seq_new = []
            for s in seq:
                s_list = list(s)
                for i, char in enumerate(s_list):
                    if char not in self.record.alphabet.letters:
                        s_list[i] = "-"
                s_new = "".join(s_list)
                seq_new.append(s_new)
            mod_motif_list.append(zip(heads, seq_new))
        return mod_motif_list

    def _get_motives_list(self, nmotifs, headers, seqs):
        motives_list = []
        for n in xrange(nmotifs):
            width = self.record.widths[n]
            nsites = self.record.nsites[n]
            pm = self.record.prob_matrices[n]
            alpha = self.record.alphabet.letters
            alpha = sorted(alpha)

            scores = self._get_scores_list(seqs, width, alpha, pm)

            scores_sort = sorted(scores, reverse=True)
            scores_top = scores_sort[:nsites]

            motif_data = self._get_motif_seqs(scores_top, headers, seqs, width)
            motives_list.append(motif_data)

        motives_list = self._modify_chars(motives_list)
        return motives_list


class Record(list):
    """A class for holding the results of a DREME run."""

    def __init__(self, threshold=None):
        """init."""
        self.threshold = threshold

        self.version = ""
        self.alphabet = None
        self.consensus_seqs = list()
        self.widths = list()
        self.nsites = list()
        self.e_values = list()
        self.prob_matrices = list()

    def _read_version(self, handle):
        for line in handle:
            if line.startswith('# DREME'):
                break
        else:
            raise ValueError("Improper input file. File should contain a line starting DREME.")
        line = line.strip()
        ls = line.split()
        self.version = ls[2]

    def _read_alphabet(self, handle):
        for line in handle:
            if line.startswith('ALPHABET'):
                break
        if not line.startswith('ALPHABET'):
            raise ValueError("Line does not start with 'ALPHABET':\n%s" % line)

        line = line.strip()

        if 'DNA' in line:
            al = IUPAC.unambiguous_dna
        elif 'RNA' in line:
            al = IUPAC.unambiguous_rna
        else:
            al = IUPAC.protein

        self.alphabet = al

    def _read_motifs(self, handle):
        # insts = []
        consensus = []
        lengths = []
        num_sites = []
        evalues = []
        matrices = []
        for line in handle:
            if '# Stopping reason:' in line:
                break

            if line.startswith('MOTIF'):
                # instance = Instance()
                consensus.append(line.split(' ')[1])

            if line.startswith('letter-probability matrix'):
                line = line.split()
                width = int(line[5])
                lengths.append(width)
                num_sites.append(int(line[7]))
                evalues.append(float(line[9]))

                pwm = []
                for i in range(width):
                    line = next(handle)
                    data = line.strip()
                    data = data.split(' ')
                    data = [float(x) for x in data]
                    pwm.append(data)

                matrices.append(pwm)

        # select motives below E-value threshold
        if self.threshold:
            indexes = [i for i, x in enumerate(evalues) if x <= self.threshold]
            consensus = [consensus[ind] for ind in indexes]
            lengths = [lengths[ind] for ind in indexes]
            num_sites = [num_sites[ind] for ind in indexes]
            evalues = [evalues[ind] for ind in indexes]
            matrices = [matrices[ind] for ind in indexes]

        self.consensus_seqs = consensus[:]
        self.widths = lengths[:]
        self.nsites = num_sites[:]
        self.e_values = evalues[:]
        self.prob_matrices = matrices[:]

    def _get_data(self, handle):
        self._read_version(handle)
        self._read_alphabet(handle)
        self._read_motifs(handle)
