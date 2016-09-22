"""Wrapper for motif discovery tool DREME."""
import os

from subprocess import PIPE, Popen

from utilities import MotifWrapper

from Bio.Alphabet import IUPAC

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Dreme(MotifWrapper):
    """
    Wrapper for DREME 4.11.0.

    Usage: dreme [options] -p <primary sequence file> [-n <control sequence file>]
    """

    def __init__(self,
                 output_dir='dreme_out',

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
        # self.n_seqs = 0
        # to store the names of sequences as given in input file to fit()
        # self.seq_names = list()
        # over-rides same attribute of MotifWrapper class
        # self.pseudocounts = pseudocounts
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

        logger.info(stdout)

    def fit(self, fasta_file='', control_file=None):
        """Save the output of DREME and parse it."""
        if not fasta_file:
            return NameError('Input fasta file not specified')

        cmd_params = self._make_param_string()

        self._command_exec(fasta_file, control_file, cmd_params)

        filename = os.path.join(self.output_dir, 'dreme.txt')

        # record = self._parse_output(filename)
        # parsing
        record = Record()
        with open(filename) as handle:
            record._get_data(handle)
        return record


class Record(list):
    """A class for holding the results of a DREME run."""

    def __init__(self):
        """init."""
        self.version = ""
        self.alphabet = None
        # self.instances = []
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
        self.consensus_seqs = consensus[:]
        self.widths = lengths[:]
        self.nsites = num_sites[:]
        self.e_values = evalues[:]
        self.prob_matrices = matrices[:]

    def _get_data(self, handle):
        self._read_version(handle)
        self._read_alphabet(handle)
        self._read_motifs(handle)
