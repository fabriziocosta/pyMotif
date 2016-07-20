"""A wrapper of the motif discovery tool GLAM2."""
from subprocess import PIPE, Popen

import commands

from utilities import MotifWrapper

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Glam2(MotifWrapper):
    """Python wrapper for glam2 v=4.11.0 motif discovery tool.

    Original usage: Usage: glam2 [options] alphabet my_seqs.fa

    Attributes:
        alphabet: "p" for proteins, "n" for nucleotides (n)
        fasta_file: Input file with sequences in FASTA format
        outptu_dir: Directory to save the output (glam2_out); -O TODO: check path compatibility
        number_alignment_runs: number of alignment runs (10); -r
        number_iterations: end each run after this many iterations without improvement (10000); -n
        both_strands: examine both strands - forward and reverse complement; -2
        min_sequences: minimum number of sequences in the alignment (2); -z
        min_aligned_columns: minimum number of aligned columns (2); -a
        max_aligned_columns: maximum number of aligned columns (50); -b
        initial_aligned_columns: initial number of aligned columns (20); -w
        threshold:
    """

    def __init__(self,
                 alphabet="dna",    # ["dna", protein"]
                 gap_in_alphabet=True,
                 scoring_criteria="pwm",    # ["pwm", "hmm"]
                 output_dir="glam2_out",
                 number_alignment_runs=10,
                 number_iterations=10000,
                 both_strands=False,
                 min_sequences=2,
                 min_aligned_columns=2,
                 max_aligned_columns=50,
                 initial_aligned_columns=20,

                 muscle_obj=None,
                 weblogo_obj=None
                 ):
        """Return a Glam2 object with specified attribute values."""
        self.alphabet = alphabet
        self.gap_in_alphabet = gap_in_alphabet
        self.scoring_criteria = scoring_criteria
        # threshold for scoring sequences
        if scoring_criteria == 'pwm':
            self.threshold = 1.0e-9
        else:
            self.threshold = -200    # TODO: examine threshold for hmm
        self.output_dir = output_dir    # -O
        self.number_alignment_runs = number_alignment_runs    # -r
        self.number_iterations = number_iterations     # -n
        self.both_strands = both_strands    # -2
        self.min_sequences = min_sequences    # -z
        self.min_aligned_columns = min_aligned_columns    # -a
        self.max_aligned_columns = max_aligned_columns    # -b
        self.initial_aligned_columns = initial_aligned_columns    # -w

        self.muscle_obj = muscle_obj
        self.weblogo_obj = weblogo_obj

        # no. of seqs in input file, to be set by fit()
        self.n_seqs = 0
        # names of sequences as given in input file to fit()
        self.seq_names = list()
        # list-of-strings representation of motifs
        self.motives_list = list()
        # aligned list-of-strings of motifs
        self.aligned_motives_list = list()
        # list of sequence logos created with WebLogo
        self.logos = list()
        # threshold for scoring sequences
        self.threshold

    def _make_param_string(self):
        # Creates a string of parameters
        params = " -O " + self.output_dir

        if self.number_alignment_runs != 10:
            params += " -r " + str(self.number_alignment_runs)

        if self.number_iterations != 10000:
            params += " -n " + str(self.number_iterations)

        if self.both_strands is True:
            params += " -2"

        if self.min_sequences != 2:
            params += " -z " + str(self.min_sequences)

        if self.min_aligned_columns != 2:
            params += " -a " + str(self.min_aligned_columns)

        if self.max_aligned_columns != 50:
            params += " -b" + str(self.max_aligned_columns)

        if self.initial_aligned_columns != 20:
            params += " -w " + str(self.initial_aligned_columns)

        return params

    def _command_exec(self, fasta_file, params):
        if self.alphabet == "dna":
            alpha = ' n '
        elif self.alphabet == "rna":
            raise ValueError('Glam2 does not support RNA alphabet, use "dna" or "protein"')
        else:
            alpha = ' p '

        cmd = "glam2" + params + alpha + fasta_file
        io = Popen(cmd.split(" "), stdout=PIPE, stderr=PIPE)
        (stderr, stdout) = io.communicate()

        logger.info(stdout)

    def fit(self, fasta_file=''):
        """Save the output of Glam2 and parse it."""
        if not fasta_file:
            return NameError('Input fasta file not specified')

        cmd_params = self._make_param_string()
        self._command_exec(fasta_file, cmd_params)

        """
        command = "glam2" + " -O " + self.output_dir + " -r " + str(self.number_alignment_runs) + " -n " + str(self.number_iterations)

        if self.both_strands == 1:
            command = command + " -2 "

        command += " -z " + str(self.min_sequences) + " -a " + str(self.min_aligned_columns) + " -b " + str(self.max_aligned_columns) + " -w " + str(self.initial_aligned_columns) + " " + self.alphabet + " " + fasta_file

        status, output = commands.getstatusoutput(command)
        if status != 0:
            print "Error: Command not executed on terminal."
            return output

        return output
        """

    def predict(self, return_list=False):
        """Output 1 integer per sequence indicating the motif id present in the sequence."""
        pass

    def transform(self, return_match=False):
        """Transform."""
        pass

    def fit_predict(self, fasta_file, return_list=False):
        """Run fit and predict."""
        self.fit(fasta_file)
        return self.predict(return_list=return_list)

    def display_glam2_help(self):
        """Display GLAM2 command line help."""
        command = "glam2 -h"
        print commands.getoutput(command)
