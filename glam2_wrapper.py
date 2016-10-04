"""A wrapper of the motif discovery tool GLAM2."""
import os

from subprocess import PIPE, Popen

import commands

from utilities import MotifWrapper

from Bio.Alphabet import IUPAC

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Glam2(MotifWrapper):
    """Python wrapper for glam2 v=4.11.0 motif discovery tool.

    Original usage: Usage: glam2 [options] alphabet my_seqs.fa
    """

    def __init__(self,
                 alphabet="dna",    # ["dna", protein"]
                 gap_in_alphabet=True,
                 scoring_criteria="pwm",    # ["pwm", "hmm"]
                 pseudocounts=0,
                 threshold=None,
                 k=1,

                 # Input / Output
                 output_dir="glam2_out",

                 # Alignment Options
                 alignment_runs=10,
                 iterations=10000,
                 both_strands=False,
                 min_sequences=2,
                 min_aligned_columns=2,
                 max_aligned_columns=50,
                 initial_aligned_columns=20,

                 # Scoring Scheme
                 del_pseudocount=0.1,
                 no_del_pseudocount=2.0,
                 ins_pseudocount=0.02,
                 no_ins_pseudocount=1.0,
                 weight=1e+99,

                 # Search Algorithm
                 temperature=1.2,
                 cooling_factor=1.44,
                 temp_lower_bound=0.1,
                 seed=1,

                 muscle_obj=None,
                 weblogo_obj=None
                 ):
        """Return a Glam2 object with specified attribute values."""
        self.alphabet = alphabet
        self.gap_in_alphabet = gap_in_alphabet
        self.scoring_criteria = scoring_criteria
        self.pseudocounts = pseudocounts
        if threshold is None:
            if scoring_criteria == 'pwm':
                self.threshold = 1.0e-9
            else:
                self.threshold = 0.8
        else:
            self.threshold = threshold
        self.k = k
        self.output_dir = output_dir    # -O
        self.number_alignment_runs = alignment_runs    # -r
        self.number_iterations = iterations     # -n
        self.both_strands = both_strands    # -2
        self.min_sequences = min_sequences    # -z
        self.min_aligned_columns = min_aligned_columns    # -a
        self.max_aligned_columns = max_aligned_columns    # -b
        self.initial_aligned_columns = initial_aligned_columns    # -w
        self.del_pseudocount = del_pseudocount    # -D, deletion pseudocount
        self.no_del_pseudocount = no_del_pseudocount    # -E, no-deletion pseudocount
        self.ins_pseudocount = ins_pseudocount    # -I, insertion pseudocount
        self.no_ins_pseudocount = no_ins_pseudocount    # -J, no-insertion pseudocount
        self.weight = weight    # -q, weight for generic versus sequence-set-specific residue abundances
        self.temperature = temperature    # -t, initial temperature
        self.cooling_factor = cooling_factor    # -c, cooling factor per n iterations
        self.temp_lower_bound = temp_lower_bound    # -u, temperature lower bound
        self.seed = seed    # -s, seed for pseudo-random numbers

        self.muscle_obj = muscle_obj
        self.weblogo_obj = weblogo_obj

        # Parsing data as a class
        self.record = None
        # no. of seqs in input file, to be set by fit()
        self.n_seqs = 0
        # number of motives found = number of alignment runs
        self.nmotifs = 0
        # list-of-strings representation of motifs
        self.motives_list = list()
        # aligned list-of-strings of motifs
        self.aligned_motives_list = list()
        # list of sequence logos created with WebLogo
        self.logos = list()

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

        if self.del_pseudocount != 0.1:
            params += " -D " + str(self.del_pseudocount)

        if self.no_del_pseudocount != 2.0:
            params += " -E " + str(self.no_del_pseudocount)

        if self.ins_pseudocount != 0.02:
            params += " -I " + str(self.ins_pseudocount)

        if self.no_ins_pseudocount != 1.0:
            params += " -J " + str(self.no_ins_pseudocount)

        if self.weight != 1e+99:
            params += " -q " + str(self.weight)

        if self.temperature != 1.2:
            params += " -t " + str(self.temperature)

        if self.cooling_factor != 1.44:
            params += " -c " + str(self.cooling_factor)

        if self.temp_lower_bound != 0.1:
            params += " -u " + str(self.temp_lower_bound)

        if self.seed != 1:
            params += " -s " + str(self.seed)

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

        # logger.info(stdout)

    def _parse_glam2_outfile(self, output_dir):
        record = Record()
        filename = os.path.join(output_dir, 'glam2.meme')
        with open(filename) as f1:
            record._read_version(f1)
            record._read_alphabet(f1)

        filename = os.path.join(output_dir, 'glam2.txt')
        record._parse_output(filename)
        return record

    def fit(self, fasta_file=''):
        """Save the output of Glam2 and parse it."""
        if not fasta_file:
            return NameError('Input fasta file not specified')

        cmd_params = self._make_param_string()
        self._command_exec(fasta_file, cmd_params)

        record = self._parse_glam2_outfile(self.output_dir)
        self.record = record
        self.nmotifs = record.nmotifs

        self.motives_list = record.motives_list[:]

        # create PWMs
        motives_list = self.motives_list[:]
        super(Glam2, self).fit(motives=motives_list)

    def fit_predict(self, fasta_file='', return_list=False):
        """Run fit and predict."""
        self.fit(fasta_file=fasta_file)
        return self.predict(input_seqs=fasta_file,
                            return_list=return_list)

    def fit_transform(self, fasta_file='', return_match=False):
        """Run fit and transform."""
        self.fit(fasta_file=fasta_file)
        return self.transform(input_seqs=fasta_file,
                              return_match=return_match)

    def display_glam2_help(self):
        """Display GLAM2 command line help."""
        command = "glam2 -h"
        print commands.getoutput(command)


class Record(object):
    """A class to store information after parsing output files."""

    def __init__(self):
        """init."""
        self.version = ""
        self.alphabet = None
        self.n_seqs = None
        self.nmotifs = None
        self.motives_list = []
        self.skeletons = []

    def _read_version(self, handle):
        for line in handle:
            if line.startswith('MEME'):
                break
        else:
            raise ValueError("Improper input files. One of the files should contain a line starting 'MEME'.")
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
        ls = line.split(' ')
        alph_letters = ls[1]

        if alph_letters == 'ACGT':
            self.alphabet = IUPAC.unambiguous_dna
        elif alph_letters == 'ACGU':
            self.alphabet = IUPAC.unambiguous_dna
        else:
            self.alphabet = IUPAC.protein

    def _read_nseqs(self, lines):
        for line in lines:
            if line.startswith('Sequences'):
                break
        else:
            raise ValueError("Improper input files. One of the files should contain a line starting 'Sequences'.")

        line = line.strip()
        ls = line.split(' ')

        n_seqs = ls[1]
        return n_seqs

    def _read_scores(self, lines):
        line_nums = []
        scores = []
        columns = []
        sequences = []

        for line in lines:
            if line.startswith('Score:'):
                line_nums.append(lines.index(line))
                l = line.strip()
                ls = l.split()
                scores.append(float(ls[1]))
                columns.append(int(ls[3]))
                sequences.append(int(ls[5]))

        return line_nums, scores, columns, sequences

    def _read_motives(self, n_motifs, lines, line_nums, seqs):
        motives = []
        motif_skeleton = []
        n_motifs = len(line_nums)
        for i in range(n_motifs):
            for j in range(line_nums[i], line_nums[i] + seqs[i] + 10):
                if '*' in lines[j]:
                    line = lines[j].replace('.', '-')
                    line = line.strip()
                    ls = line.split()
                    motif_skeleton.append(ls)

                    motif = []
                    for k in range(j + 1, j + seqs[i] + 1):    # TODO: missing last motif
                        line = lines[k].replace('.', '-')
                        line = line.strip()
                        ls = line.split()
                        mot = ls[2].upper()
                        motif.append((ls[0], mot))
                    motives.append(motif)
                    break
        return motives, motif_skeleton

    def _parse_output(self, filename):
        with open(filename) as f2:
            lines = f2.readlines()

        self.n_seqs = self._read_nseqs(lines)
        line_nums, scores, columns, seqs = self._read_scores(lines)

        n_motifs = len(scores)
        self.nmotifs = n_motifs

        motives_list, motif_skeles = self._read_motives(n_motifs, lines, line_nums, seqs)
        self.motives_list = motives_list[:]
        self.skeletons = motif_skeles[:]
