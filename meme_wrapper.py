"""A wrapper of the motif discovery tool MEME."""
import os

from subprocess import PIPE, Popen

from Bio import motifs

from utilities import MotifWrapper

import logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)


class Meme(MotifWrapper):
    """
    Wrapper for MEME 4.11.0.

    Usage: meme <sequence file> [options]
    To see MEME help, use MEME.display_meme_help()
    """

    def __init__(self,
                 # Output
                 output_dir="meme_out",
                 text=False,

                 # Alphabet
                 alphabet="dna",  # ["dna", "rna", "protein"]
                 gap_in_alphabet=True,
                 scoring_criteria="pwm",    # ["pwm", "hmm"]

                 # Contributing Site Distribution
                 mod="zoops",

                 # Number of Motifs
                 nmotifs=1,    # 1 is the default value
                 evt=None,    # E-value
                 time=None,    # CPU time before timeout

                 # Number of Motif Occurences
                 nsites=None,    # ignored with OOPS
                 minsites=None,    # ignored with OOPS
                 maxsites=None,    # ignored with OOPS
                 wnsites=None,

                 # Motif Width
                 w=None,
                 minw=8,    # default value used by the meme tool
                 maxw=50,    # default value used by the meme tool
                 nomatrim=False,    # No Motif Alignment Trim
                 wg=None,
                 ws=None,
                 noendgaps=False,

                 # Markov Background Model
                 bfile="",    # Markov background model file

                 # Position-Specific Priors
                 psp="",    # Position Specific File

                 # Palindromes & Strands
                 revcomp=False,
                 pal=False,

                 # EM Algorithm
                 maxiter=None,
                 distance=None,
                 prior="",    # ["dirichlet","dmix","mega","megap","addone"]
                 b=None,
                 plib="",    # depends on alphabet of sequence

                 # Selecting Starts for EM
                 spfuzz=None,    # depends on spmap
                 spmap="",    # depends on alphabet, {uni, pam}
                 cons=None,

                 # Branching Search on EM Starts
                 heapsize=None,
                 x_branch=False,
                 w_branch=False,
                 bfactor=None,

                 # Miscellaneous
                 maxsize=None,
                 v=False,  # extensive message

                 # h = False # display usage message

                 # parameters for (parent) MotifWrapper
                 pseudocounts=0,    # alphabet pseudocount

                 muscle_obj=None,
                 weblogo_obj=None
                 ):
        """Initialize a MemeWrapper Object."""
        self.output_dir = output_dir
        self.text = text
        self.alphabet = alphabet
        self.gap_in_alphabet = gap_in_alphabet
        self.scoring_criteria = scoring_criteria
        if scoring_criteria == 'pwm':
            self.threshold = 1.0e-9
        else:
            self.threshold = -200    # TODO: examine threshold for hmm
        self.mod = mod
        self.nmotifs = nmotifs
        self.evt = evt
        self.time = time
        self.nsites = nsites
        self.minsites = minsites
        self.maxsites = maxsites
        self.wnsites = wnsites
        self.w = w
        self.minw = minw
        self.maxw = maxw
        self.nomatrim = nomatrim
        self.wg = wg
        self.ws = ws
        self.noendgaps = noendgaps
        self.bfile = bfile
        self.psp = psp
        self.revcomp = revcomp
        self.pal = pal
        self.maxiter = maxiter
        self.distance = distance
        self.prior = prior
        self.b = b
        self.plib = plib
        self.spmap = spmap
        self.spfuzz = spfuzz
        self.cons = cons
        self.heapsize = heapsize
        self.x_branch = x_branch
        self.w_branch = w_branch
        self.bfactor = bfactor
        self.maxsize = maxsize
        self.v = v

        self.muscle_obj = muscle_obj
        self.weblogo_obj = weblogo_obj

        # no. of seqs in input file, to be set by fit()
        self.n_seqs = 0
        # to store the names of sequences as given in input file to fit()
        self.seq_names = list()
        # list of motives, each represented by an object
        self.motives_db = list()
        # widths given in summary headers; length of each motif
        self.widths = list()
        # sites given in summary headers; num occurences
        self.sites = list()
        # e_values given in summary headers
        self.e_values = list()
        # over-rides same attribute of MotifWrapper class
        self.pseudocounts = pseudocounts
        # list-of-strings representation of motifs
        self.motives_list = list()
        # aligned list-of-strings of motifs
        self.aligned_motives_list = list()
        # list of sequence logos created with WebLogo
        self.logos = list()

    def _make_param_string(self):
        # Creates a string of parameters
        params = "-oc " + self.output_dir

        if self.nmotifs != 1:
            params += " -nmotifs " + str(self.nmotifs)

        if self.text is True:
            params += " -text"

        if self.alphabet == "dna":
            params += " -dna"
        elif self.alphabet == "rna":
            params += " -rna"
        else:  # self.alphabet == "protein":
            params += " -protein"
        # else:
        #    params += " -alph " + self.alphabet

        if self.mod == "zoops":
            params += " -mod zoops"
        elif self.mod == "anr":
            params += " -mod anr"
        else:
            params += " -mod oops"

        if self.evt is not None:
            params += " -evt " + str(self.evt)

        if self.time is not None:
            params += " -time " + str(self.time)

        if self.nsites is not None and self.mod != "oops":
            params += " -nsites " + str(self.nsites)

        if self.minsites is not None and self.mod != "oops":
            params += " -minsites " + str(self.minsites)

        if self.maxsites is not None and self.mod != "oops":
            params += " -maxsites " + str(self.maxsites)

        if self.wnsites is not None:
            params += " -wnsites " + str(self.wnsites)

        if self.w is not None:
            params += " -w " + str(self.w)
        elif self.minw != 8 or self.maxw != 50:
            params += " -minw " + str(self.minw) + " -maxw " + str(self.maxw)

        if self.nomatrim is True:
            params += " -nomatrim"

        if self.wg is not None:
            params += " -wg " + str(self.wg)

        if self.ws is not None:
            params += " -ws " + str(self.ws)

        if self.noendgaps is not False:
            params += " -noendgaps"

        if self.bfile:
            params += " -bfile " + self.bfile

        if self.psp:
            params += " -psp " + self.psp

        if self.revcomp is True:
            params += " -revcomp"

        if self.pal is True:
            params += " -pal"

        if self.maxiter is not None:
            params += " -maxiter" + str(self.maxiter)

        if self.distance is not None:
            params += " -distance " + str(self.distance)

        if self.prior:
            params += " -prior " + self.prior

        if self.b is not None:
            params += " -b " + str(self.b)

        if self.plib:
            params += " -plib " + self.plib

        if self.spfuzz is not None:
            params += " -spfuzz " + self.spfuzz

        if self.spmap:
            params += " -spmap " + str(self.spmap)

        if self.cons is not None:
            params += " -cons " + str(self.cons)

        if self.heapsize is not None:
            params += " -heapsize " + str(self.heapsize)

        if self.x_branch is True:
            params += " -x_branch"

        if self.w_branch is True:
            params += " -w_branch"

        if self.bfactor is not None:
            params += " -bfactor " + str(self.bfactor)

        if self.maxsize is not None:
            params += " -maxsize " + str(self.maxsize)

        if self.v is True:
            params += " -V"

        return params

    def _command_exec(self, fasta_file="", params=""):
        # execute the meme command and return output
        cmd = "meme " + fasta_file + " " + params
        io = Popen(cmd.split(" "), stdout=PIPE, stderr=PIPE)
        (stderr, stdout) = io.communicate()

        logger.info(stdout)

    def _get_stats(self, n_motifs):
        widths = list()
        sites = list()
        e_values = list()

        for i in range(n_motifs):
            widths.append(self.motives_db[i].length)
            sites.append(self.motives_db[i].num_occurrences)
            e_values.append(self.motives_db[i].evalue)

        self.widths = widths[:]
        self.sites = sites[:]
        self.e_values = e_values[:]

    def _get_motives_list(self):

        # list-of-strings representation of motives
        motives = []
        for i in range(self.nmotifs):
            motif_i = self.motives_db[i]
            motif_as_str = str(motif_i.instances).split('\n')[:-1]
            headers = []
            for i in range(len(motif_i.instances)):
                headers.append(motif_i.instances[i].sequence_name)
            motives.append(zip(headers, motif_as_str))
        return motives

    def fit(self, fasta_file=''):
        """Save the output of MEME and parse it."""
        if not fasta_file:
            raise NameError('Input fasta file not specified')

        cmd_params = self._make_param_string()
        self._command_exec(fasta_file, cmd_params)

        # parsing meme output file with Biopython
        filename = os.path.join(self.output_dir, 'meme.txt')
        handle = open(filename)
        record = motifs.parse(handle, 'meme')
        handle.close()

        # store names of sequences given as input to fit()
        self.seq_names = record.sequences[:]
        self.n_seqs = len(record.sequences)

        # create a list of motives, each represented by an object
        motives_db = list()
        for i in range(self.nmotifs):
            motives_db.append(record[i])
        self.motives_db = motives_db[:]

        # store length, number of occurences and e-value of each motif
        self._get_stats(self.nmotifs)

        # get string representation of motives
        self.motives_list = list(self._get_motives_list())

        # create PWMs
        motives_list = self.motives_list[:]
        super(Meme, self).fit(motives=motives_list)

    def _get_seq_header_list(self, fasta_file=''):

        header = list()
        with open(fasta_file, 'r') as f:
            for line in f.readlines():
                if line[0] == ">":
                    header.append(line.split(">")[1].split(" ")[0])
        return header

    def _get_motif_occurences(self):
        # returns a list of lists containing motif occurence sequences' names
        occurences = [[] for i in range(self.nmotifs)]

        for i in range(self.nmotifs):
            for j in range(len(self.motives_db[i].instances)):
                occurences[i].append(
                    self.motives_db[i].instances[j].sequence_name)
        return occurences

    def _get_match_list(self, header=None, motif_occurences=None):
        match_list = [
            [
                [] for j in range(self.nmotifs)
            ] for i in range(len(header))
        ]

        # for all seqs in transform() input
        for seq_id, seq_name in enumerate(header):
            # for all motives found
            for i in range(self.nmotifs):
                if seq_name in motif_occurences[i]:
                    # for multiple occurences of same motif in a sequence
                    repetitions = motif_occurences[i].count(seq_name)

                    if repetitions == 1:
                        seq_index = [motif_occurences[i].index(seq_name)]
                    else:
                        seq_index = [
                            k for k, val in enumerate(motif_occurences[i])
                            if val == seq_name]

                    for j in seq_index:
                        start_pos = self.motives_db[i].instances[j].start
                        motif_len = self.widths[i]
                        motif_location = zip(
                            [start_pos], [start_pos + motif_len])[0]
                        match_list[seq_id][i].append(motif_location)
        return match_list

    def fit_predict(self, fasta_file="", return_list=False):
        """Build a model and find motif predictions in training data."""
        self.fit(fasta_file=fasta_file)
        return self.predict(input_seqs=fasta_file,
                            return_list=return_list)

    def fit_transform(self, fasta_file="", return_match=False):
        """Build a model and find matches in input sequences."""
        self.fit(fasta_file=fasta_file)
        return self.transform(input_seqs=fasta_file, return_match=return_match)

    def display_meme_help(self):
        """Display MEME's CLI help."""
        cmd = "meme --help"
        io = Popen(cmd.split(" "), stdout=PIPE, stderr=PIPE)
        (error, output) = io.communicate()
        print output
