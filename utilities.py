"""Utility classes and functions used by MotifWrapper class."""
from StringIO import StringIO

from Bio import SeqIO, motifs
from Bio.Align.Applications import MuscleCommandline
from Bio.Alphabet import Gapped, IUPAC

from Bio.Seq import Seq


from Bio.SeqRecord import SeqRecord

from IPython.display import Image, display

from corebio.seq import Alphabet, SeqList

import numpy as np

import weblogolib as wbl

import random

from Bio import MarkovModel

import math


def _fasta_to_fasta(lines):
    seq = ""
    for line in lines:
        if line:
            if line[0] == '>':
                if seq:
                    yield seq
                    seq = ""
                line_str = str(line)
                yield line_str.strip()
            else:
                line_str = line.split()
                if line_str:
                    seq += str(line_str[0]).strip()
    if seq:
        yield seq


class MuscleAlignWrapper(object):
    """A wrapper to perform Muscle Alignment on sequences."""

    def __init__(self,
                 diags=False,
                 maxiters=16,
                 maxhours=None,
                 # TODO: check if this alphabet is required
                 # it over-rides tool.alphabet
                 alphabet='dna',  # ['dna', 'rna', 'protein']
                 ):
        """Initialize an instance."""
        self.diags = diags
        self.maxiters = maxiters
        self.maxhours = maxhours

        if alphabet == 'protein':
            self.alphabet = IUPAC.protein
        elif alphabet == 'rna':
            self.alphabet = IUPAC.unambiguous_rna
        else:
            self.alphabet = IUPAC.unambiguous_dna

    def _seq_to_stdin_fasta(self, seqs):
        # seperating headers
        headers, instances = [list(x) for x in zip(*seqs)]

        instances_seqrecord = []
        for i, j in enumerate(instances):
            instances_seqrecord.append(
                SeqRecord(Seq(j, self.alphabet), id=str(i)))

        handle = StringIO()
        SeqIO.write(instances_seqrecord, handle, "fasta")
        data = handle.getvalue()
        return headers, data

    def _perform_ma(self, data):
        params = {'maxiters': 7}
        if self.diags is True:
            params['diags'] = True
        if self.maxhours is not None:
            params['maxhours'] = self.maxhours

        muscle_cline = MuscleCommandline(**params)
        stdout, stderr = muscle_cline(stdin=data)
        return stdout

    def _fasta_to_seqs(self, headers, stdout):
        out = list(_fasta_to_fasta(stdout.split('\n')))
        motif_seqs = [''] * len(headers)
        for i in range(len(out[:-1]))[::2]:
            id = int(out[i].split(' ')[0].split('>')[1])
            motif_seqs[id] = out[i + 1]
        return zip(headers, motif_seqs)

    def transform(self, seqs=[]):
        """Carry out alignment."""
        headers, data = self._seq_to_stdin_fasta(seqs)
        stdout = self._perform_ma(data)
        aligned_seqs = self._fasta_to_seqs(headers, stdout)
        return aligned_seqs


class Weblogo(object):
    """A wrapper of weblogolib for creating sequence."""

    def __init__(self,

                 output_format='png',  # ['eps','png','png_print','jpeg']
                 stacks_per_line=40,
                 sequence_type='dna',  # ['protein','dna','rna']
                 ignore_lower_case=False,
                 units='bits',
                 # ['bits','nats','digits','kT','kJ/mol','kcal/mol','probability']
                 first_position=1,
                 logo_range=list(),
                 # composition = 'auto',
                 scale_stack_widths=True,
                 error_bars=True,
                 title='',
                 figure_label='',
                 show_x_axis=True,
                 x_label='',
                 show_y_axis=True,
                 y_label='',
                 y_axis_tic_spacing=1.0,
                 show_ends=False,
                 color_scheme='classic',
                 # ['auto','base','pairing','charge','chemistry','classic','monochrome']
                 resolution=200,
                 fineprint='',
                 ):
        """Initialize an instance."""
        options = wbl.LogoOptions()

        options.stacks_per_line = stacks_per_line
        options.sequence_type = sequence_type
        options.ignore_lower_case = ignore_lower_case
        options.unit_name = units
        options.first_index = first_position
        if logo_range:
            options.logo_start = logo_range[0]
            options.logo_end = logo_range[1]
        options.scale_width = scale_stack_widths
        options.show_errorbars = error_bars
        if title:
            options.title = title
        if figure_label:
            options.logo_label = figure_label
        options.show_xaxis = show_x_axis
        if x_label:
            options.xaxis_label = x_label
        options.show_yaxis = show_y_axis
        if y_label:
            options.yaxis_label = y_label
        options.yaxis_tic_interval = y_axis_tic_spacing
        options.show_ends = show_ends
        options.color_scheme = wbl.std_color_schemes[color_scheme]
        options.resolution = resolution
        if fineprint:
            options.fineprint = fineprint

        self.options = options
        self.output_format = output_format

    def create_logo(self, seqs=[]):
        """Create sequence logo for input sequences."""
        # seperate headers
        headers, instances = [list(x) for x in zip(*seqs)]

        if self.options.sequence_type is 'rna':
            alphabet = Alphabet('ACGU')
        elif self.options.sequence_type is 'protein':
            alphabet = Alphabet('ACDEFGHIKLMNPQRSTVWY')
        else:
            alphabet = Alphabet('ACGT')
        motif_corebio = SeqList(alist=instances, alphabet=alphabet)
        data = wbl.LogoData().from_seqs(motif_corebio)

        format = wbl.LogoFormat(data, self.options)

        if self.output_format == 'png':
            return wbl.png_formatter(data, format)
        elif self.output_format == 'png_print':
            return wbl.png_print_formatter(data, format)
        elif self.output_format == 'jpeg':
            return wbl.jpeg_formatter(data, format)
        else:
            return wbl.eps_formatter(data, format)


class MotifWrapper(object):
    """A generic wrapper for motif discovery tools."""

    def __init__(self,
                 alphabet='dna',  # ['dna', 'rna', 'protein']
                 gap_in_alphabet=True,
                 scoring_criteria='pwm',    # ['pwm','hmm']
                 pseudocounts=0,    # integer or dictionary {'A':0, 'C': 0, 'G': 0, 'T': 0}
                 threshold=None,    # scoring threshold
                 k=1,    # top-k scores returned for hmm score

                 # parameters for Muscle Alignment
                 ma_diags=False,
                 ma_maxiters=16,
                 ma_maxhours=None,

                 muscle_obj=None,
                 weblogo_obj=None
                 ):
        """Initialize an instance of MotifWrapper."""
        self.pseudocounts = pseudocounts
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

        self.muscle_obj = muscle_obj
        self.weblogo_obj = weblogo_obj

        self.ma_diags = ma_diags
        self.ma_maxiters = ma_maxiters
        self.ma_maxhours = ma_maxhours

        # list of PWMs corresponding to each motif
        self.pwms_list = list()
        # list-of-strings representation of motifs
        self.motives_list = list()
        # aligned list-of-strings of motifs, created by display_logo method
        self.aligned_motives_list = list()
        # number of motives, over-ridden by inherited class
        self.nmotifs = 0
        # list of sequence logos created with WebLogo
        self.logos = list()
        # list of Hidden Markov Model for each motif
        self.hmms_list = list()

    def _get_pwm(self, input_motif=list()):
        # seperate headers from sequences
        headers, instances = [list(x) for x in zip(*input_motif)]
        motif_seq = list()

        if self.alphabet == 'protein':
            alphabet = IUPAC.protein
        elif self.alphabet == 'rna':
            alphabet = IUPAC.unambiguous_rna
        else:
            alphabet = IUPAC.unambiguous_dna

        if self.gap_in_alphabet is True:
            alphabet = Gapped(alphabet, "-")

        for i in instances:
            # motif as Bio.Seq instance
            motif_seq.append(Seq(i, alphabet))

        motif_obj = motifs.create(motif_seq)
        return motif_obj.counts.normalize(self.pseudocounts)

    def fit(self, motives=list()):
        """Compute PWMs or HMMs for each found motif."""
        if self.scoring_criteria == 'pwm':
            pwms = list()
            for i in range(len(motives)):
                pwms.append(self._get_pwm(input_motif=motives[i]))
            self.pwms_list = pwms[:]
        else:
            if self.alphabet == 'protein':
                alphabet = 'ACDEFGHIKLMNPQRSTVWY'
            elif self.alphabet == 'rna':
                alphabet = 'ACGU'
            else:
                alphabet = 'ACGT'

            if self.gap_in_alphabet is True:
                alphabet += '-'

            hmms = list()
            for i in range(len(motives)):
                hmms.append(self._create_mm(motif_num=i + 1, alphabet=alphabet))
            self.hmms_list = hmms[:]

    def display(self, motif_num=None):
        """Display PWM of motives as dictionaries."""
        if motif_num is None:
            for i in range(len(self.pwms_list)):
                print self.pwms_list[i]
        else:
            print self.pwms_list[motif_num - 1]

    def _create_matrix(self, motif_num):
        m = []
        pwm_i = self.pwms_list[motif_num]
        symbols = pwm_i.alphabet.letters
        for s in sorted(symbols):
            m.append(pwm_i[s])
        return np.array(m)

    def matrix(self, motif_num=None):
        """If motif_num not specified, returns a list of numpy arrays."""
        if motif_num is None:
            matrix_list = []
            for i in range(len(self.pwms_list)):
                matrix_list.append(self._create_matrix(i))
            return matrix_list
        else:
            return self._create_matrix(motif_num - 1)

    def _eval_pwm(self, motif_num=1, seq=''):
        pwm_i = self.pwms_list[motif_num - 1]
        seq_len = len(seq)
        motif_len = len(pwm_i.itervalues().next())
        if seq_len < motif_len:
            raise ValueError('Sequence must be at least as long as the motif')
        scores = list()
        for i in range(seq_len - motif_len + 1):
            segment_score = 1
            for j in range(motif_len):
                letter = seq[i + j]
                segment_score *= pwm_i[letter][j]
            scores.append(segment_score)
        # zero padding
        for i in range(seq_len - len(scores)):
            scores.append(0)
        return scores

    def _parse_fasta_file(self, fasta_file):
        headers = list()
        seqs = list()
        fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')
        for fasta in fasta_sequences:
            seq_name, sequence = fasta.id, str(fasta.seq)
            headers.append(seq_name)
            seqs.append(sequence)
        return zip(headers, seqs)

    def _predict_pwm(self, sequences, seq_lists):
        for i, s in enumerate(sequences):
            for j in range(len(self.pwms_list)):
                score = self._eval_pwm(motif_num=j + 1, seq=s)
                score_max = max(score)
                for scr in score:
                    if scr > self.threshold:
                        if score_max > self.threshold:
                            seq_lists[i].append(j)
        return seq_lists

    def _eval_mm(self, motif_num=1, seq=''):
        """Return log_score_list of a sequence according to motif's HMM."""
        mm = self.hmms_list[motif_num - 1]
        hidden_states = len(mm.states)
        seq_len = len(seq)

        if seq_len < hidden_states:
            raise ValueError('Sequence must be at least as long as the motif')
        score = list()
        for i in range(seq_len - hidden_states + 1):
            seq_segment = seq[i:i + hidden_states - 1]
            result = MarkovModel.find_states(mm, seq_segment)
            score.append(result[0][1])

        eps = 1e-100
        log_score = [math.log(x + eps) for x in score]
        # zero padding
        for i in range(len(seq) - len(score)):
            log_score.append(0)
        return log_score

    def _predict_mm(self, sequences, seq_lists):
        for i, s in enumerate(sequences):
            for j in range(len(self.hmms_list)):
                score = self._eval_mm(motif_num=j + 1, seq=s)
                for scr in score:
                    if scr < self.threshold:
                        if max(score) < self.threshold:
                            seq_lists[i].append(j)
        return seq_lists

    def predict(self, input_seqs='', return_list=True):
        """Score each motif found with input sequence according to PWM."""
        if '.fa' in input_seqs:
            input_seqs = self._parse_fasta_file(fasta_file=input_seqs)
        headers, sequences = [list(x) for x in zip(*input_seqs)]

        # motives list for every sequence
        seq_lists = [[] for i in range(len(headers))]

        if self.scoring_criteria == 'pwm':
            seq_lists = self._predict_pwm(sequences=sequences, seq_lists=seq_lists)
        else:
            seq_lists = self._predict_mm(sequences=sequences, seq_lists=seq_lists)

        if return_list is True:
            return seq_lists
        return [len(i) for i in seq_lists]

    def _get_occurence_indexandscore_pwm(self, seq, motif_num):
        pwm_i = self.pwms_list[motif_num]
        seq_len = len(seq)
        motif_len = len(pwm_i.itervalues().next())

        scores = list()
        start_indexes = list()

        for i in range(seq_len - motif_len + 1):
            segment_score = 1
            for j in range(motif_len):
                letter = seq[i + j]
                segment_score *= pwm_i[letter][j]
            if segment_score > self.threshold:
                scores.append(segment_score)
                start_indexes.append(i + 1)
        last_indexes = [i + motif_len for i in start_indexes]
        return zip(start_indexes, last_indexes, scores)

    def _get_key(self, item):
        return item[2]

    def _get_occurence_indexandscore_mm(self, seq, motif_num):
        mm_i = self.hmms_list[motif_num]
        seq_len = len(seq)
        motif_len = len(mm_i.states)

        scores = list()
        start_indexes = list()

        for i in range(seq_len - motif_len + 1):
            segment_score = 0
            for j in range(motif_len):
                letter = seq[i + j]
                segment_score += MarkovModel.find_states(mm_i, letter)[0][1]
            if segment_score > self.threshold:
                scores.append(segment_score)
                start_indexes.append(i + 1)

        last_indexes = [i + motif_len for i in start_indexes]
        data = zip(start_indexes, last_indexes, scores)
        sorted_data = sorted(data, key=self._get_key, reverse=True)

        top_result = sorted_data[:self.k]
        return top_result

    def _transform_pwm(self, sequences, match_list):
        for i, s in enumerate(sequences):
            for j in range(len(self.pwms_list)):
                occs = self._get_occurence_indexandscore_pwm(seq=s, motif_num=j)
                match_list[i][j] = occs
        return match_list

    def _transform_mm(self, sequences, match_list):
        for i, s in enumerate(sequences):
            for j in range(len(self.hmms_list)):
                occs = self._get_occurence_indexandscore_mm(seq=s, motif_num=j)
                match_list[i][j] = occs
        return match_list

    def transform(self, input_seqs='', return_match=True):
        """Return summary of matches of each motif.

        return_match: True returns (start, end, score)-tuple for each
        motif found, False returns a summary of motives found
        """
        if '.fa' in input_seqs:
            input_seqs = self._parse_fasta_file(fasta_file=input_seqs)
        headers, sequences = [list(x) for x in zip(*input_seqs)]

        try:
            n_motives = len(self.pwms_list)
        except AttributeError:
            n_motives = len(self.hmms_list)

        match_list = [[[] for j in range(n_motives)]for i in range(len(headers))]

        if self.scoring_criteria == 'pwm':
            match_list = self._transform_pwm(sequences=sequences, match_list=match_list)
        else:
            match_list = self._transform_mm(sequences=sequences, match_list=match_list)

        if return_match is False:
            match_nums = [[] for i in range(len(headers))]
            for i in range(len(headers)):
                for j in range(n_motives):
                    if match_list[i][j]:
                        match_nums[i].append(1)
                    else:
                        match_nums[i].append(0)
            return match_nums
        return match_list

    def align_motives(self):
        """Perform Muscle Alignment on motives."""
        motives = list(self.motives_list)
        aligned_motives = list()

        if self.muscle_obj is None:
            self.muscle_obj = MuscleAlignWrapper(alphabet=self.alphabet)
        ma = self.muscle_obj

        for i in range(self.nmotifs):
            aligned_motives.append(ma.transform(seqs=motives[i]))

        self.aligned_motives_list = aligned_motives[:]

    def _get_logos_list(self, do_alignment=True):
        logos_list = list()

        motives = list(self.motives_list)
        if do_alignment is True:
            if not self.aligned_motives_list:
                self.align_motives()
            motives = list(self.aligned_motives_list)

        if self.weblogo_obj is None:
            self.weblogo_obj = Weblogo()
        wl = self.weblogo_obj
        wl.sequence_type = self.alphabet

        for i in range(self.nmotifs):
            logo = wl.create_logo(seqs=motives[i])
            logos_list.append(logo)
        return logos_list

    def display_logo(self, motif_num=None, do_alignment=True):
        """Display sequence logos of motives.

        motif_num: None displays all motives, i displays i-th motif
        do_alignment: Perform Muscle Alignment on motives before creating logos
        """
        self.logos = self._get_logos_list(do_alignment=do_alignment)[:]

        if motif_num is not None:
            display(Image(self.logos[motif_num - 1]))
        else:
            for i in range(self.nmotifs):
                display(Image(self.logos[i]))

    def _seq_to_columns(self, motif):
        motif_len = len(motif[0])
        cols = list()
        for i in range(motif_len):
            col = list()
            for j in range(len(motif)):
                col.append(motif[j][i])
            cols.append(col)
        return cols

    def _delete_gaps(self, columns):
        # List to store columns with no gaps
        new_columns = list()
        for col in columns:
            count = dict()
            for i in col:
                if i not in count.keys():
                    count[i] = 1
                else:
                    count[i] += 1

            freq_letter = max(count, key=count.get)
            non_gap_letters = [x for x in col if x != 'a']
            # If gap is one of the most frequent letters,
            # then discard the column
            # Discards even if half of the column has gap letter
            if freq_letter == '-':
                continue
            # Else, replace all gaps in column with one of the other letters
            else:
                for i, j in enumerate(col):
                    if j == '-':
                        col[i] = random.choice(non_gap_letters)
                new_columns.append(col)
        return new_columns

    def _columns_to_seqs(self, columns):
        n_seqs = len(columns[0])
        seqs = [[] for i in range(n_seqs)]
        for col in columns:
            for i in range(n_seqs):
                seqs[i].append(col[i])
        # concatenation of single letters into strings
        for i, s in enumerate(seqs):
            seqs[i] = ''.join(s)
        return seqs

    def _get_new_seqs(self, motif):
        columns = self._seq_to_columns(motif)
        new_columns = self._delete_gaps(columns)
        new_seqs = self._columns_to_seqs(new_columns)
        return new_seqs

    def adapt_motives(self, motives):
        """Perform adaption for motives of different lengths.

        If a single motif consists of instances of different lengths,
        then adaption trims the motif by removing columns with more
        gaps than characters.
        """
        modified_motives_list = list()
        for m in motives:
            heads = list()
            seqs = list()
            for j, k in m:
                heads.append(j)
                seqs.append(k)
                new_seqs = self._get_new_seqs(seqs)
            modified_motives_list.append(zip(heads, new_seqs))
        return modified_motives_list

    def _create_mm(self, motif_num, alphabet):
        try:
            # Only EDeN has original_motives_list
            input_motif = self.original_motives_list[motif_num - 1]
        except AttributeError:
            input_motif = self.motives_list[motif_num - 1]

        headers, instances = [list(x) for x in zip(*input_motif)]

        lengths = [len(instances[i]) for i in range(len(instances))]
        median_len = int(math.ceil(np.median(lengths)))

        # Hidden states for Markov Model
        states = [str(i + 1) for i in range(median_len)]

        print "original samples: %d" % len(instances)
        print "states:", len(states)
        # under sampling
        if (len(instances) * len(states)) > 500:
            samples = 500 / len(states)
            # samples = 50    # fixed sampling
            print 'sample size = %d' % samples
            instances = random.sample(instances, samples)

        instances = random.sample(instances, samples)

        mm = MarkovModel.train_bw(states=states,
                                  alphabet=alphabet,
                                  training_data=instances)
        return mm

    def _score_mm(self, motif_num=1, seq=''):
        return self._eval_mm(motif_num=motif_num, seq=seq)

    def _score_pwm(self, motif_num=1, seq=''):
        pwm_i = self.pwms_list[motif_num - 1]
        seq_len = len(seq)
        motif_len = len(pwm_i.itervalues().next())
        if seq_len < motif_len:
            raise ValueError('Sequence must be at least as long as the motif')
        score_mat = np.zeros((seq_len - motif_len + 1, seq_len))
        for i in range(seq_len - motif_len + 1):
            segment_score = 1
            for j in range(motif_len):
                letter = seq[i + j]
                segment_score *= pwm_i[letter][j]
            for j in range(i, i + motif_len):
                score_mat[i][j] = segment_score
        # voting
        max_score = [max(score_mat[:, i]) for i in range(seq_len)]
        return max_score

    def score(self, motif_num=1, seq=''):
        """Return score of the test seq according to the scoring criteria."""
        if self.scoring_criteria == "pwm":
            return self._score_pwm(motif_num=motif_num, seq=seq)
        return self._score_mm(motif_num=motif_num, seq=seq)
