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

    def __init__(self,
                 diags=False,
                 maxiters=16,
                 maxhours=None,
                 ):
        self.diags = diags
        self.maxiters = maxiters
        self.maxhours = maxhours

    def _seq_to_stdin_fasta(self, seqs):
        # seperating headers
        headers, instances = [list(x) for x in zip(*seqs)]
        # num_instances = len(instances)

        instances_seqrecord = []
        for i, j in enumerate(instances):
            instances_seqrecord.append(
                SeqRecord(Seq(j, IUPAC.unambiguous_dna), id=str(i)))
            # TODO: alphabet according to data

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
        headers, data = self._seq_to_stdin_fasta(seqs)
        stdout = self._perform_ma(data)
        aligned_seqs = self._fasta_to_seqs(headers, stdout)
        return aligned_seqs


class Weblogo(object):

    def __init__(self,

                 output_format='png',  # ['eps','png','png_print','jpeg']
                 stacks_per_line=40,
                 sequence_type='dna',  # ['protein','dna','rna']
                 ignore_lower_case=False,
                 # ['bits','nats','digits','kT','kJ/mol','kcal/mol','probability']
                 units='bits',
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
                 # ['auto','base','pairing','charge','chemistry','classic','monochrome']
                 color_scheme='classic',
                 resolution=96,
                 fineprint='',
                 ):

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
        # seperating headers
        headers, instances = [list(x)
                              for x in zip(*seqs)]

        if self.options.sequence_type is 'rna':
            alphabet = Alphabet('ACGU')
        elif self.options.sequence_type is 'protein':
            alphabet = Alphabet('ACDEFGHIKLMNPQRSTVWY')
        else:
            alphabet = Alphabet('AGCT')
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

    def __init__(self,
                 alphabet='dna',  # ['dna', 'rna', 'protein']
                 pseudocounts=0,    # {'A':0, 'C': 0, 'G': 0, 'T': 0}

                 # parameters for Muscle Alignment
                 ma_diags=False,
                 ma_maxiters=16,
                 ma_maxhours=None,

                 # parameters for WebLogo
                 wl_output_format='png',  # ['eps', 'png', 'png_print', 'jpeg']
                 wl_stacks_per_line=40,
                 wl_ignore_lower_case=False,
                 wl_units='bits',
                 # ['bits','nats','digits','kT','kJ/mol','kcal/mol','probability']
                 wl_first_position=1,
                 wl_logo_range=list(),
                 wl_scale_stack_widths=True,
                 wl_error_bars=True,
                 wl_title='',
                 wl_figure_label='',
                 wl_show_x_axis=True,
                 wl_x_label='',
                 wl_show_y_axis=True,
                 wl_y_label='',
                 wl_y_axis_tic_spacing=1.0,
                 wl_show_ends=False,
                 wl_color_scheme='classic',
                 # ['auto','base','pairing','charge','chemistry','classic','monochrome']
                 wl_resolution=96,
                 wl_fineprint='',
                 ):

        self.pseudocounts = pseudocounts
        self.alphabet = alphabet

        self.ma_diags = ma_diags
        self.ma_maxiters = ma_maxiters
        self.ma_maxhours = ma_maxhours

        self.wl_output_format = wl_output_format
        self.wl_stacks_per_line = wl_stacks_per_line
        self.wl_ignore_lower_case = wl_ignore_lower_case
        self.wl_units = wl_units
        self.wl_first_position = wl_first_position
        self.wl_logo_range = wl_logo_range
        self.wl_scale_stack_widths = wl_scale_stack_widths
        self.wl_error_bars = wl_error_bars
        self.wl_title = wl_title
        self.wl_figure_label = wl_figure_label
        self.wl_show_x_axis = wl_show_x_axis
        self.wl_x_label = wl_x_label
        self.wl_show_y_axis = wl_show_y_axis
        self.wl_y_label = wl_y_label
        self.wl_y_axis_tic_spacing = wl_y_axis_tic_spacing
        self.wl_show_ends = wl_show_ends
        self.wl_color_scheme = wl_color_scheme
        self.wl_resolution = wl_resolution
        self.wl_fineprint = wl_fineprint
        if self.alphabet == "dna":
            self.wl_sequence_type = "dna"
        elif self.alphabet == "rna":
            self.wl_sequence_type = "rna"
        elif self.alphabet == "protein":
            self.wl_sequence_type = "protein"
        else:
            self.wl_sequence_type = "auto"

        self.pwms_list = list()
        # list-of-strings representation of motifs
        self.motives_list = list()
        # aligned list-of-strings of motifs, created by display_logo method
        self.aligned_motives_list = list()
        # number of motives, over-ridden by inherited class
        self.nmotifs = 0
        # list of sequence logos created with WebLogo
        self.logos = list()

    def _get_pwm(self, input_motif=list()):
        # seperate headers from sequences
        headers, instances = [list(x) for x in zip(*input_motif)]
        motif_seq = list()

        if self.alphabet == 'protein':
            alphabet = Gapped(IUPAC.protein, "-")
        elif self.alphabet == 'rna':
            alphabet = Gapped(IUPAC.unambiguous_rna, "-")
        else:
            alphabet = Gapped(IUPAC.unambiguous_dna, "-")

        for i in instances:
            # motif as Bio.Seq instance
            motif_seq.append(Seq(i, alphabet))

        motif_obj = motifs.create(motif_seq)
        return motif_obj.counts.normalize(self.pseudocounts)

    def fit(self, motives=list()):
        pwms = list()
        for i in range(len(motives)):
            pwms.append(self._get_pwm(input_motif=motives[i]))
        self.pwms_list = pwms[:]

    def display(self, motif_num=None):
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
        """
        if motif_num not specified, returns a list of numpy arrays
        """
        if motif_num is None:
            matrix_list = []
            for i in range(len(self.pwms_list)):
                matrix_list.append(self._create_matrix(i))
            return matrix_list
        else:
            return self._create_matrix(motif_num - 1)

    def score(self, motif_num=1, seq='', zero_padding=False):
        """
        Scores a single sequence according to specified motif
        """
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
        if zero_padding is True:
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

    def predict(self,
                fasta_file='',
                return_list=True,
                threshold=1.0e-9,
                append_score=False):
        # TODO: remove append_score after debugging is complete
        # TODO: also remove append_score from Meme.predict and fit_predict
        input_seqs = self._parse_fasta_file(fasta_file)
        headers, sequences = [list(x) for x in zip(*input_seqs)]

        # motives list for every sequence
        seq_lists = [[] for i in range(len(headers))]

        for i, s in enumerate(sequences):
            for j in range(len(self.pwms_list)):
                score = self.score(motif_num=j + 1, seq=s)
                for scr in score:
                    if scr > threshold:
                        # seq_lists[i].append(j)
                        if sum(score) > threshold:
                            # TODO: remove append score, only else's line
                            # remains
                            if append_score is True:
                                seq_lists[i].append(score)
                            else:
                                seq_lists[i].append(j)
        if return_list is True:
            return seq_lists
        return [len(i) for i in seq_lists]

    def _get_occurence_indexandscore(self, seq, motif_num, threshold):
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
            if segment_score > threshold:
                scores.append(segment_score)
                start_indexes.append(i + 1)
        last_indexes = [i + motif_len for i in start_indexes]
        return zip(start_indexes, last_indexes, scores)

    def transform(self, fasta_file='', return_match=True, threshold=1.0e-9):
        input_seqs = self._parse_fasta_file(fasta_file)
        headers, sequences = [list(x) for x in zip(*input_seqs)]
        match_list = [
                     [
                         [] for j in range(len(self.pwms_list))
                     ]
            for i in range(len(headers))
        ]

        for i, s in enumerate(sequences):
            for j in range(len(self.pwms_list)):
                occs = self._get_occurence_indexandscore(
                    seq=s, motif_num=j, threshold=threshold)
                match_list[i][j] = occs

        if return_match is False:
            match_nums = [[] for i in range(len(headers))]
            for i in range(len(headers)):
                for j in range(len(self.pwms_list)):
                    if match_list[i][j]:
                        match_nums[i].append(1)
                    else:
                        match_nums[i].append(0)
            return match_nums
        return match_list

    def align_motives(self):
        motives = list(self.motives_list)
        aligned_motives = list()
        ma = MuscleAlignWrapper(diags=self.ma_diags,
                                maxiters=self.ma_maxiters,
                                maxhours=self.ma_maxhours,
                                )
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

        wl = Weblogo(output_format=self.wl_output_format,
                     stacks_per_line=self.wl_stacks_per_line,
                     sequence_type=self.wl_sequence_type,
                     ignore_lower_case=self.wl_ignore_lower_case,
                     units=self.wl_units,
                     first_position=self.wl_first_position,
                     logo_range=self.wl_logo_range,
                     # composition =
                     scale_stack_widths=self.wl_scale_stack_widths,
                     error_bars=self.wl_error_bars,
                     title=self.wl_title,
                     figure_label=self.wl_figure_label,
                     show_x_axis=self.wl_show_x_axis,
                     x_label=self.wl_x_label,
                     show_y_axis=self.wl_show_y_axis,
                     y_label=self.wl_y_label,
                     y_axis_tic_spacing=self.wl_y_axis_tic_spacing,
                     show_ends=self.wl_show_ends,
                     color_scheme=self.wl_color_scheme,
                     resolution=self.wl_resolution,
                     fineprint=self.wl_fineprint,
                     )
        for i in range(self.nmotifs):
            logo = wl.create_logo(seqs=motives[i])
            logos_list.append(logo)
        return logos_list

    def display_logo(self, motif_num=None, do_alignment=True):
        # Displays logos of all motifs if motif_num is not specified

        self.logos = self._get_logos_list(do_alignment=do_alignment)[:]

        if motif_num is not None:
            display(Image(self.logos[motif_num - 1]))
        else:
            for i in range(self.nmotifs):
                display(Image(self.logos[i]))
