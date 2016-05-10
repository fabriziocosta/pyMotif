from StringIO import StringIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline

from weblogolib import *
from corebio.seq import SeqList, Alphabet
from IPython.display import Image

from Bio import motifs


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


class MuscleAlign_wrapper(object):

    def __init__(self,
                 diags=False,
                 maxiters=16,
                 maxhours=None,
                 ):
        self.diags = diags
        self.maxiters = maxiters
        self.maxhours = maxhours

    def transform(self, seqs=[]):

        # seperating headers
        headers, instances = [list(x) for x in zip(*seqs)]
        num_instances = len(instances)

        instances_SeqRecord = []
        for i, j in enumerate(instances):
            instances_SeqRecord.append(SeqRecord(Seq(j, IUPAC.unambiguous_dna), id=str(i)))

        handle = StringIO()
        SeqIO.write(instances_SeqRecord, handle, "fasta")
        data = handle.getvalue()

        params = {'maxiters': 7}
        if self.diags is True:
            params['diags'] = True
        if self.maxhours is not None:
            params['maxhours'] = self.maxhours

        muscle_cline = MuscleCommandline(**params)
        stdout, stderr = muscle_cline(stdin=data)

        out = list(_fasta_to_fasta(stdout.split('\n')))
        motif_seqs = [''] * num_instances
        for i in range(len(out[:-1]))[::2]:
            ID = int(out[i].split(' ')[0].split('>')[1])
            motif_seqs[ID] = out[i + 1]

        return zip(headers, motif_seqs)


class Weblogo(object):

    def __init__(self,

                 output_format='png',  # [eps, png, png_print, jpeg]
                 stacks_per_line=40,
                 sequence_type='dna',  # [protein, dna, rna]
                 ignore_lower_case=False,
                 units='bits',  # ['bits', 'nats', 'digits', 'kT', 'kJ/mol', 'kcal/mol', 'probability']
                 first_position=1,
                 logo_range=list(),
                 #composition = 'auto',
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
                 color_scheme='classic',  # [auto, base, pairing, charge, chemistry, classic, monochrome]
                 resolution=96,
                 fineprint='',
                 ):

        options = LogoOptions()

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
        options.color_scheme = std_color_schemes[color_scheme]
        options.resolution = resolution
        if fineprint:
            options.fineprint = fineprint

        self.options = options
        self.output_format = output_format

    def create_logo(self, seqs=[]):
        headers, instances = [list(x) for x in zip(*seqs)]  # seperating headers

        if self.options.sequence_type is 'rna':
            alphabet = Alphabet('ACGU')
        elif self.options.sequence_type is 'protein':
            alphabet = Alphabet('ACDEFGHIKLMNPQRSTVWY')
        else:
            alphabet = Alphabet('AGCT')
        motif_corebio = SeqList(alist=instances, alphabet=alphabet)
        data = LogoData().from_seqs(motif_corebio)

        format = LogoFormat(data, self.options)

        if self.output_format == 'png':
            return png_formatter(data, format)
        elif self.output_format == 'png_print':
            return png_print_formatter(data, format)
        elif self.output_format == 'jpeg':
            return jpeg_formatter(data, format)
        else:
            return eps_formatter(data, format)


class PWM(object):

    def __init__(self,
                 alphabet='dna',  # ['dna', 'rna', 'protein']
                 pseudocounts=0,    # {'A':0, 'C': 0, 'G': 0, 'T': 0}
                 ):

        self.pseudocounts = pseudocounts

        if alphabet == 'protein':
            self.alphabet = IUPAC.protein
        elif alphabet == 'rna':
            self.alphabet = IUPAC.unambiguous_rna
        else:
            self.alphabet = IUPAC.unambiguous_dna

        self.pwms_list = list()

    def _get_pwm(self, input_motif=list()):

        headers, instances = [list(x) for x in zip(*input_motif)]  # seperating headers

        motif_seq = list()
        for i in instances:
            motif_seq.append(Seq(i, self.alphabet))  # motif as Bio.Seq instance

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

    def score(self, motif_num=1, seq=''):
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
        return scores

    def predict(self, motif_num=1, threshold=0.0001, seq=''):
        pwm_i = self.pwms_list[motif_num - 1]
        seq_len = len(seq)
        motif_len = len(pwm_i.itervalues().next())

        scores = list()
        indexes = list()

        for i in range(seq_len - motif_len + 1):
            segment_score = 1
            for j in range(motif_len):
                letter = seq[i + j]
                segment_score *= pwm_i[letter][j]
            if segment_score > threshold:
                scores.append(segment_score)
                indexes.append(i)

        return zip(indexes, scores)

    def _parse_input_file(self, fasta_file):
        headers = list()
        seqs = list()

        fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')
        for fasta in fasta_sequences:
            seq_name, sequence = fasta.id, str(fasta.seq)

            headers.append(seq_name)
            seqs.append(sequence)
        return zip(headers, seqs)

    def predict_new(self, fasta_file='', return_list=True, threshold=1.0e-9, append_score=False):
        input_seqs = self._parse_input_file(fasta_file)
        headers, sequences = [list(x) for x in zip(*input_seqs)]

        seq_lists = [[] for i in range(len(headers))]  # motifs list for every sequence

        for i, s in enumerate(sequences):
            for j in range(len(self.pwms_list)):  # TODO: adjust for self
                score = self.score(motif_num=j + 1, seq=s)  # TODO: adjust for self
                for scr in score:
                    if scr > threshold:
                        # seq_lists[i].append(j)
                        if sum(score) > threshold:
                            if append_score == True:  # TODO: remove append score, only else's line remains
                                seq_lists[i].append(score)
                            else:
                                seq_lists[i].append(j)

        if return_list == True:
            return seq_lists
        return [len(i) for i in seq_lists]

    def _get_occurence_indexandscore(self, seq, motif_num, threshold):
        pwm_i = self.pwms_list[motif_num]  # TODO: adjust for self, pwm_list[motif_num]
        seq_len = len(seq)
        motif_len = len(pwm_i.itervalues().next())

        scores = list()
        indexes = list()

        for i in range(seq_len - motif_len + 1):
            segment_score = 1
            for j in range(motif_len):
                letter = seq[i + j]
                segment_score *= pwm_i[letter][j]
            if segment_score > threshold:
                scores.append(segment_score)
                indexes.append(i + 1)
        indexes = [(i, i + motif_len) for i in indexes]
        return indexes, scores

    def transform_new(self, fasta_file='', return_match=True, threshold=1.0e-9):
        input_seqs = self._parse_input_file(fasta_file)
        headers, sequences = [list(x) for x in zip(*input_seqs)]

        match_list = [[[] for j in range(len(self.pwms_list))] for i in range(len(headers))]  # TODO: adjust for self

        for i, s in enumerate(sequences):
            for j in range(len(self.pwms_list)):  # TODO: adjust for self, nmotifs
                occs, scores = self._get_occurence_indexandscore(seq=s, motif_num=j, threshold=threshold)
                match_list[i][j] = occs

        if return_match is False:
            match_nums = [[] for i in range(len(headers))]
            for i in range(len(headers)):
                for j in range(len(self.pwms_list)):  # TODO: adjust for self, nmotifs
                    if match_list[i][j]:
                        match_nums[i].append(1)
                    else:
                        match_nums[i].append(0)
            return match_nums
        return match_list
    # TODO: replace predict() with predict_new()
    # TODO: Make one Abstract class
    # TODO: attributes of WL and MA classes
