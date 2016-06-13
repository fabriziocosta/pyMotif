"""A eden.motif.SequenceMotif wrapper for comparison against Meme."""

from utilities import MotifWrapper, MuscleAlignWrapper
from eden.motif import SequenceMotif


class EdenWrapper(MotifWrapper):
    """Wrapper for EDeN Sequence Motif."""

    def __init__(self,
                 alphabet='dna',
                 min_subarray_size=7,
                 max_subarray_size=10,
                 min_motif_count=1,
                 min_cluster_size=1,
                 training_size=None,
                 negative_ratio=1,
                 shuffle_order=2,
                 n_iter_search=1,
                 complexity=4,
                 # radius=None,    #TODO: check radius
                 # distance=None,    #TODO: check distance
                 nbits=20,
                 clustering_algorithm=None,
                 n_jobs=4,
                 n_blocks=8,
                 block_size=None,
                 pre_processor_n_jobs=4,
                 pre_processor_n_blocks=8,
                 pre_processor_block_size=None,
                 random_state=1,
                 pseudocounts=0,
                 threshold=1.0e-9
                 ):
        """Initialize a EdenWrapper object."""
        self.alphabet = alphabet
        self.min_subarray_size = min_subarray_size
        self.max_subarray_size = max_subarray_size
        self.min_motif_count = min_motif_count
        self.min_cluster_size = min_cluster_size
        self.training_size = training_size
        self.negative_ratio = negative_ratio
        self.shuffle_order = shuffle_order
        self.n_iter_search = n_iter_search
        self.complexity = complexity
        # self.radius = radius
        # self.distance = distance
        self.nbits = nbits
        self.clustering_algorithm = clustering_algorithm
        self.n_jobs = n_jobs
        self.n_blocks = n_blocks
        self.block_size = block_size
        self.pre_processor_n_jobs = pre_processor_n_jobs
        self.pre_processor_n_blocks = pre_processor_n_blocks
        self.pre_processor_block_size = pre_processor_block_size
        self.random_state = random_state

        # Number of motives found
        self.nmotifs = 0
        # over-rides same attribute of MotifWrapper class
        self.pseudocounts = pseudocounts
        # list-of-strings representation of motifs
        self.original_motives_list = list()
        # aligned list-of-strings of motifs;
        # also created by display_logo method
        self.aligned_motives_list = list()
        # modified motives with no gaps
        self.adapted_motives_list = list()
        # list of sequence logos created with WebLogo
        self.logos = list()
        # threshold for scoring sequences
        self.threshold = threshold

    def _get_motives_list(self, db):
        motives = list()
        for i in db.keys():
            motives.append(db[i])
        self.original_motives_list = motives[:]

    def _get_aligned_motives_list(self, motives):
        aligned_motives = []
        ma = MuscleAlignWrapper()
        for i in range(len(motives)):
            aligned_motives.append(ma.transform(seqs=motives[i]))
        self.aligned_motives_list = aligned_motives[:]

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
            # Get the most frequent letter of column
            freq_letter = max(count, key=count.get)
            # If gap is one of the most frequent letters,
            # then discard the column
            # Discards even if half of the column has gap letter
            if freq_letter == '-':
                continue
            # Else, replace all gaps in column with most frequent letter
            else:
                for i, j in enumerate(col):
                    if j == '-':
                        col[i] = freq_letter
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
            seqs[i] = [''.join(s)]
        return seqs

    def _get_new_seqs(self, motif):
        columns = self._seq_to_columns(motif)
        new_columns = self._delete_gaps(columns)
        new_seqs = self._columns_to_seqs(new_columns)
        return new_seqs

    def _adapt_motives(self, motives):
        modified_motives_list = list()
        for m in motives:
            heads = list()
            seqs = list()
            for j, k in m:
                heads.append(j)
                seqs.append(k)
                new_seqs = self._get_new_seqs(seqs)
            modified_motives_list.append(zip(heads, new_seqs))
        self.adapted_motives_list = modified_motives_list[:]

    def fit(self, seqs, neg_seqs=None):
        """Find motives with EDeN Sequence Motif."""
        sm = SequenceMotif(min_subarray_size=self.min_subarray_size,
                           max_subarray_size=self.max_subarray_size,
                           min_motif_count=self.min_motif_count,
                           min_cluster_size=self.min_cluster_size,
                           training_size=self.training_size,
                           negative_ratio=self.negative_ratio,
                           shuffle_order=self.shuffle_order,
                           n_iter_search=self.n_iter_search,
                           complexity=self.complexity,
                           # radius=self.radius,
                           # distance=self.distance,
                           nbits=self.nbits,
                           clustering_algorithm=self.clustering_algorithm,
                           n_jobs=self.n_jobs,
                           n_blocks=self.n_blocks,
                           block_size=self.block_size,
                           pre_processor_n_jobs=self.pre_processor_n_jobs,
                           pre_processor_n_blocks=self.pre_processor_n_blocks,
                           # TODO: check pre_processor_block_size
                           # pre_processor_block_size=self.pre_processor_block_size,
                           random_state=self.random_state,
                           )
        sm.fit(seqs=seqs, neg_seqs=neg_seqs)

        self.nmotifs = len(sm.motives_db.keys())
        self._get_motives_list(sm.motives_db)
        self._get_aligned_motives_list(self.original_motives_list)
        self._adapt_motives(self.aligned_motives_list)

        # create PWMs
        aligned_motives_list = self.aligned_motives_list[:]
        super(EdenWrapper, self).fit(motives=aligned_motives_list)

    def fit_predict(self,
                    seqs,
                    neg_seqs=None,
                    return_list=False,
                    ):
        """Find motives and return motif occurence list."""
        self.fit(seqs=seqs, neg_seqs=neg_seqs)
        return self.predict(input_seqs=seqs, return_list=return_list)

    def fit_transform(self,
                      seqs,
                      neg_seqs=None,
                      return_match=False,
                      ):
        """Find motives and return motif match list."""
        self.fit(seqs=seqs, neg_seqs=neg_seqs)
        return self.transform(input_seqs=seqs, return_match=return_match)
