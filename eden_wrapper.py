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
                 # radius=None,
                 # distance=None,
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
        self.motives_list = list()
        # aligned list-of-strings of motifs, created by display_logo method
        self.aligned_motives_list = list()
        # list of sequence logos created with WebLogo
        self.logos = list()
        # threshold for scoring sequences
        self.threshold = threshold
        # TODO: remove
        self.motives_db = None

    def _get_motives_list(self, db):
        motives = list()
        for i in db.keys():
            motives.append(db[i])
        self.motives_list = motives[:]

        aligned_motives = []
        ma = MuscleAlignWrapper()
        for i in range(len(motives)):
            aligned_motives.append(ma.transform(seqs=motives[i]))
        self.aligned_motives = aligned_motives[:]

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
                           # pre_processor_block_size=self.pre_processor_block_size,
                           random_state=self.random_state,
                           )
        sm.fit(seqs=seqs, neg_seqs=neg_seqs)

        # TODO: remove db
        self.motives_db = sm.motives_db
        self.nmotifs = len(sm.motives_db.keys())
        self._get_motives_list(sm.motives_db)

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
