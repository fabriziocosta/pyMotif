from utilities import MotifWrapper
from eden.motif import SequenceMotif


class EdenWrapper(MotifWrapper):
    """Wrapper for EDeN Sequence Motif."""

    def __init__(self,
                 min_subarray_size=7,
                 max_subarray_size=10,
                 min_motif_count=1,
                 min_cluster_size=1,
                 training_size=None,
                 negative_ratio=1,
                 shuffle_order=2,
                 n_iter_search=1,
                 complexity=4,
                 radius=None,
                 distance=None,
                 nbits=20,
                 clustering_algorithm=None,
                 n_jobs=4,
                 n_blocks=8,
                 block_size=None,
                 pre_processor_n_jobs=4,
                 pre_processor_n_blocks=8,
                 pre_processor_block_size=None,
                 random_state=1,
                 ):
        """Initialize a EdenWrapper object."""
        self.min_subarray_size = min_subarray_size
        self.max_subarray_size = max_subarray_size
        self.min_motif_count = min_motif_count
        self.min_cluster_size = min_cluster_size
        self.training_size = training_size
        self.negative_ratio = negative_ratio
        self.shuffle_order = shuffle_order
        self.n_iter_search = n_iter_search
        self.complexity = complexity
        self.radius = radius
        self.distance = distance
        self.nbits = nbits
        self.clustering_algorithm = clustering_algorithm
        self.n_jobs = n_jobs
        self.n_blocks = n_blocks
        self.block_size = block_size
        self.pre_processor_n_jobs = self.pre_processor_n_jobs
        self.pre_processor_n_blocks = pre_processor_n_blocks
        self.pre_processor_block_size = pre_processor_block_size
        self.random_state = random_state

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
                           radius=self.radius,
                           distance=self.distance,
                           nbits=self.nbits,
                           clustering_algorithm=self.clustering_algorithm,
                           n_jobs=self.n_jobs,
                           n_blocks=self.n_blocks,
                           block_size=self.block_size,
                           pre_processor_n_jobs=self.pre_processor_n_jobs,
                           pre_processor_n_blocks=self.pre_processor_n_blocks,
                           pre_processor_block_size=self.pre_processor_block_size,
                           random_state=self.random_state,
                           )
        sm.fit(seqs=seqs, neg_seqs=neg_seqs)

        for cluster_id in sm.motives_db:
            cluster_seqs = [(i, motif) for i, (count, motif) in
                            enumerate(sm.motives_db[cluster_id])]
        self.motives_list = cluster_seqs[:]

        # create PWMs
        super(EdenWrapper, self).fit(motives=cluster_seqs[:])

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
