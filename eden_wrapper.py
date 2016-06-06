# from eden.motif import SequenceMotif


class Eden_wrapper(MotifWrapper):

    """
    Wrapper for EDeN Sequence Motif.
    """

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
