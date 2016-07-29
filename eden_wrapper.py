"""A eden.motif.SequenceMotif wrapper for comparison against Meme."""

from utilities import MotifWrapper, MuscleAlignWrapper
from eden.motif import SequenceMotif


class EdenWrapper(MotifWrapper):
    """Wrapper for EDeN Sequence Motif."""

    def __init__(self,
                 alphabet='dna',
                 gap_in_alphabet=True,
                 scoring_criteria='pwm',    # ["pwm","hmm"]
                 pseudocounts=0,    # integer or dictionary {'A':0, 'C': 0, 'G': 0, 'T': 0}
                 threshold=None,
                 k=1,    # top-k scores returned for hmm score

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

                 muscle_obj=None,
                 weblogo_obj=None
                 ):
        """Initialize a EdenWrapper object."""
        self.sm = SequenceMotif(min_subarray_size=min_subarray_size,
                                max_subarray_size=max_subarray_size,
                                min_motif_count=min_motif_count,
                                min_cluster_size=min_cluster_size,
                                training_size=training_size,
                                negative_ratio=negative_ratio,
                                shuffle_order=shuffle_order,
                                n_iter_search=n_iter_search,
                                complexity=complexity,
                                # radius=radius,
                                # distance=distance,
                                nbits=nbits,
                                clustering_algorithm=clustering_algorithm,
                                n_jobs=n_jobs,
                                n_blocks=n_blocks,
                                block_size=block_size,
                                pre_processor_n_jobs=pre_processor_n_jobs,
                                pre_processor_n_blocks=pre_processor_n_blocks,
                                pre_processor_block_size=pre_processor_block_size,
                                random_state=random_state,
                                )
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

        # Number of motives found
        self.nmotifs = 0
        # over-rides same attribute of MotifWrapper class
        self.pseudocounts = pseudocounts
        # list-of-strings representation of motifs
        self.original_motives_list = list()
        # aligned list-of-strings of motifs;
        # also created by display_logo method
        self.aligned_motives_list = list()
        # adapted motives with no gaps
        self.motives_list = list()
        # list of sequence logos created with WebLogo
        self.logos = list()

    def _get_motives_list(self, db):
        motives = list()
        for i in db.keys():
            motives.append(db[i])
        return motives

    def _get_aligned_motives_list(self, motives):
        aligned_motives = []
        ma = MuscleAlignWrapper()
        for i in range(len(motives)):
            aligned_motives.append(ma.transform(seqs=motives[i]))
        return aligned_motives

    def fit(self, seqs, neg_seqs=None):
        """Find motives with EDeN.SequenceMotif."""
        self.sm.fit(seqs=seqs, neg_seqs=neg_seqs)

        self.nmotifs = len(self.sm.motives_db.keys())
        self.original_motives_list = self._get_motives_list(
            self.sm.motives_db)[:]
        self.aligned_motives_list = self._get_aligned_motives_list(
            self.original_motives_list)[:]
        self.motives_list = self.adapt_motives(
            self.aligned_motives_list)[:]

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
