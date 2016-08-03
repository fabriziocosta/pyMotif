"""A eden.sequence_motif_decomposer.SequenceMotifDecomposer wrapper for comparison against Meme."""

from utilities import MotifWrapper, MuscleAlignWrapper
from eden.sequence_motif_decomposer import SequenceMotifDecomposer as SMoD
from sklearn.linear_model import SGDClassifier
from sklearn.cluster import MiniBatchKMeans


class SMoDWrapper(MotifWrapper):
    """Wrapper for EDeN SequenceMotifDecomposer."""

    def __init__(self,
                 alphabet='dna',
                 gap_in_alphabet=True,
                 scoring_criteria='pwm',    # ["pwm","hmm"]
                 pseudocounts=0,
                 threshold=None,    # scoring threshold
                 k=1,

                 complexity=3,
                 n_clusters=10,
                 min_subarray_size=5,
                 max_subarray_size=10,
                 estimator=SGDClassifier(warm_start=True),
                 clusterer=MiniBatchKMeans(),
                 pos_block_size=300,
                 neg_block_size=300,
                 n_jobs=-1,
                 similarity_threshold=0.5,
                 min_score=4,
                 min_freq=0.6,
                 min_cluster_size=10,
                 regex_th=.3,
                 freq_threshold=0.05,

                 muscle_obj=None,
                 weblogo_obj=None
                 ):
        """Initialize a SMoDWrapper object."""
        self.smd = SMoD(complexity=complexity,
                        n_clusters=n_clusters,
                        min_subarray_size=min_subarray_size,
                        max_subarray_size=max_subarray_size,
                        estimator=estimator,
                        clusterer=clusterer,
                        pos_block_size=pos_block_size,
                        neg_block_size=neg_block_size,
                        n_jobs=n_jobs)
        self.similarity_threshold = similarity_threshold
        self.min_score = min_score
        self.min_freq = min_freq
        self.min_cluster_size = min_cluster_size
        self.regex_th = regex_th
        self.freq_threshold = freq_threshold

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
        # list of sequence logos created by WebLogo
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
        """Find motives with SequenceMotifDecomposer."""
        if neg_seqs is None:
            from eden.modifier.seq import seq_to_seq, shuffle_modifier
            neg_seqs = seq_to_seq(seqs, modifier=shuffle_modifier, times=1, order=2)
            neg_seqs = list(neg_seqs)

        self.smd = self.smd.fit(pos_seqs=seqs, neg_seqs=neg_seqs)
        orig_clusters = self.smd.predict(seqs)
        clusters, motives = self.smd.merge(orig_clusters,
                                           similarity_threshold=self.similarity_threshold,
                                           min_score=self.min_score,
                                           min_freq=self.min_freq,
                                           min_cluster_size=self.min_cluster_size,
                                           regex_th=self.regex_th)
        clusters, motives = self.smd.frequency_filter(seqs=seqs,
                                                      motives=motives,
                                                      freq_threshold=self.freq_threshold)

        self.nmotifs = len(clusters.keys())
        if self.nmotifs == 0:
            raise AttributeError("0 motives (clusters) found.")
        self.original_motives_list = self._get_motives_list(clusters)[:]
        self.aligned_motives_list = self._get_aligned_motives_list(
            self.original_motives_list)[:]
        self.motives_list = self.adapt_motives(
            self.aligned_motives_list)[:]

        # create PWMs
        super(SMoDWrapper, self).fit(motives=self.aligned_motives_list)

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
