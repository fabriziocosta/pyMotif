import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
import PIL
import re
import os
from Bio import motifs
from IPython.display import Image, display
from utilities import MuscleAlign_wrapper, Weblogo, PWM

class Meme(PWM):
    """
    Wrapper for MEME 4.11.0
    Usage: meme <sequence file> [options]
    To see MEME help, use MEME.display_meme_help()
    """
    
    def __init__(self,
                 # Output
                 output_dir = "meme_out",
                 text = False,
                 
                 # Alphabet
                 alphabet = "protein", # ["dna", "rna", "protein"]
                 # TODO: meme can also take alphabet file as input
                 
                 # Contributing Site Distribution
                 mod = "zoops",
                 
                 # Number of Motifs
                 nmotifs = 1,    # 1 is the default value
                 evt = None,    # E-value
                 time = None,    # CPU time before timeout
                 
                 # Number of Motif Occurences
                 nsites = None,    # ignored with OOPS
                 minsites = None,    # ignored with OOPS
                 maxsites = None,    # ignored with OOPS
                 wnsites = None, 
                 
                 # Motif Width
                 w = None,
                 minw = 8,    # default value used by the meme tool
                 maxw = 50,    # default value used by the meme tool
                 nomatrim = False,    # No Motif Alignment Trim
                 wg = None,
                 ws = None,
                 noendgaps = False,
                 
                 # Markov Background Model
                 bfile = "",    # Markov background model file
                 
                 # Position-Specific Priors
                 psp = "",    # Position Specific File
                 
                 # Palindromes & Strands
                 revcomp = False,
                 pal = False,
                 
                 # EM Algorithm
                 maxiter = None,
                 distance = None,
                 # type of prior to use from {"dirichlet, dmix, mega, megap, addone}
                 prior = "",
                 b = None,
                 plib = "",    # depends on alphabet of sequence
                 
                 # Selecting Starts for EM
                 spfuzz = None,    # depends on spmap
                 spmap = "",    # depends on alphabet, {uni, pam}
                 cons = None,
                 
                 # Branching Search on EM Starts
                 heapsize = None, 
                 x_branch = False,
                 w_branch = False,
                 bfactor = None, 
                 
                 # Miscellaneous
                 maxsize = None,
                 V = False,    #extensive message
                 
                 #h = False # display usage message

                 # parameters for (parent) PWM
                 pseudocounts = 0,    # alphabet pseudocount

                 # parameters for Muscle Alignment
                 ma_diags = False,
                 ma_maxiters = 16,
                 ma_maxhours = None,

                 #parameters for WebLogo
                 wl_output_format = 'png',  # ['eps', 'png', 'png_print', 'jpeg']
                 wl_stacks_per_line = 40,
                 wl_ignore_lower_case = False,
                 wl_units = 'bits',  # ['bits', 'nats', 'digits', 'kT', 'kJ/mol', 'kcal/mol', 'probability']
                 wl_first_position = 1,
                 wl_logo_range = list(),
                 wl_scale_stack_widths = True,
                 wl_error_bars = True,
                 wl_title = '',
                 wl_figure_label = '',
                 wl_show_x_axis = True,
                 wl_x_label = '',
                 wl_show_y_axis = True,
                 wl_y_label = '',
                 wl_y_axis_tic_spacing = 1.0,
                 wl_show_ends = False,
                 wl_color_scheme = 'classic',  # ['auto', 'base', 'pairing', 'charge', 'chemistry', 'classic', 'monochrome']
                 wl_resolution = 96,
                 wl_fineprint = '',
                ):
        
        self.output_dir = output_dir
        self.text = text
        self.alphabet = alphabet
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
        self.V = V

        self.MA = MuscleAlign_wrapper(diags = ma_diags,
        							  maxiters = ma_maxiters,
        							  maxhours = ma_maxhours,
        							  )

        self.WL = Weblogo(output_format = wl_output_format,
						  stacks_per_line = wl_stacks_per_line,
						  #sequence_type = 
						  ignore_lower_case = wl_ignore_lower_case, 
						  units = wl_units,
						  first_position = wl_first_position,
						  logo_range = wl_logo_range,
						  #composition = 
						  scale_stack_widths = wl_scale_stack_widths,
						  error_bars = wl_error_bars,
						  title = wl_title,
						  figure_label = wl_figure_label,
						  show_x_axis = wl_show_x_axis,
						  x_label = wl_x_label,
						  show_y_axis = wl_show_y_axis,
						  y_label = wl_y_label,
						  y_axis_tic_spacing = wl_y_axis_tic_spacing,
						  show_ends = wl_show_ends,
						  color_scheme = wl_color_scheme,
						  resolution = wl_resolution,
						  fineprint = wl_fineprint,
						  )

        if self.alphabet == "dna":
        	self.WL.sequence_type = "dna"
        elif self.alphabet == "rna":
        	self.WL.sequence_type = "rna"
        elif self.alphabet == "protein":
        	self.WL.sequence_type = "protein"
        else:
        	self.WL.sequence_type = "auto"


        self.cmd_params = ""    # parameters for command string
        self.n_seqs = 0    # no. of seqs in input file, to be set by fit()
        self.seq_names = list()    # to store the names of sequences as given in input file to fit()
        self.motives_db = list()    # list of motives, each represented by an object
        self.motives_list = list()    #list-of-strings representation of motifs
        self.aligned_motives_list = list()    #aligned list-of-strings of motifs, created by display_logo method
        self.widths = list()    # widths given in summary headers; length of each motif
        self.sites = list()    # sites given in summary headers; num occurences
        self.e_values = list()    # e_values given in summary headers
        self.pseudocounts = pseudocounts    # over-rides same attribute of PWM class
        self.seq_scores = list()    #list of p-value of sequences given as input to predict()
        self.logos = list()    #list of sequence logos created with WebLogo
        
    def _make_param_string(self):
        # Creates a string of parameters
        params = "-oc " + self.output_dir
        
        if self.nmotifs!=1:
            params += " -nmotifs " + str(self.nmotifs)
        
        if self.text == True:
            params += " -text"
        
        if self.alphabet == "dna":
            params += " -dna"
        elif self.alphabet == "rna":
            params += " -rna"
        else: # self.alphabet == "protein":
            params += " -protein"
        #else:
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
        
        if self.minsites is not None and self.mod !="oops":
            params += " -minsites " + str(self.minsites)
        
        if self.maxsites is not None and self.mod!="oops":
            params += " -maxsites " + str(self.maxsites)
        
        if self.wnsites is not None:
            params += " -wnsites " + str(self.wnsites)
        
        if self.w is not None:
            params += " -w " + str(self.w)
        elif self.minw!=8 or self.maxw!=50:
            params += " -minw " + str(self.minw)\
                                + " -maxw " + str(self.maxw)
        
        if self.nomatrim == True:
            params += " -nomatrim"
        
        if self.wg is not None:
            params += " -wg " + str(self.wg)
        
        if self.ws is not None:
            params += " -ws " + str(self.ws)
        
        if self.noendgaps != False:
            params += " -noendgaps"
        
        if self.bfile:
            params += " -bfile " + self.bfile
        
        if self.psp:
            params += " -psp " + self.psp
        
        if self.revcomp == True:
            params += " -revcomp"
        
        if self.pal == True:
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
        
        if self.x_branch == True:
            params += " -x_branch"
        
        if self.w_branch == True:
            params += " -w_branch"
        
        if self.bfactor is not None:
            params += " -bfactor " + str(self.bfactor)
        
        if self.maxsize is not None:
            params += " -maxsize " + str(self.maxsize)
        
        if self.V == True:
            params += " -V"
        
        self.cmd_params = params
    
    def _command_exec(self, fasta_file=""):
        # execute the meme command and returns output
        cmd = "meme "+ fasta_file + " " + self.cmd_params
        io = Popen(cmd.split(" "), stdout=PIPE, stderr=PIPE)
        (stderr, stdout) = io.communicate()
        
        if re.search('error', stdout): 
        	raise NameError(stdout.split('\n')[0])
        elif stderr: raise NameError(stdout)
        
        return stdout

    def _get_stats(self):
    	widths = list()
    	sites = list()
    	e_values = list()

    	for i in range(self.nmotifs):
    		widths.append(self.motives_db[i].length)
    		sites.append(self.motives_db[i].num_occurrences)
    		e_values.append(self.motives_db[i].evalue)

    	self.widths = widths[:]
    	self.sites = sites[:]
    	self.e_values = e_values[:]
    
    
    def fit(self, fasta_file=""):
        if not fasta_file:
        	raise NameError('Input fasta file not specified')
        
        self._make_param_string()
        meme_output = self._command_exec(fasta_file)
        
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
        
        # get string representation of motives
        self.motives_list = list(self._get_motives_list())
        
        # store length, number of occurences and e-value of each motif
        self._get_stats()

        # create PWMs
        motives_list = self.motives_list[:]
        super(Meme, self).fit(motives = motives_list)

    def _get_seq_header_list(self, fasta_file=''):

        header=list()
        with open(fasta_file, 'r') as f:
            for line in f.readlines():
                    if line[0]==">":
                        header.append(line.split(">")[1].split(" ")[0])
        return header
    
    def _get_motif_occurences(self):
        #returns a list of lists containing motif occurence sequences' names
        occurences = [[] for i in range(self.nmotifs)]
        
        for i in range(self.nmotifs):
            for j in range(len(self.motives_db[i].instances)):
                occurences[i].append( self.motives_db[i].instances[j].sequence_name)
        return occurences
        
    def _make_score_list(self, header=None, motif_occurences=None):
        pvalue_lists = [
        				[
        				[] for j in range(self.nmotifs)
        				] for i in range(len(header))
        				]
        
        for seq_id, seq_name in enumerate(header):    #for all seqs in predict() input
            for i in range(self.nmotifs):    #for all motives found
                if seq_name in motif_occurences[i]:
                    repetitions = motif_occurences[i].count(seq_name)
                    #for multiple occurences of same motif in a sequence
                    
                    if repetitions == 1: 
                    	seq_index = [motif_occurences[i].index(seq_name)]
                    else: 
                    	seq_index = [k for k,val in enumerate(motif_occurences[i]) if val==seq_name]
                    
                    #seq_pvalue=list()
                    for j in seq_index:
                        seq_pvalue = self.motives_db[i].instances[j].pvalue
                        pvalue_lists[seq_id][i].append(seq_pvalue)
                else:    #motif[i] does not occur in a sequence
                    pvalue_lists[seq_id][i].append(0.0)
        self.seq_scores = pvalue_lists[:]
    
    def _get_match_list(self, header=None, motif_occurences=None):
        match_list = [
        			 [
        			 [] for j in range(self.nmotifs)
        			 ] for i in range(len(header))
        			 ]
        
        for seq_id, seq_name in enumerate(header):    # for all seqs in transform() input
            for i in range(self.nmotifs):    # for all motives found
                if seq_name in motif_occurences[i]:
                    repetitions = motif_occurences[i].count(seq_name)
                    # for multiple occurences of same motif in a sequence
                    
                    if repetitions == 1: 
                    	seq_index = [motif_occurences[i].index(seq_name)]
                    else: 
                    	seq_index = [k for k,val in enumerate(motif_occurences[i]) if val==seq_name]
                    
                    for j in seq_index:
                        start_pos = self.motives_db[i].instances[j].start
                        motif_len = self.widths[i]
                        motif_location = zip([start_pos],[start_pos + motif_len])[0]
                        match_list[seq_id][i].append(motif_location)
        return match_list
    
    def predict(self, 
    			fasta_file = '', 
    			return_list = True,
    			threshold = 1.0e-9,
    			append_score = False):
        header = self._get_seq_header_list(fasta_file=fasta_file)
        motif_occurences=self._get_motif_occurences()
        seq_lists = [[] for i in range(len(header))]    #motifs list for every sequence
        self._make_score_list(header=header,
        					  motif_occurences=motif_occurences)
        
        return super(Meme, self).predict(fasta_file = fasta_file, 
        								return_list = return_list,
        								threshold = threshold,
        								append_score = append_score,
        								)

    def fit_predict(self,
    				fasta_file = "",
    				return_list = False,
    				threshold = 1.0e-9,
    				append_score = False):
        self.fit(fasta_file = fasta_file)
        return self.predict(fasta_file = fasta_file,
        					return_list = return_list,
        					threshold = threshold,
        					append_score = append_score)
    
    def fit_transform(self,
    				  fasta_file = "",
    				  return_match = False,
    				  threshold = 1.0e-9):
        self.fit(fasta_file = fasta_file)
        return self.transform(fasta_file = fasta_file, 
        					  return_match = return_match,
        					  threshold = threshold)
    
    def display_meme_help(self):
        cmd = "meme --help"
        io = Popen(cmd.split(" "), stdout=PIPE, stderr=PIPE)
        (error, output) = io.communicate()
        print output
        
    def _get_motives_list(self):

        #list-of-strings representation of motives
        motives=[]
        for i in range(self.nmotifs):
            motif_i = self.motives_db[i]
            motif_as_str = str(motif_i.instances).split('\n')[:-1]
            headers=[]
            for i in range(len(motif_i.instances)):
                headers.append(motif_i.instances[i].sequence_name)
            motives.append(zip(headers,motif_as_str))
        return motives
        
    def align_motives(self):
        motives=list(self.motives_list)
        aligned_motives = list()
        for i in range(self.nmotifs):
            aligned_motives.append( self.MA.transform(seqs=motives[i]) )
        
        self.aligned_motives_list = aligned_motives[:]
    
    def _get_logos_list(self, do_alignment=True):
        logos_list = list()
        
        motives=list(self.motives_list)
        if do_alignment is True:
            if not self.aligned_motives_list:
                self.align_motives()
            motives=list(self.aligned_motives_list)
            
        for i in range(self.nmotifs):
            logo = self.WL.create_logo(seqs=motives[i])
            logos_list.append(logo)
        return logos_list

    def display_logo(self, motif_num=None, do_alignment=True):
        """Displays logos of all motifs if motif_num is not specified"""
        
        self.logos = self._get_logos_list(do_alignment=do_alignment)[:]

        if motif_num is not None:
            display(Image(self.logos[motif_num-1]))
        else:
            for i in range(self.nmotifs):
                display(Image(self.logos[i]))
