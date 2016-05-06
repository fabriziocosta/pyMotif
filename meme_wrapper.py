import matplotlib.pyplot as plt
from subprocess import Popen, PIPE
import PIL
import re
import os
from Bio import motifs
#from weblogo_wrapper import Weblogo
from IPython.display import Image, display
from utilities import MuscleAlign_wrapper, Weblogo

class Meme(object):
    """
    Wrapper for MEME 4.11.0
    Usage: meme <sequence file> [options]
    To see MEME help, use MEME.display_meme_help()
    """
    
    def __init__(self,
                 #Output
                 output_dir = "meme_out",
                 text = False,
                 
                 #Alphabet
                 alphabet = "protein", # ["dna", "rna", "protein", file]
                 
                 #Contributing Site Distribution
                 mod = "zoops",
                 
                 #Number of Motifs
                 nmotifs = 1,# 1 is the default value
                 evt = None, #E-value
                 time = None, #CPU time before timeout
                 
                 #Number of Motif Occurences
                 nsites = None, #ignored with OOPS
                 minsites = None, #ignored with OOPS
                 maxsites = None, #ignored with OOPS
                 wnsites = None, 
                 
                 #Motif Width
                 w = None,
                 minw = 8, #default value used by the meme tool
                 maxw = 50, #default value used by the meme tool
                 nomatrim = False, #No Motif Alignment Trim
                 wg = None,
                 ws = None,
                 noendgaps = False,
                 
                 #Markov Background Model
                 bfile = "", #Markov background model file
                 
                 #Position-Specific Priors
                 psp = "", #Position Specific File
                 
                 #Palindromes & Strands
                 revcomp = False,
                 pal = False,
                 
                 #EM Algorithm
                 maxiter = None,
                 distance = None,
                 prior = "", #type of prior to use from {"dirichlet, dmix, mega, megap, addone}
                 b = None,
                 plib = "", #depends on alphabet of sequence
                 
                 #Selecting Starts for EM
                 spfuzz = None, #depends on spmap
                 spmap = "", #depends on alphabet, {uni, pam}
                 cons = None,
                 
                 #Branching Search on EM Starts
                 heapsize = None, 
                 x_branch = False,
                 w_branch = False,
                 bfactor = None, 
                 
                 #Miscellaneous
                 maxsize = None,
                 V = False, #extensive message # 
                 
                 #h = False #display usage message
                 pseudocounts = 0    # alphabet pseudocount for PWM matrix
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
        
        self.cmd_params = ""    # parameters for command string
        self.n_seqs = 0    # no. of seqs in input file, to be set by fit()
        self.seq_names = list()    # to store the names of sequences as given in input file to fit()
        self.motives_db = list()    # list of motives, each represented by an object
        self.motives_list = list()    #list-of-strings representation of motifs
        self.aligned_motives_list = list()    #aligned list-of-strings of motifs, created by display_logo method
        self.widths = list()    # widths given in summary headers; length of each motif
        self.sites = list()    # sites given in summary headers; num occurences
        self.e_values = list()    # e_values given in summary headers
        self.pseudocounts = pseudocounts
        self.pwms = list()    #list of PWM matrices for each motif found
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
        elif self.alphabet == "protein":
            params += " -protein"
        else:
            params += " -alph " + self.alphabet
        
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
        #executes the meme command and returns output
        cmd = "meme "+ fasta_file + " " + self.cmd_params
        io = Popen(cmd.split(" "), stdout=PIPE, stderr=PIPE)
        (stderr, stdout) = io.communicate()
        
        if re.search('error', stdout): raise NameError(stdout.split('\n')[0])
        elif stderr: raise NameError(stdout)
        
        return stdout
    
    
    def fit(self, fasta_file=""):
        if not fasta_file: raise NameError('Input fasta file not specified')
        
        self._make_param_string()
        meme_output = self._command_exec(fasta_file)
        
        # parsing meme output file with Biopython
        filename = os.path.join(self.output_dir, 'meme.txt')
        handle = open(filename)
        record = motifs.parse(handle, 'meme')
        handle.close()
        
        #storing names of sequences given as input to fit()
        self.seq_names = record.sequences[:]
        self.n_seqs = len(record.sequences)
        
        motives_db = list()    #a list of motives, each represented by an object
        for i in range(self.nmotifs):
            motives_db.append(record[i])
        self.motives_db = motives_db[:]
        
        # get string representation of motives
        self.motives_list = list(self._get_motives_list())
        
        # storing length, number of occurences and e-value of each motif
        widths=list()
        sites=list()
        e_values=list()
        
        for i in range(self.nmotifs):
            widths.append(motives_db[i].length)
            sites.append(motives_db[i].num_occurrences)
            e_values.append(motives_db[i].evalue)
            
        self.widths=widths[:]
        self.sites=sites[:]
        self.e_values=e_values[:]
        
        # storing a list of PWMs
        pwm_list=list()
        for i in range(self.nmotifs):
            motif_i = motives_db[i]
            motif_i.pseudocounts = self.pseudocounts
            pwm_list.append(motif_i.pwm)
        self.pwms = pwm_list[:]
        #return output

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
        pvalue_lists = [[[] for j in range(self.nmotifs)] for i in range(len(header))]
        
        for seq_id, seq_name in enumerate(header):    #for all seqs in predict() input
            for i in range(self.nmotifs):    #for all motives found
                if seq_name in motif_occurences[i]:
                    repetitions = motif_occurences[i].count(seq_name)
                    #for multiple occurences of same motif in a sequence
                    
                    if repetitions == 1: seq_index = [motif_occurences[i].index(seq_name)]
                    else: seq_index = [k for k,val in enumerate(motif_occurences[i]) if val==seq_name]
                    
                    #seq_pvalue=list()
                    for j in seq_index:
                        seq_pvalue = self.motives_db[i].instances[j].pvalue
                        pvalue_lists[seq_id][i].append(seq_pvalue)
                else:    #motif[i] does not occur in a sequence
                    pvalue_lists[seq_id][i].append(0.0)
        self.seq_scores = pvalue_lists[:]
        
    
    def _get_match_list(self, header=None, motif_occurences=None):
        match_list = [[[] for j in range(self.nmotifs)] for i in range(len(header))]
        
        for seq_id, seq_name in enumerate(header):    #for all seqs in transform() input
            for i in range(self.nmotifs):    #for all motives found
                if seq_name in motif_occurences[i]:
                    repetitions = motif_occurences[i].count(seq_name)
                    #for multiple occurences of same motif in a sequence
                    
                    if repetitions == 1: seq_index = [motif_occurences[i].index(seq_name)]
                    else: seq_index = [k for k,val in enumerate(motif_occurences[i]) if val==seq_name]
                    
                    for j in seq_index:
                        start_pos = self.motives_db[i].instances[j].start
                        motif_len = self.widths[i]
                        motif_location = zip([start_pos],[start_pos + motif_len])[0]
                        match_list[seq_id][i].append(motif_location)
        return match_list
    
    def predict(self, fasta_file='', return_list=True):
        header = self._get_seq_header_list(fasta_file=fasta_file)
        motif_occurences=self._get_motif_occurences()
        seq_lists = [[] for i in range(len(header))]    #motifs list for every sequence
        
        for seq_id, seq_name in enumerate(header):    #for all seqs in input
            for i in range(self.nmotifs):    #for all motives found
                if seq_name in motif_occurences[i]:
                    for j in range(motif_occurences[i].count(seq_name)):
                        #for multiple occurences of same motif in a seq
                        seq_lists[seq_id].append(i)
                        
        self._make_score_list(header=header, motif_occurences=motif_occurences)
        
        if return_list==True:
            return seq_lists
        return [len(i) for i in seq_lists]
    
    def transform(self, fasta_file='', return_match=True):
        header = self._get_seq_header_list(fasta_file=fasta_file)
        motif_occurences=self._get_motif_occurences()
        
        match_list = self._get_match_list(header=header, motif_occurences=motif_occurences)
        
        if return_match is False:
            match_nums = [[] for i in range(len(header))]
            for i in range(len(header)):
                for j in range(self.nmotifs):
                    if match_list[i][j]:
                        match_nums[i].append(1)
                    else:
                        match_nums[i].append(0)
            return match_nums
        return match_list
    
    def fit_predict(self, fasta_file="", return_list=False):
        output = self.fit(fasta_file=fasta_file)
        return self.predict(fasta_file=fasta_file, return_list=return_list)
    
    def fit_transform(self, fasta_file="", return_match=False):
        output = self.fit(fasta_file=fasta_file)
        return self.transform(fasta_file=fasta_file, return_match=return_match)
    
    def display_meme_help(self):
        cmd = "meme -h"
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
        ma = MuscleAlign_wrapper()
        for i in range(self.nmotifs):
            aligned_motives.append( ma.transform(seqs=motives[i]) )
        
        self.aligned_motives_list = aligned_motives[:]
    
    def _get_logos_list(self, do_alignment=True):
        logos_list = list()
        
        motives=list(self.motives_list)
        if do_alignment is True:
            if not self.aligned_motives_list:
                self.align_motives()
            motives=list(self.aligned_motives_list)
            
        wb = Weblogo(output_format='png', alphabet = 'AGCT', resolution=200, fineprint=' ')
        for i in range(self.nmotifs):
            logo = wb.create_logo(seqs=motives[i])
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
