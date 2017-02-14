# -*- encoding: utf8 -*-

"""
Once a Linguistica object (such as ``lxa_object`` below with the Brown corpus)
is initialized, various methods and attributes are available for automatic
linguistic analysis:

.. code-block:: python

   >>> import linguistica as lxa
   >>> lxa_object = lxa.read_corpus('path/to/english-brown.txt')
   >>> words = lxa_object.wordlist()  # using wordlist()

Basic information
-----------------

.. currentmodule:: linguistica.lexicon.Lexicon

.. autosummary::

   number_of_word_tokens
   number_of_word_types

Word ngrams
-----------

Parameter: ``max_word_tokens``

.. currentmodule:: linguistica.lexicon.Lexicon

.. autosummary::

   wordlist
   word_unigram_counter
   word_bigram_counter
   word_trigram_counter

Morphological signatures
------------------------

Parameters: ``min_stem_length``, ``max_affix_length``, ``min_sig_count``, ``suffixing``

.. currentmodule:: linguistica.lexicon.Lexicon

.. autosummary::

   signatures
   stems
   affixes

   signatures_to_stems
   signatures_to_words
   affixes_to_signatures
   stems_to_signatures
   stems_to_words
   words_in_signatures
   words_to_signatures
   words_to_sigtransforms

Word manifolds and syntactic word neighborhood
----------------------------------------------

Parameters: ``max_word_types``, ``min_context_count``, ``n_neighbors``, ``n_eigenvectors``

.. currentmodule:: linguistica.lexicon.Lexicon

.. autosummary::

   words_to_neighbors
   neighbor_graph
   words_to_contexts
   contexts_to_words

Phonology
---------

.. currentmodule:: linguistica.lexicon.Lexicon

.. autosummary::

   phone_unigram_counter
   phone_bigram_counter
   phone_trigram_counter

Tries
-----

Parameter: ``min_stem_length``

.. currentmodule:: linguistica.lexicon.Lexicon

.. autosummary::

   broken_words_left_to_right
   broken_words_right_to_left
   successors
   predecessors

Other methods and attributes
----------------------------

.. currentmodule:: linguistica.lexicon.Lexicon

.. autosummary::

   parameters
   change_parameters
   use_default_parameters
   reset

"""

#######Importance of _make_all_trie_objects ??

import sys
import os
from io import StringIO

from linguistica import (ngram, signature, manifold, phon, trie)
from linguistica.util import (ENCODING, PARAMETERS, SEP_SIG, SEP_SIGTRANSFORM,
                              double_sorted, fix_punctuations,
                              output_latex, vprint)


"""
FindSuffixesFlag=True	




datafolder 				= "../../data/" + language + "/"
outfolder     			= datafolder    + "lxa/"
infolder 				= datafolder    + 'dx1/'	

graphicsfolder= outfolder + "graphics/"
if not os.path.exists(graphicsfolder):
	os.makedirs(graphicsfolder)

infilename 						= infolder  + shortfilename + ".dx1"
stemfilename 					= infolder  + shortfilename + "_stems.txt"
outfile_Signatures_name 		= outfolder + shortfilename + "_Signatures.txt"  
outfile_SigTransforms_name 		= outfolder + shortfilename + "_sigtransforms.txt"
outfile_SigExtensions_name 		= outfolder + shortfilename + "_sigextensions.txt"
outfile_signature_exemplar_name = outfolder + shortfilename + "_signature_exemplars.txt"
outfile_stemtowords_name 		= outfolder + shortfilename + "_stemtowords.txt"
outfile_FSA_name				= outfolder + shortfilename + "_FSA.txt"
outfile_log_name 				= outfolder + shortfilename + "_log.txt"
outfile_WordToSig_name			= outfolder + shortfilename + "_WordToSig.txt"
outfile_wordparses_name 		= outfolder + shortfilename + "_WordParses.txt"
outfile_wordlist_name 			= outfolder + shortfilename + "_WordList.txt"
outfile_suffix_name 			= outfolder + shortfilename + "_suffixes.txt"
outfile_rebalancing_name 			= outfolder + shortfilename + "_rebalancing_signatures.txt"
#outfile_FSA_graphics_name		= graphicsfolder + shortfilename + "_FSA_graphics.png"
 
print("\n\n" + "-"*100)
if FindSuffixesFlag:
	print("Finding suffixes.")
else:
	print("Finding prefixes.")
print("{:40s}{:>15s}".format("Reading dx file: ", infilename))
print("{:40s}{:>15s}".format("Logging to: ", outfile_log_name))
print("-"* 100)
lxalogfile = open(outfile_log_name, "w")
 

if len(sys.argv) > 1:
	print(sys.argv[1])
	infilename = sys.argv[1] 
if not os.path.isfile(infilename):
	print("Warning: ", infilename, " does not exist.")
if g_encoding == "utf8":
	infile = codecs.open(infilename, g_encoding = 'utf-8')
else:
	infile = open(infilename) 


#----------------------------------------------------------#
 
 
 
if g_encoding == "utf8":
	print("yes utf8")
else:
	Signatures_outfile = open (outfile_Signatures_name, mode = 'w')

outfile_Signatures		= open (outfile_Signatures_name,"w")
outfile_FSA 			= open (outfile_FSA_name,"w")
outfile_SigExemplars 	= open (outfile_signature_exemplar_name,"w")
outfile_WordToSig		= open (outfile_WordToSig_name,"w")
outfile_SigTransforms   = open (outfile_SigTransforms_name,"w")
outfile_StemToWords   	= open (outfile_stemtowords_name,"w") 
outfile_WordParses   	= open (outfile_wordparses_name,"w") 
outfile_WordList   		= open (outfile_wordlist_name,"w") 
outfile_SigExtensions   = open (outfile_SigExtensions_name,"w") 
outfile_Suffixes   		= open (outfile_suffix_name,"w") 
outfile_Rebalancing_Signatures    		= open (outfile_rebalancing_name,"w") 

#outfile_FSA_graphics   		= open (outfile_FSA_graphics_name,"w")  
#-------------------------------------------------------------------------------------------------------# 
#-------------------------------------------------------------------------------------------------------#
# 					Main part of program		   			   	#
#-------------------------------------------------------------------------------------------------------#
#-------------------------------------------------------------------------------------------------------#


	# This is just part of documentation:
	# A signature is a tuple of strings (each an affix).
	# Signatures is a map: its keys are signatures.  Its values are *sets* of stems. 
	# StemToWord is a map; its keys are stems.       Its values are *sets* of words.
	# StemToSig  is a map; its keys are stems.       Its values are individual signatures.
	# WordToSig  is a Map. its keys are words.       Its values are *lists* of signatures.
	# StemCounts is a map. Its keys are words. 	 Its values are corpus counts of stems.
	# SignatureToStems is a dict: its keys are tuples of strings, and its values are dicts of stems. We don't need both this and Signatures!

 

Lexicon = CLexicon( )
 
#--------------------------------------------------------------------##
#		read wordlist (dx1)
#--------------------------------------------------------------------##

filelines= infile.readlines()
 
 #move this function into jackson's code and get it to print out a wordlist
 #it should put the data into the data structure that Jackson sets up
 #get a dx1 file from the github page that Jackson's set up (under Data)

#this is for a dx1 file:
for line in filelines:
	#print line
	pieces = line.split()	
	word=pieces[0] 	
	if word == '#' :
		continue
	if word.isalpha() == False:
		continue
	if len(pieces) > 1:
		count = int(pieces[1])
	else:
		count =1 
	word = word.lower()
	if (BreakAtHyphensFlag and '-' in word):
		words = word.split('-')
		for word in words:
			if word.isalpha() == False:
				continue
			if word not in WordCounts:
				WordCounts[word]=0
			Lexicon.WordCounts[word]+=count
			Lexicon.TotalLetterCountInWords  += len(word)
	else:	
		if word not in Lexicon.WordCounts:
			Lexicon.WordCounts[word]=0 
		Lexicon.WordCounts[word]  = count
		Lexicon.TotalLetterCountInWords += len(word)
	if len(Lexicon.WordCounts) >= numberofwords:
		break
 
	Lexicon.WordList.AddWord(word) 
 
Lexicon.WordList.sort()
print("\n1. Finished reading word list.\n")
Lexicon.PrintWordList(outfile_WordList)

 

print(>>outfile_Signatures, "# ", language, numberofwords)



print(>>lxalogfile, "{:40s}{:>15s}".format("Language: ", language))
print(>>lxalogfile, "{:40s}{:10,d}".format("Total words:", len(Lexicon.WordList.mylist)))
print(>>lxalogfile, "{:40s}{:>10,}".format("Minimum Stem Length", Lexicon.MinimumStemLength))
print(>>lxalogfile, "{:40s}{:>10,}".format("Maximum Affix Length", Lexicon.MaximumAffixLength ))
print(>>lxalogfile, "{:40s}{:>10,}".format("Minimum Number of stems in signature: ", Lexicon.MinimumStemsInaSignature))
print(>>lxalogfile, "{:40s}{:10,d}".format("Total letter count in words: ", Lexicon.TotalLetterCountInWords))
print(>>lxalogfile, "{:40s}{:10.2f}".format("Average letters per word: ", float(Lexicon.TotalLetterCountInWords)/len(Lexicon.WordList.mylist)))

print("{:40s}{:>15s}".format("Language: ", language)
print("{:40s}{:10,d}".format("Total words:", len(Lexicon.WordList.mylist)))
print("{:40s}{:>10,}".format("Minimum Stem Length", Lexicon.MinimumStemLength))
print("{:40s}{:>10,}".format("Maximum Affix Length", Lexicon.MaximumAffixLength ))
print("{:40s}{:>10,}".format("Minimum Number of stems in signature: ", Lexicon.MinimumStemsInaSignature))
print("{:40s}{:10,d}".format("Total letter count in words: ", Lexicon.TotalLetterCountInWords))
print("{:40s}{:10.2f}".format("Average letters per word: ", float(Lexicon.TotalLetterCountInWords)/len(Lexicon.WordList.mylist)))
 
 

splitEndState = True
morphology= FSA_lxa(splitEndState)


if True: 	
	print ("2. Make Signatures.")
	MakeSignatures(Lexicon, lxalogfile,outfile_Rebalancing_Signatures,FindSuffixesFlag,Lexicon.MinimumStemLength)
	
 
if True:
	print ("3. Printing signatures.")
	printSignatures(Lexicon, lxalogfile, outfile_Signatures, outfile_WordToSig, outfile_StemToWords, outfile_SigExtensions,outfile_Suffixes ,g_encoding, FindSuffixesFlag)
 
if False:
	print ("4. Printing signature transforms for each word.")
	printWordsToSigTransforms(Lexicon.SignatureToStems, Lexicon.WordToSig, Lexicon.StemCounts, outfile_SigTransforms, g_encoding, FindSuffixesFlag)
 
if False:
	print ("5. Slicing signatures.")
	SliceSignatures(Lexicon,  g_encoding, FindSuffixesFlag, lxalogfile	)

if False:
	print ("6. Adding signatures to the FSA.")
	AddSignaturesToFSA(Lexicon, Lexicon.SignatureToStems, morphology,FindSuffixesFlag) 

if False:
	print ("7. Printing the FSA.")
	print (>>outfile_FSA, "#", language, shortfilename, numberofwords)
	morphology.printFSA(outfile_FSA) 
 
if False :
	print ("8. Printing signatures.")
	printSignatures(Lexicon, outfile_Signatures, outfile_WordToSig, outfile_StemToWords,outfile_Suffixes, g_encoding, FindSuffixesFlag,letterCostOfStems, letterCostOfSignatures)

if False:
	print ("9. Finding robust peripheral pieces on edges in FSA.")
	for loopno in range( NumberOfCorrections):
		morphology.find_highest_weight_affix_in_an_edge ( lxalogfile, FindSuffixesFlag)

if False:
	print ("10. Printing graphs of the FSA.")
	for state in morphology.States:	
		# do not print the entire graph:
		#if state == morphology.startState:
		#	continue
		###
		graph = morphology.createPySubgraph(state) 
		# if the graph has 3 edges or fewer, do not print it:	

	 	if len(graph.edges()) < 4:
	 		continue
 
 
		graph.layout(prog='dot')
		filename = graphicsfolder + 'morphology' + str(state.index) + '.png'
		graph.draw(filename) 
		if (False):
			filename = graphicsfolder + 'morphology' + str(state.index) + '.dot'
			graph.write(filename)
  	
#---------------------------------------------------------------------------------#	
#	5d. Print FSA again, with these changes.
#---------------------------------------------------------------------------------# 

if False:
	print "11. Printing the FSA."
	morphology.printFSA(outfile_FSA)
 
if False:
	print "12. Parsing all words through FSA."
	morphology.parseWords(Lexicon.WordToSig.keys(), outfile_WordParses)
	
if False:	
	print "13. Printing all the words' parses."
	morphology.printAllParses(outfile_WordParses)


if False:
	print "14. Now look for common morphemes across different edges."
	morphology.findCommonMorphemes(lxalogfile)

if False:
 
	# we will remove this. It will be replaced by a function that looks at all cross-edge sharing of morphemes. 
	print >>outfile_FSA, "Finding common stems across edges."
	HowManyTimesToCollapseEdges = 0
	for loop in range(HowManyTimesToCollapseEdges): 
	 	print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		print  "Loop number", loop
	 	print "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
		(commonEdgePairs,  EdgeToEdgeCommonMorphs) = morphology.findCommonStems(lxalogfile)
	 
		if len( commonEdgePairs ) == 0:
			print "There are no more pairs of edges to consider."
			break
		edge1, edge2 = commonEdgePairs[0]
		state1 = edge1.fromState
		state2 = edge2.fromState
		state3 = edge1.toState
		state4 = edge2.toState
		print "\n\nWe are considering merging edge ", edge1.index,"(", edge1.fromState.index, "->", edge1.toState.index, ") and  edge", edge2.index, "(", edge2.fromState.index, "->", edge2.toState.index , ")"
		 
		print "Printed graph", str(loop), "before_merger"
		graph = morphology.createDoublePySubgraph(state1,state2) 	
		graph.layout(prog='dot')
		filename = graphicsfolder +  str(loop) + '_before_merger' + str(state1.index) + "-" + str(state2.index) + '.png'
		graph.draw(filename) 

		if state1 == state2:
			print "The from-States are identical"
			state_changed_1 = state1
			state_changed_2 = state2
			morphology.mergeTwoStatesCommonMother(state3,state4)
			morphology.EdgePairsToIgnore.append((edge1, edge2))

		elif state3 == state4:
			print "The to-States are identical"
			state_changed_1 = state3
			state_changed_2 = state4	 
			morphology.mergeTwoStatesCommonDaughter(state1,state2) 
			morphology.EdgePairsToIgnore.append((edge1, edge2))

		elif morphology.mergeTwoStatesCommonMother(state1,state2):
			print "Now we have merged two sister edges from line 374 **********"
			state_changed_1 = state1
			state_changed_2 = state2
			morphology.EdgePairsToIgnore.append((edge1, edge2))

	
		elif   morphology.mergeTwoStatesCommonDaughter((state3,state4))  : 
			print "Now we have merged two daughter edges from line 377 **********"
			state_changed_1 = state3
			state_changed_2 = state4
			morphology.EdgePairsToIgnore.append((edge1, edge2))
			 
		graph = morphology.createPySubgraph(state1) 	
		graph.layout(prog='dot')
		filename = graphicsfolder + str(loop) +  '_after_merger_' + str(state_changed_1.index) +  "-" + str(state_changed_2.index) + '.png'
		print "Printed graph", str(loop), "after_merger"
		graph.draw(filename) 
 

#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#
#		User inquiries about morphology
#------------------------------------------------------------------------------------------#
#------------------------------------------------------------------------------------------#

morphology_copy = morphology.MakeCopy()


initialParseChain = list()
CompletedParses = list()
IncompleteParses = list()
word = "" 
while False:
	word = raw_input('Inquire about a word: ')
	if word == "exit":
		break
	if word == "State":
		while True:
			stateno = raw_input("State number:")
			if stateno == "" or stateno == "exit":
				break
			stateno = int(stateno)	
			for state in morphology.States:
				if state.index == stateno:
					break	
			state = morphology.States[stateno]
			for edge in state.getOutgoingEdges():
				print "Edge number", edge.index 
				i = 0
				for morph in edge.labels:
					print "%12s" % morph,
					i+=1
					if i%6 == 0: print 
			print "\n\n"		
			continue
	if word == "Edge":
		while True:
			edgeno = raw_input("Edge number:")
			if edgeno == "" or edgeno == "exit":
				break
			edgeno = int(edgeno)
			for edge in morphology.Edges:
				if edge.index == edgeno:
					break
			print "From state", morphology.Edges[edgeno].fromState.index, "To state", morphology.Edges[edgeno].toState.index
			for edge in morphology.Edges:
				if edge.index == int(edgeno):
					morphlist = list(edge.labels)
			for i in range(len( morphlist )):
				print "%12s" % morphlist[i],
				if i%6 == 0:
					print	
			print "\n\n"
			continue
	if word == "graph":
		while True:
			stateno = raw_input("Graph state number:")
			
	del CompletedParses[:]
	del IncompleteParses[:]
	del initialParseChain[:]
	startingParseChunk = parseChunk("", word)
	startingParseChunk.toState = morphology.startState

	initialParseChain.append(startingParseChunk)
	IncompleteParses.append(initialParseChain)
	while len(IncompleteParses) > 0 :
		CompletedParses, IncompleteParses = morphology.lparse(CompletedParses, IncompleteParses)
	if len(CompletedParses) == 0: print "no analysis found." 
	 
	for parseChain in CompletedParses:
		for thisParseChunk in  parseChain:			
			if (thisParseChunk.edge):				 
				print "\t",thisParseChunk.morph,  
		print 
	print

	for parseChain in CompletedParses:
		print "\tStates: ",
		for thisParseChunk in  parseChain:			
			if (thisParseChunk.edge):				 
				print "\t",thisParseChunk.fromState.index, 
		print "\t",thisParseChunk.toState.index 	 
	print 

	for parseChain in CompletedParses:
		print "\tEdges: ",
		for thisParseChunk in  parseChain:			
			if (thisParseChunk.edge):				 
				print "\t",thisParseChunk.edge.index,
		print
	print "\n\n"



#---------------------------------------------------------------------------------------------------------------------------#  
#---------------------------------------------------------------------------------#	
#	Close output files
#---------------------------------------------------------------------------------# 
  
outfile_FSA.close()
Signatures_outfile.close() 
outfile_SigTransforms.close() 
lxalogfile.close()
outfile_SigExtensions.close()
outfile_Suffixes.close()
#---------------------------------------------------------------------------------#	
#	Logging information
#---------------------------------------------------------------------------------# 

localtime = time.asctime( time.localtime(time.time()) )
print "Local current time :", localtime

 
lxalogfile.close()
#--------------------------------------------------#

"""

def getrobustness(sig, stems):
	# ----------------------------------------------------------------------------------------------------------------------------#
	countofsig = len(sig)
	countofstems = len(stems)
	lettersinstems = 0
	lettersinaffixes = 0
	for stem in stems:
		lettersinstems += len(stem)
	for affix in sig:
		lettersinaffixes += len(affix)
	# ----------------------------------------------------------------------------------------------------------------------------#
	return lettersinstems * (countofsig - 1) + lettersinaffixes * (countofstems - 1)

class CWordList:
    def __init__(self):
        self.mylist = list()
       
    def GetCount(self):
        return len(self.mylist)
    def AddWord(self, word):
        self.mylist.append(CWord(word))

    def at(self, n):
        return self.mylist[n]

    def sort(self):
        self.mylist.sort(key=lambda item: item.Key)
        # for item in self.mylist:
        #	print item.Key
        for i in range(len(self.mylist)):
            word = self.mylist[i]
            word.leftindex = i
        templist = list()
        for word in self.mylist:
            thispair = (word.Key[::-1], word.leftindex)
            templist.append(thispair)
        templist.sort(key=lambda item: item[0])
        for i in range(len(self.mylist)):
            (drow, leftindex) = templist[i]
            self.mylist[leftindex].rightindex = i

    def PrintWordList(self, outfile):
        Size = float(len(self.mylist))
        for word in self.mylist:
            x = word.leftindex / Size
            y = word.rightindex / Size
            print >> outfile, "{:20s} {:9.5} {:9.5}".format(word.Key, x, y)


class CWord:
    def __init__(self, key):
        self.Key = key
        self.leftindex = -1
        self.rightindex = -1

    def makeword(stem, affix, sideflag):
        if sideflag == True:
                return stem + affix
        else:
                return affix + stem

    def byWordKey(word):
        return word.Key

class Lexicon:
    """
    A class for a Linguistica object. It is called "Lexicon" for the historical
    reason that the same element in the C++ version of Linguistica 4 is also
    called as such.
    """

    def __init__(self, file_path=None, wordlist_file=False, corpus_object=None,
                 wordlist_object=None, encoding=ENCODING, **kwargs):
        self.file_abspath = self._check_file_path(file_path)

        if self.file_abspath is None:
            self.directory = None
        else:
            self.directory = os.path.dirname(self.file_abspath)

        self.file_is_wordlist = wordlist_file
        self.encoding = encoding
        self.corpus_object = corpus_object
        self.wordlist_object = wordlist_object
        self.parameters_ = self._determine_parameters(**kwargs)

        self._initialize()

    @staticmethod
    def _check_file_path(file_path):
        """
        Return the absolute path of *file_path*.
        """
        if file_path is None:
            return None

        if type(file_path) != str:
            raise TypeError('file path must be a str -- ' + file_path)

        if sys.platform.startswith('win'):
            file_path = file_path.replace('/', os.sep)
        else:
            file_path = file_path.replace('\\', os.sep)

        file_abspath = os.path.abspath(file_path)
        if not os.path.isfile(file_abspath):
            raise FileNotFoundError(file_path)
        else:
            return file_abspath

    @staticmethod
    def _determine_parameters(**kwargs):
        """
        Determine the parameter dict.
        """
        temp_parameters = dict(PARAMETERS)

        for parameter in kwargs.keys():
            if parameter not in PARAMETERS:
                raise KeyError('unknown parameter -- ' + parameter)
            else:
                temp_parameters[parameter] = kwargs[parameter]

        return temp_parameters

    def PrintWordList(self, outfile):
        self._wordlist.PrintWordList(outfile)

    def parameters(self):
        """
        Return the parameter dict.

        :rtype: dict(str: int)
        """
        return self.parameters_

    def change_parameters(self, **kwargs):
        """
        Change parameters specified by *kwargs*.

        :param kwargs: keyword arguments for parameters and their new values
        """
        for parameter, new_value in kwargs.items():
            if parameter not in self.parameters_:
                raise KeyError('unknown parameter -- ' + parameter)

            self.parameters_[parameter] = new_value

    def use_default_parameters(self):
        """
        Reset parameters to their default values.
        """
        self.parameters_ = dict(PARAMETERS)

    def _initialize(self):
        # number of word types and tokens
        self._number_of_word_types = None
        self._number_of_word_tokens = None

        # word ngrams
        self._word_unigram_counter = None #WordCounts
        self._word_bigram_counter = None
        self._word_trigram_counter = None

        # wordlist
        self._wordlist = CWordList()
        """if self.wordlist_object is not None:
            # self.wordlist_object is
            # either an iterable or a dict of word-count pairs

            # 11/7 now only for a dict
            
            if type(self.wordlist_object) is dict:
                word_count_dict = dict()
                if self.parameters_['keep_case']:
                    word_count_dict = self.wordlist_object
                else:
                    for word, count in self.wordlist_object:
                        word = word.lower()
                        if word not in word_count_dict:
                            word_count_dict[word] = 0
                        word_count_dict[word] += count

                #self._wordlist = [word for word, _ in
                                  #double_sorted(word_count_dict.items(),
                                                key=lambda x: x[1],
                                                reverse=True)]
                self._word_unigram_counter = word_count_dict

            elif hasattr(self.wordlist_object, '__iter__'):
                if self.parameters_['keep_case']:
                    self._wordlist = sorted(set(self.wordlist_object))
                else:
                    self._wordlist = sorted(
                        set(w.lower() for w in self.wordlist_object))
                self._word_unigram_counter = {w: 1 for w in self._wordlist}

            else:
                raise TypeError('wordlist object must be a dict of word-count'
                                'pairs or an iterable of words')
"""
        # signature-related objects
        self._stems_to_words = {} #StemToWord
        self._signatures_to_stems = {} #SignaturetoStem
        self._stems_to_signatures = {} #StemToSignature
        self._words_to_signatures = {} #WordToSig
        self._signatures_to_words = {}
        self._words_to_sigtransforms = {}

        self._signatures = {} #Signatures
        self._affixes_to_signatures = {}
        self._words_in_signatures = {}
        self._affixes = {} #StemToAffix
        self._stems = {}
        self._suffixes = {} #Suffixes
        self._prefixes = {} #Prefixes
        self.StemCounts = {}

# don't distinguish between a corpus file and a wordlist file
        # corpus file object
        if self.corpus_object is not None:
            # self.corpus_object is either a list of strings or a long str
            if type(self.corpus_object) is list:
                corpus_str = fix_punctuations(' '.join(self.corpus_object))
            elif type(self.corpus_object) is str:
                corpus_str = fix_punctuations(self.corpus_object)
            else:
                raise TypeError('corpus object must be either a str or a list')
            self.corpus_file_object = StringIO(corpus_str)
        elif self.file_abspath and not self.file_is_wordlist:
            self.corpus_file_object = open(self.file_abspath,
                                           encoding=self.encoding)
        else:
            self.corpus_file_object = None

        # wordlist file object
        if self.file_is_wordlist:
            self.wordlist_file_object = open(self.file_abspath,
                                             encoding=self.encoding)
        else:
            self.wordlist_file_object = StringIO()

        self._wordlist = CWordList()

        #set numbers from ClassLexicon file
        self.MinimumStemsInaSignature =2 #means that a signature needs more than 1 stem
        self.MinimumStemLength = 5 
        self.MaximumAffixLength =3
        self.MaximumNumberOfAffixesInASignature = 10
        self.NumberOfAnalyzedWords = 0 #increases as you go through the wordlist
        self.LettersInAnalyzedWords = 0 #total number of letters
        self.NumberOfUnanalyzedWords = 0 
        self.LettersInUnanalyzedWords = 0
        self.TotalLetterCountInWords = 0
        self.LettersInStems = 0 #overall sum of the letters in the stems
        self.AffixLettersInSignatures = 0 #count NULL as 1

        self.TotalRobustInSignatures = 0 #robustmess is how many letters you're saving by putting the signature together


        # manifold-related objects
        self._words_to_neighbors = None
        self._words_to_contexts = None
        self._contexts_to_words = None
        self._neighbor_graph = None

        # phon objects
        self._phone_unigram_counter = None
        self._phone_bigram_counter = None
        self._phone_trigram_counter = None

        self._phone_dict = None
        self._biphone_dict = None
        self._word_dict = None
        self._words_to_phones = None

        # trie objects
        self._broken_words_left_to_right = {}
        self._broken_words_right_to_left = {}
        self._successors = {}
        self._predecessors = {}

    def sort_wordlist(self):
        self._wordlist.sort()

    def reset(self):
        """
        Reset the Linguistica object. While the file path information is
        retained, all computed objects (ngrams, signatures, word neighbors, etc)
        are reset to ``NULL``; if they are called again, they are re-computed.
        """
        self._initialize()

    def run_all_modules(self, verbose=False):
        """
        Run all modules.
        """
        self.run_ngram_module(verbose=verbose)
        self.run_phon_module(verbose=verbose)
        self.run_signature_module(verbose=verbose)
        self.run_trie_module(verbose=verbose)

        if self.corpus_file_object:
            self.run_manifold_module(verbose=verbose)

    def ComputeRobustness(self):
        self.NumberOfAnalyzedWords= len(self._words_to_signatures)
        self.NumberOfUnanalyzedWords= len(self._wordlist.mylist) - self.NumberOfAnalyzedWords
        for sig in self._signatures_to_stems:
                numberofaffixes = len(sig)
                mystems = self._signatures_to_stems[sig]
                numberofstems = len(mystems)
                AffixListLetterLength = 0
                for affix in sig:
                        if affix == "NULL":
                                continue
                        AffixListLetterLength += len(affix)
                StemListLetterLength = 0
                for stem in mystems:
                        StemListLetterLength += len(stem)

                self.TotalRobustnessInSignatures +=  getrobustness(mystems,sig)
                self.AffixLettersInSignatures += AffixListLetterLength


    def StemsToSignatures(self,FindSuffixesFlag):
        if FindSuffixesFlag:
                Affixes = self._suffixes
        else:
                Affixes = self._prefixes
        for stem in self._stems_to_words:
                self.LettersInStems += len(stem)
                signature = list(self._affixes[stem])
                signature.sort()
                signature_tuple = tuple(signature)
                #print "A", stem, signature_tuple
                for affix in signature:
                        if affix not in Affixes:
                                Affixes[affix] = 0
                                Affixes[affix] += 1
                        if signature_tuple not in self._signatures_to_stems:
                                self._signatures_to_stems[signature_tuple] = dict()
                                for affix in signature:
                                        self.TotalLetterCostOfAffixesInSignatures += len(affix)
                        self._signatures_to_stems[signature_tuple][stem] = 1
                        self._stems_to_signatures[stem] = signature_tuple
                        for word in self._stems_to_words[stem]:
                                if word not in self._words_to_signatures:
                                        self._words_to_signatures[word] = list()
                                self._words_to_signautres[word].append(signature_tuple)
                for sig in self._signatures_to_stems:
                        if len(self._signatures_to_stems[sig]) < self.MinimumStemsInaSignature:
                                for stem in self._signatures_to_stems[sig]:
                                        del self._stems_to_signatures[stem]
                                        for word in self._stems_to_words[stem]:
                                                if len( self._words_to_signatures[word] ) == 1:
                                                        del self._words_to_signatures[word]
                                                else:
                                                        self._words_to_signatures[word].remove(sig)
                                        del self._stems_to_words[stem]


    def output_all_results(self, directory=None, verbose=False, test=False):
        """
        Output all Linguistica results to *directory*.

        :param directory: output directory. If not specified, it defaults to
            the current directory given by ``os.getcwd()``.
        """
        if not directory:
            output_dir = os.getcwd()
        else:
            output_dir = os.path.abspath(directory)

        # ----------------------------------------------------------------------
        if self.corpus_file_object:
            vprint('ngram objects', verbose=verbose)
            fname = 'word_bigrams.txt'
            obj = double_sorted(self.word_bigram_counter().items(),
                                key=lambda x: x[1], reverse=True)
            f_path = os.path.join(output_dir, fname)
            print (output_dir, "is the output directory.")
            output_latex(obj, f_path,
                         title='Word bigrams',
                         headers=['Word bigram', 'Count'],
                         row_functions=[lambda x: ' '.join(x[0]),
                                        lambda x: x[1]],
                         column_widths=[50, 10],
                         lxa_parameters=self.parameters(),
                         test=test, encoding=self.encoding,
                         number_of_word_types=self.number_of_word_types(),
                         number_of_word_tokens=self.number_of_word_tokens(),
                         input_file_path=self.file_abspath)
            vprint('\t' + fname, verbose=verbose)

            fname = 'word_trigrams.txt'
            obj = double_sorted(self.word_trigram_counter().items(),
                                key=lambda x: x[1], reverse=True)
            f_path = os.path.join(output_dir, fname)
            output_latex(obj, f_path,
                         title='Word trigrams',
                         headers=['Word trigram', 'Count'],
                         row_functions=[lambda x: ' '.join(x[0]),
                                        lambda x: x[1]],
                         column_widths=[75, 10],
                         lxa_parameters=self.parameters(),
                         test=test, encoding=self.encoding,
                         number_of_word_types=self.number_of_word_types(),
                         number_of_word_tokens=self.number_of_word_tokens(),
                         input_file_path=self.file_abspath)
            vprint('\t' + fname, verbose=verbose)

        # ----------------------------------------------------------------------
        vprint('morphological signature objects', verbose=verbose)
        fname = 'stems_to_words.txt'
        obj = double_sorted(self.stems_to_words().items(),
                            key=lambda x: len(x[1]), reverse=True)
        print(obj)
        f_path = os.path.join(output_dir, fname)
        output_latex(obj, f_path,
                     title='Stems to words '
                           '(descending order of word count)',
                     headers=['Stem', 'Word count', 'Words'],
                     row_functions=[lambda x: x[0],
                                    lambda x: len(x[1]),
                                    lambda x: ', '.join(sorted(x[1]))],
                     column_widths=[15, 15, 0],
                     lxa_parameters=self.parameters(),
                     test=test, encoding=self.encoding,
                     number_of_word_types=self.number_of_word_types(),
                     number_of_word_tokens=self.number_of_word_tokens(),
                     input_file_path=self.file_abspath)
        vprint('\t' + fname, verbose=verbose)

        fname = 'stems_to_words.txt'
        obj = double_sorted(self.stems_to_words().items(),
                            key=lambda x: x[0], reverse=False)
        f_path = os.path.join(output_dir, fname)
        output_latex(obj, f_path,
                     title='Stems to words '
                           '(alphabetical order of stems)',
                     headers=['Stem', 'Word count', '1st 10 words'],
                     row_functions=[lambda x: x[0],
                                    lambda x: len(x[1]),
                                    lambda x: ', '.join(sorted(x[1]))],
                     column_widths=[15, 15, 0],
                     lxa_parameters=self.parameters(),
                     test=test, encoding=self.encoding,
                     number_of_word_types=self.number_of_word_types(),
                     number_of_word_tokens=self.number_of_word_tokens(),
                     input_file_path=self.file_abspath)
        vprint('\t' + fname, verbose=verbose)

        fname = 'signatures_to_stems.txt'
        obj = double_sorted(self.signatures_to_stems().items(),
                            key=lambda x: len(x[1]), reverse=True)
        f_path = os.path.join(output_dir, fname)
        output_latex(obj, f_path,
                     title='Signatures to stems',
                     headers=['Signature', 'Stem count', 'Stems'],
                     row_functions=[lambda x: SEP_SIG.join(x[0]),
                                    lambda x: len(x[1]),
                                    lambda x: ', '.join(sorted(x[1]))],
                     column_widths=[30, 15, 0],
                     lxa_parameters=self.parameters(),
                     test=test, encoding=self.encoding,
                     number_of_word_types=self.number_of_word_types(),
                     number_of_word_tokens=self.number_of_word_tokens(),
                     input_file_path=self.file_abspath)
        vprint('\t' + fname, verbose=verbose)

        fname = 'signatures_to_stems_truncated.txt'
        obj = double_sorted(self.signatures_to_stems().items(),
                            key=lambda x: len(x[1]), reverse=True)
        f_path = os.path.join(output_dir, fname)
        output_latex(obj, f_path,
                     title='Signatures to stems '
                           '(first 10 stems for each sig)',
                     headers=['Signature', 'Stem count', '1st 10 stems'],
                     row_functions=[lambda x: SEP_SIG.join(x[0]),
                                    lambda x: len(x[1]),
                                    lambda x:
                                    ' '.join(sorted(x[1])[:10])],
                     column_widths=[30, 15, 0],
                     lxa_parameters=self.parameters(),
                     test=test, encoding=self.encoding,
                     number_of_word_types=self.number_of_word_types(),
                     number_of_word_tokens=self.number_of_word_tokens(),
                     input_file_path=self.file_abspath)
        vprint('\t' + fname, verbose=verbose)

        fname = 'stems_to_signatures.txt'
        obj = double_sorted(self.stems_to_signatures().items(),
                            key=lambda x: len(x[1]), reverse=True)
        f_path = os.path.join(output_dir, fname)
        output_latex(obj, f_path,
                     title='Stems to signatures',
                     headers=['Stems', 'Signatures'],
                     row_functions=[lambda x: x[0],
                                    lambda x:
                                    ', '.join(SEP_SIG.join(sig)
                                              for sig in sorted(x[1]))],
                     column_widths=[15, 0],
                     lxa_parameters=self.parameters(),
                     test=test, encoding=self.encoding,
                     number_of_word_types=self.number_of_word_types(),
                     number_of_word_tokens=self.number_of_word_tokens(),
                     input_file_path=self.file_abspath)
        vprint('\t' + fname, verbose=verbose)

        fname = 'words_to_signatures.txt'
        obj = double_sorted(self.words_to_signatures().items(),
                            key=lambda x: len(x[1]), reverse=True)
        f_path = os.path.join(output_dir, fname)
        output_latex(obj, f_path,
                     title='Words to signatures',
                     headers=['Word', 'Sig count', 'Signatures'],
                     row_functions=[lambda x: x[0],
                                    lambda x: len(x[1]),
                                    lambda x:
                                    ', '.join(SEP_SIG.join(sig)
                                              for sig in sorted(x[1]))],
                     column_widths=[25, 15, 0],
                     lxa_parameters=self.parameters(),
                     test=test, encoding=self.encoding,
                     number_of_word_types=self.number_of_word_types(),
                     number_of_word_tokens=self.number_of_word_tokens(),
                     input_file_path=self.file_abspath)
        vprint('\t' + fname, verbose=verbose)

        fname = 'signatures_to_words.txt'
        obj = double_sorted(self.signatures_to_words().items(),
                            key=lambda x: len(x[1]), reverse=True)
        f_path = os.path.join(output_dir, fname)
        output_latex(obj, f_path,
                     title='Signatures to words',
                     headers=['Signature', 'Word count', 'Words'],
                     row_functions=[lambda x: SEP_SIG.join(x[0]),
                                    lambda x: len(x[1]),
                                    lambda x: ', '.join(sorted(x[1]))],
                     column_widths=[20, 15, 0],
                     lxa_parameters=self.parameters(),
                     test=test, encoding=self.encoding,
                     number_of_word_types=self.number_of_word_types(),
                     number_of_word_tokens=self.number_of_word_tokens(),
                     input_file_path=self.file_abspath)
        vprint('\t' + fname, verbose=verbose)

        fname = 'signatures_to_words_truncated.txt'
        obj = double_sorted(self.signatures_to_words().items(),
                            key=lambda x: len(x[1]), reverse=True)
        f_path = os.path.join(output_dir, fname)
        output_latex(obj, f_path,
                     title='Signatures to words '
                           '(first 10 words for each sig)',
                     headers=['Signature', 'Word count', '1st 10 words'],
                     row_functions=[lambda x: SEP_SIG.join(x[0]),
                                    lambda x: len(x[1]),
                                    lambda x:
                                    ', '.join(sorted(x[1])[:10])],
                     column_widths=[20, 15, 0],
                     lxa_parameters=self.parameters(),
                     test=test, encoding=self.encoding,
                     number_of_word_types=self.number_of_word_types(),
                     number_of_word_tokens=self.number_of_word_tokens(),
                     input_file_path=self.file_abspath)
        vprint('\t' + fname, verbose=verbose)

        fname = 'words_to_sigtransforms.txt'
        obj = double_sorted(self.words_to_sigtransforms().items(),
                            key=lambda x: len(x[1]), reverse=True)
        f_path = os.path.join(output_dir, fname)
        output_latex(obj, f_path,
                     title='Words to sigtransforms',
                     headers=['Word', 'Signature transforms'],
                     row_functions=[lambda x: x[0],
                                    lambda x:
                                    ', '.join(SEP_SIG.join(sig) +
                                              SEP_SIGTRANSFORM + affix
                                              for sig, affix in sorted(x[1]))],
                     column_widths=[20, 0],
                     lxa_parameters=self.parameters(),
                     test=test, encoding=self.encoding,
                     number_of_word_types=self.number_of_word_types(),
                     number_of_word_tokens=self.number_of_word_tokens(),
                     input_file_path=self.file_abspath)
        vprint('\t' + fname, verbose=verbose)

        fname = 'affixes_to_signatures.txt'
        obj = double_sorted(self.affixes_to_signatures().items(),
                            key=lambda x: len(x[1]), reverse=True)
        f_path = os.path.join(output_dir, fname)
        output_latex(obj, f_path,
                     title='Affixes to signatures',
                     headers=['Affix', 'Sig count', 'Signatures'],
                     row_functions=[lambda x: x[0],
                                    lambda x: len(x[1]),
                                    lambda x:
                                    ', '.join(SEP_SIG.join(sig)
                                              for sig in sorted(x[1]))],
                     column_widths=[15, 15, 0],
                     lxa_parameters=self.parameters(),
                     test=test, encoding=self.encoding,
                     number_of_word_types=self.number_of_word_types(),
                     number_of_word_tokens=self.number_of_word_tokens(),
                     input_file_path=self.file_abspath)
        vprint('\t' + fname, verbose=verbose)

        # ----------------------------------------------------------------------
        if self.corpus_file_object:
            vprint('manifold objects', verbose=verbose)

            fname = 'words_to_neighbors.txt'
            obj = list()  # list of tuple(word, list of neighbor words)
            for word in self.wordlist()[: self.parameters()['max_word_types']]:
                obj.append((word, self.words_to_neighbors()[word]))
            f_path = os.path.join(output_dir, fname)
            output_latex(obj, f_path,
                         title='Words to neighbors',
                         headers=['Word', 'Neighbors'],
                         row_functions=[lambda x: x[0],
                                        lambda x: ' '.join(x[1])],
                         column_widths=[25, 0],
                         lxa_parameters=self.parameters(),
                         test=test, encoding=self.encoding,
                         number_of_word_types=self.number_of_word_types(),
                         number_of_word_tokens=self.number_of_word_tokens(),
                         input_file_path=self.file_abspath)
            vprint('\t' + fname, verbose=verbose)

        # ----------------------------------------------------------------------
        vprint('phon objects', verbose=verbose)

        def output_latex_for_phon_words(obj_, f_path_, title_, lxa_parameters_,
                                        test_, encoding_, number_of_word_types_,
                                        number_of_word_tokens_,
                                        input_file_path_):
            output_latex(obj_, f_path_,
                         title=title_,
                         headers=['Word', 'Count', 'Frequency', 'Phones',
                                  'Unigram plog', 'Avg unigram plog',
                                  'Bigram plog', 'Avg bigram plog'],
                         row_functions=[lambda x: x[0],
                                        lambda x: x[1].count,
                                        lambda x:
                                        '%.6f' % x[1].frequency,
                                        lambda x:
                                        ' '.join(x[1].phones),
                                        lambda x:
                                        '%8.3f' % x[1].unigram_plog,
                                        lambda x:
                                        '%8.3f' % x[1].avg_unigram_plog,
                                        lambda x:
                                        '%8.3f' % x[1].bigram_plog,
                                        lambda x:
                                        '%8.3f' % x[1].avg_bigram_plog,
                                        ],
                         column_widths=[35, 10, 15, 60, 15, 15, 15, 15],
                         lxa_parameters=lxa_parameters_,
                         test=test_, encoding=encoding_,
                         number_of_word_types=number_of_word_types_,
                         number_of_word_tokens=number_of_word_tokens_,
                         input_file_path=input_file_path_)

        fname = 'wordlist.txt'
        obj_word_phon = list()  # list of tuple(word, list of neighbor words)
        for word in self._wordlist.mylist:
            obj_word_phon.append((word, self.word_phonology_dict()[word]))
        f_path = os.path.join(output_dir, 'wordlist.txt')
        output_latex_for_phon_words(obj_word_phon, f_path,
                                    'Wordlist sorted by word count',
                                    self.parameters(), test, self.encoding,
                                    self.number_of_word_types(),
                                    self.number_of_word_tokens(),
                                    self.file_abspath)
        vprint('\t' + fname, verbose=verbose)

        fname = 'wordlist_by_avg_unigram_plog.txt'
        obj_unigram_plog = double_sorted(obj_word_phon,
                                         key=lambda x: x[1].avg_unigram_plog,
                                         reverse=False)
        f_path = os.path.join(output_dir, fname)
        output_latex_for_phon_words(obj_unigram_plog, f_path,
                                    'Wordlist sorted by avg unigram plog',
                                    self.parameters(), test, self.encoding,
                                    self.number_of_word_types(),
                                    self.number_of_word_tokens(),
                                    self.file_abspath)
        vprint('\t' + fname, verbose=verbose)

        fname = 'wordlist_by_avg_bigram_plog.txt'
        obj_bigram_plog = double_sorted(obj_word_phon,
                                        key=lambda x: x[1].avg_bigram_plog,
                                        reverse=False)
        f_path = os.path.join(output_dir, fname)
        output_latex_for_phon_words(obj_bigram_plog, f_path,
                                    'Wordlist sorted by avg bigram plog',
                                    self.parameters(), test, self.encoding,
                                    self.number_of_word_types(),
                                    self.number_of_word_tokens(),
                                    self.file_abspath)
        vprint('\t' + fname, verbose=verbose)

        fname = 'phones.txt'
        obj = double_sorted(self.phone_dict().items(),
                            key=lambda x: x[1].count, reverse=True)
        f_path = os.path.join(output_dir, fname)
        output_latex(obj, f_path,
                     title='Phones',
                     headers=['Phone', 'Count', 'Frequency', 'Plog'],
                     row_functions=[lambda x: x[0],
                                    lambda x: x[1].count,
                                    lambda x: '%.6f' % x[1].frequency,
                                    lambda x: '%8.3f' % x[1].plog,
                                    ],
                     column_widths=[10, 10, 15, 15],
                     lxa_parameters=self.parameters(),
                     test=test, encoding=self.encoding,
                     number_of_word_types=self.number_of_word_types(),
                     number_of_word_tokens=self.number_of_word_tokens(),
                     input_file_path=self.file_abspath)
        vprint('\t' + fname, verbose=verbose)

        fname = 'biphones.txt'
        obj = double_sorted(self.biphone_dict().items(),
                            key=lambda x: x[1].count, reverse=True)
        f_path = os.path.join(output_dir, fname)
        output_latex(obj, f_path,
                     title='Biphones',
                     headers=['Biphone', 'Count', 'Frequency',
                              'MI', 'Weighted MI'],
                     row_functions=[lambda x: ' '.join(x[0]),
                                    lambda x: x[1].count,
                                    lambda x:
                                    '%.6f' % x[1].frequency,
                                    lambda x:
                                    '%8.3f' % x[1].MI,
                                    lambda x:
                                    '%8.3f' % x[1].weighted_MI,
                                    ],
                     column_widths=[10, 10, 15, 15, 15],
                     lxa_parameters=self.parameters(),
                     test=test, encoding=self.encoding,
                     number_of_word_types=self.number_of_word_types(),
                     number_of_word_tokens=self.number_of_word_tokens(),
                     input_file_path=self.file_abspath)
        vprint('\t' + fname, verbose=verbose)

        fname = 'triphones.txt'
        obj = double_sorted(self.phone_trigram_counter().items(),
                            key=lambda x: x[1], reverse=True)
        f_path = os.path.join(output_dir, fname)
        output_latex(obj, f_path,
                     title='Triphones',
                     headers=['Triphone', 'Count'],
                     row_functions=[lambda x: ' '.join(x[0]),
                                    lambda x: x[1],
                                    ],
                     column_widths=[15, 10],
                     lxa_parameters=self.parameters(),
                     test=test, encoding=self.encoding,
                     number_of_word_types=self.number_of_word_types(),
                     number_of_word_tokens=self.number_of_word_tokens(),
                     input_file_path=self.file_abspath)
        vprint('\t' + fname, verbose=verbose)

        # ----------------------------------------------------------------------
        vprint('trie objects', verbose=verbose)

        fname = 'words_as_tries.txt'
        obj = list()
        for word in self._wordlist.mylist:
            obj.append((word,
                        self.broken_words_left_to_right()[word],
                        self.broken_words_right_to_left()[word]))
        f_path = os.path.join(output_dir, fname)
        output_latex(obj, f_path,
                     title='Words as tries',
                     headers=['Word', 'Left-to-right trie',
                              'Right-to-left trie'],
                     row_functions=[lambda x: x[0],
                                    lambda x: ' '.join(x[1]),
                                    lambda x: ' '.join(x[2]),
                                    ],
                     column_widths=[35, 50, 50],
                     lxa_parameters=self.parameters(),
                     test=test, encoding=self.encoding,
                     number_of_word_types=self.number_of_word_types(),
                     number_of_word_tokens=self.number_of_word_tokens(),
                     input_file_path=self.file_abspath)
        vprint('\t' + fname, verbose=verbose)

        fname = 'successors.txt'
        obj = double_sorted(self._successors,
                            key=lambda x: len(x[1]), reverse=False)
        f_path = os.path.join(output_dir, fname)
        output_latex(obj, f_path,
                     title='Successors',
                     headers=['String', 'Successors'],
                     row_functions=[lambda x: x[0],
                                    lambda x: ' '.join(sorted(x[1])),
                                    ],
                     column_widths=[35, 0],
                     lxa_parameters=self.parameters(),
                     test=test, encoding=self.encoding,
                     number_of_word_types=self.number_of_word_types(),
                     number_of_word_tokens=self.number_of_word_tokens(),
                     input_file_path=self.file_abspath)
        vprint('\t' + fname, verbose=verbose)

        fname = 'predecessors.txt'
        obj = double_sorted(self.predecessors().items(),
                            key=lambda x: len(x[1]), reverse=False)
        f_path = os.path.join(output_dir, fname)
        output_latex(obj, f_path,
                     title='Predecessors',
                     headers=['String', 'Predecessors'],
                     row_functions=[lambda x: x[0],
                                    lambda x: ' '.join(sorted(x[1])),
                                    ],
                     column_widths=[35, 0],
                     lxa_parameters=self.parameters(),
                     test=test, encoding=self.encoding,
                     number_of_word_types=self.number_of_word_types(),
                     number_of_word_tokens=self.number_of_word_tokens(),
                     input_file_path=self.file_abspath)
        vprint('\t' + fname, verbose=verbose)

    # --------------------------------------------------------------------------
    # for number of word types and tokens

    def number_of_word_types(self):
        """
        Return the number of word types.

        :rtype: int
        """
        if self._number_of_word_types is None:
            self._number_of_word_types = len(self.word_unigram_counter())
        return self._number_of_word_types

    def number_of_word_tokens(self):
        """
        Return the number of word tokens.

        :rtype: int
        """
        if self._number_of_word_tokens is None:
            self._number_of_word_tokens = sum(self.word_unigram_counter()
                                              .values())
        return self._number_of_word_tokens

    # --------------------------------------------------------------------------
    # for the "ngram" module

    def word_unigram_counter(self):
        """
        Return a dict of words with their counts.

        :rtype: dict(str: in)
        """
        if self._word_unigram_counter is None:
            if self.corpus_file_object:
                self._make_word_ngrams_from_corpus_file_object()
            elif self.wordlist_file_object:
                self._read_from_wordlist_file_object()

        return self._word_unigram_counter

    def word_bigram_counter(self):
        """
        Return a dict of word bigrams with their counts.

        :rtype: dict(tuple(str): int)
        """
        if self._word_bigram_counter is None:
            self._make_word_ngrams_from_corpus_file_object()
        return self._word_bigram_counter

    def word_trigram_counter(self):
        """
        Return a dict of word trigrams with their counts.

        :rtype: dict(tuple(str): int)
        """
        if self._word_trigram_counter is None:
            self._make_word_ngrams_from_corpus_file_object()
        return self._word_trigram_counter
#************************
#************************
    """def _make_wordlist(self):
        
        Return a wordlist sorted by word frequency in descending order.
        (So "the" will most likely be the first word for written English.)
        
        word_counter = self.word_unigram_counter()
        word_counter_sorted = double_sorted(word_counter.items(),
                                            key=lambda x: x[1], reverse=True)
        self._wordlist = [word for word, _ in word_counter_sorted]

        self._wordlist = CWordList()
    def wordlist(self):
        
        Return a wordlist sorted by word frequency in descending order.
        (So "the" will most likely be the first word for written English.)

        :rtype: list(str)
        
        if self._wordlist is None:
            self._make_wordlist()
        return self._wordlist
        """
    def _read_from_wordlist_file_object(self):
        word_freq_dict = dict()
        words_to_phones = dict()

        for line in self.wordlist_file_object:
            line = line.strip()
            if not line or line.startswith('#'):
                continue

            word, *rest = line.split()

            if not self.parameters_['keep_case']:
                word = word.lower()

            try:
                freq = int(rest[0])
            except (ValueError, IndexError):
                freq = 1

            phones = rest[1:]
            if not phones:
                phones = list(word)

            word_freq_dict[word] = freq
            words_to_phones[word] = phones

        self._word_unigram_counter = word_freq_dict
        self._words_to_phones = words_to_phones

    def _make_word_ngrams_from_corpus_file_object(self):
        if self.corpus_file_object is None:
            self._word_bigram_counter = dict()
            self._word_trigram_counter = dict()
            return

        unigrams, bigrams, trigrams = ngram.run(
            corpus_file_object=self.corpus_file_object,
            keep_case=self.parameters_['keep_case'],
            max_word_tokens=self.parameters_['max_word_tokens'])

        self._word_unigram_counter = unigrams
        self._word_bigram_counter = bigrams
        self._word_trigram_counter = trigrams

    def run_ngram_module(self, verbose=False):
        """
        Run the ngram module.
        """
        vprint('Extracting word ngrams...', verbose=verbose)
        if self._wordlist is None:
            self._make_wordlist()

    # --------------------------------------------------------------------------
    # for the "signature" module

    def stems_to_words(self):
        """
        Return a dict of stems to words.

        :rtype: dict(str: set(str))
        """
        if self._stems_to_words is None:
            self._make_all_signature_objects()
        return self._stems_to_words

    def signatures_to_stems(self):
        """
        Return a dict of morphological signatures to stems.

        :rtype: dict(tuple(str): set(str))
        """
        if self._signatures_to_stems is None:
            self._make_all_signature_objects()
        return self._signatures_to_stems

    def stems_to_signatures(self):
        """
        Return a dict of stems to morphological signatures.

        :rtype: dict(str: set(tuple(str)))
        """
        if self._stems_to_signatures is None:
            self._make_all_signature_objects()
        return self._stems_to_signatures

    def words_to_signatures(self):
        """
        Return a dict of words to morphological signatures.

        :rtype: dict(str: set(tuple(str)))
        """
        if self._words_to_signatures is None:
            self._make_all_signature_objects()
        return self._words_to_signatures

    def signatures_to_words(self):
        """
        Return a dict of morphological signatures to words.

        :rtype: dict(tuple(str): set(str))
        """
        if self._signatures_to_words is None:
            self._make_all_signature_objects()
        return self._signatures_to_words

    def words_to_sigtransforms(self):
        """
        Return a dict of words to signature transforms.

        :rtype: dict(str: set(tuple(tuple(str), str))
        """
        if self._words_to_sigtransforms is None:
            self._make_all_signature_objects()
        return self._words_to_sigtransforms

    def signatures(self):
        """
        Return a set of morphological signatures.

        :rtype: set(tuple(str))
        """
        if self._signatures is None:
            self._make_all_signature_objects()
        return self._signatures

    def affixes_to_signatures(self):
        """
        Return a dict of affixes to morphological signatures.

        :rtype: dict(str: set(tuple(str)))
        """
        if self._affixes_to_signatures is None:
            self._make_all_signature_objects()
        return self._affixes_to_signatures

    def words_in_signatures(self):
        """
        Return a set of words that are in at least one morphological signature.

        :rtype: set(str)
        """
        if self._words_in_signatures is None:
            self._make_all_signature_objects()
        return self._words_in_signatures

    def affixes(self):
        """
        Return a set of affixes.

        :rtype: set(str)
        """
        if self._affixes is None:
            self._make_all_signature_objects()
        return self._affixes

    def stems(self):
        """
        Return a set of stems.

        :rtype: set(str)
        """
        if self._stems is None:
            self._make_all_signature_objects()
        return self._stems




#put make all signature objects back into this file
    def MakeSignatures(self, FindSuffixesFlag):
            #changed from "Lexicon" to "self" because now within the Lexicon class
            #lxalogfile just prints what the program is doing at any time
            #you can comment out the lxalogfile lines and the outfile_Rebalancing_Signatures lines until ready to use them
        Protostems = dict()
        self.NumberOfAnalyzedWords = 0 #added these to the initialize function
        self.LettersInAnalyzedWords = 0
        self.NumberOfUnanalyzedWords = 0
        self.LettersInUnanalyzedWords = 0

        AnalyzedWords = dict()
        MinimumStemLength = self.MinimumStemLength

        self.TotalRobustnessInSignatures = 0
        self.LettersInStems = 0
        self.TotalLetterCostOfAffixesInSignatures = 0

        #using the wordlist function defined further up
        wordlist_ = self._wordlist.mylist

        #FindProtostems(wordlist, Protostems,minimum_stem_length,FindSuffixesFlag) function rewritten

        minimum_stem_length = MinimumStemLength
        previousword = ""
        if FindSuffixesFlag:
                for i in range(len(wordlist_)): #wordlist changed to the version found using wordlist(self)
                        word = wordlist_[i].Key
                        differencefoundflag = False
                        if previousword == "":  # only on first iteration
                                previousword = word
                                continue
                        span = min(len(word), len(previousword))
                        for i in range(span):
                                if word[i] != previousword[i]: #will a stem be found in the very first word?
                                        differencefoundflag = True
                                        stem = word[:i]
                                        if len(stem) >= minimum_stem_length :
                                                if stem not in Protostems:
                                                        Protostems[stem] = 1
                                                else:
                                                        Protostems[stem] += 1
                                        previousword = word
                                        break
                        if differencefoundflag:
                                continue
                        if len(previousword) > i + 1:
                                previousword = word
                                continue
                        if (len(word)) >= i:
                                if len(previousword) >= minimum_stem_length:
                                        if (previousword not in Protostems):
                                                Protostems[previousword] = 1
                                        else:
                                                Protostems[previousword] += 1
                        previousword = word

        else:
                print ("prefixes")
                ReversedList = list()
                TempList = list()
                for key, val in wordlist_.items():
                        k = key[::-1]
                        TempList.append(k)
                TempList.sort()
                for word in TempList:
                        ReversedList.append(word[::-1])
                for i in range(len(ReversedList)):
                        word = ReversedList[i]
                        differencefoundflag = False
                        if previousword == "":  # only on first iteration
                                previousword = word
                                continue
                        span = min(len(word), len(previousword))
                        for i in range(1,span,):
                                if word[-1*i] != previousword[-1*i]:
                                        differencefoundflag = True
                                        stem = word[-1*i+1:]
                                        if len(stem) >= minimum_stem_length:
                                                if stem not in Protostems:
                                                        Protostems[stem] = 1
                                                else:
                                                        Protostems[stem] += 1
                                        #print previousword, word, stem
                                        previousword = word
                                        break
                        if differencefoundflag:
                                continue
                        if len(previousword) > i + 1:
                                previousword = word
                                continue
                        if (len(word)) >= i:
                                if len(previousword) >= minimum_stem_length:
                                        if (previousword not in Protostems):
                                                Protostems[previousword] = 1
                                        else:
                                                Protostems[previousword] += 1
                        previousword = word
        #end of FindProtostems function
        #print "{:50s}{:>10,}".format("2a. Finished finding proto-stems.", len(Protostems)) 

        #FindAffixes_1(Lexicon, Protostems,  FindSuffixesFlag) function

        MaximumAffixLength = self.MaximumAffixLength  #changes Lexicon.MaximumAffixLength to "self"
        if FindSuffixesFlag:
                for i in range(len(wordlist_)): #also uses new wordlist
                        word = wordlist_[i].Key
                        WordAnalyzedFlag = False
                        for i in range(len(word)-1 , MinimumStemLength-1, -1):
                                stem = word[:i]
                                if stem in Protostems:
                                        suffix = word[i:]
                                        if len(suffix) > MaximumAffixLength:
                                                continue
                                        if stem not in self._stems_to_words:
                                                self._stems_to_words[stem] = dict()
                                        self._stems_to_words[stem][word]=1
                                        if stem not in self._affixes:
                                                self._affixes[stem] = dict()
                                        self._affixes[stem][suffix] = 1
                                        if stem in self._word_unigram_counter:
                                                self._stems_to_words[stem][word] = 1
                                                self._affixes[stem]["NULL"] = 1
                                        self._suffixes[suffix]=1
        else:
                for i in range(len(wordlist_)):
                        word = wordlist_[i].Key
                        WordAnalyzedFlag = False
                        for i in range(MinimumStemLength-1, len(word)-1):
                                stem = word[-1*i:]
                                if stem in Protostems:
                                        j = len(word) - i
                                        prefix = word[:j]
                                        if len(prefix) > MaximumAffixLength:
                                                continue
                                        if stem not in self._stems_to_words:
                                                self._stems_to_words[stem] = dict()
                                        self._stems_to_words[stem][word]=1
                                        if stem not in self._affixes:
                                                self._affixes[stem] = dict()
                                        self._affix[stem][prefix]=1
                                        if stem in self._word_unigram_counter:
                                                self._stems_to_words[stem][word] = 1
                                                self._affixes[stem]["NULL"]=1
                                        self._prefixes[prefix]=1

        #end of FindAffixes_1

        #print "{:50s}{:>10,}".format("2b. Finished finding affixes for protostems.", len(Lexicon._suffixes)+ len(Lexicon._prefixes)) 	

        # It is possible for a stem to have only one affix at this point. We must eliminate those analyses.

        ListOfStemsToRemove = list()
        for stem in self._affixes:
                if len(self._affixes[stem]) < 2:
                        #print "Removing: ", stem, Lexicon.StemToAffix[stem]
                        ListOfStemsToRemove.append(stem)
        #print "{:50s}{:>10,}".format("2c. Number of stems to remove due to mono-fixivity:", len(ListOfStemsToRemove)) 	

        for stem in ListOfStemsToRemove:
                del self._stems_to_words[stem]
                del self._affixes[stem]

        self.LettersInStems =0
        self.TotalLetterCostOfAffixesInSignatures =0
        self.StemsToSignatures(FindSuffixesFlag)

        #print "{:50s}{:>10,}".format("2d. Finished finding signatures.", len(Lexicon._signatures_to_stems))  

	# We look for a stem-final sequence that appears on all or almost all the stems, and shift it to affixes.
	# Make changes in Lexicon.SignatureToStems, and .StemToSig, and .WordToSig, and .StemToWord, and .StemToAffix  and signature_tuples....
	
        threshold = 0.80

	#start of RebalanceSignatureBreaks2 (Lexicon, threshold, outfile_Rebalancing_Signatures, FindSuffixesFlag)

        count=0
        MinimumNumberOfStemsInSignaturesCheckedForRebalancing = 5
        SortedListOfSignatures = sorted(self._signatures_to_stems.items())#, lambda x, y: cmp(len(x[1]), len(y[1])),
									#reverse=True)
        for (sig,wordlist_) in SortedListOfSignatures:
                sigstring="-".join(sig)
                numberofstems=len(Self._signatures_to_stems[sig])

                if numberofstems <MinimumNumberOfStemsInSignaturesCheckedForRebalancing:
                        #print >>outfile_Rebalancing_Signatures, "       Too few stems to shift material from suffixes", sigstring, numberofstems
                        continue
                #print >>outfile_Rebalancing_Signatures, "{:20s} count: {:4d} ".format(sigstring,   numberofstems),
                shiftingchunk, shiftingchunkcount  = TestForCommonEdge(Self._signatures_to_stems[sig], outfile, threshold, FindSuffixesFlag)
                #if shiftingchunkcount > 0:
                        #print "CC" ,sig, shiftingchunk
                        #print >>outfile,"{:5s} count: {:5d}".format(shiftingchunk,   shiftingchunkcount)
                #else:
                        #print >>outfile_Rebalancing_Signatures, "no chunk to shift"
                if len(shiftingchunk) >0:
                        count +=1
                        chunklength = len(shiftingchunk)
                        newsignature = list()
                        sigstring2 = ""
                        for affix in sig:
                                if affix == "NULL":
                                        affix = ""
                                if FindSuffixesFlag:
                                        newaffix = shiftingchunk + affix
                                        sigstring2 += newaffix + " - "
                                        shiftpair = (affix, newaffix)
                                        sigstring2 = sigstring2[:-2] ##????
                                else:
                                        newaffix = affix + shiftingchunk
                                        sigstring2 += "-" + newaffix
                                        shiftpair = (affix, newaffix)
                                        sigstring2 = sigstring2[2:] ##???
                                newsignature.append(newaffix)
                        #formatstring = "{:30s} {:10s} {:35s}  Number of stems {:5d} Number of shifters {:5d}"
                        #print >>outfile_Rebalancing_Signatures, formatstring.format(sigstring, shiftingchunk, sigstring2, numberofstems, shiftingchunkcount)
                        if shiftingchunkcount >= numberofstems * threshold:
                                ChangeFlag = True
                                stems_to_change = list(self._signatures_to_stems[sig])
                                for stem in stems_to_change:
                                        if FindSuffixesFlag:
                                                if stem[-1*chunklength:] != shiftingchunk:
                                                        continue
                                        else:
                                                if stem[:chunklength] != shiftingchunk:
                                                        continue
                                        if FindSuffixesFlag:
                                                newstem = stem[:len(stem)-chunklength]
                                        else:
                                                newstem = stem[chunklength:]
                                        if newstem not in self._stems_to_words:
                                                self._stems_to_words[newstem] = dict()
                                        for word in self._stems_to_words[stem]:
                                                self._stems_to_words[newstem][word] = 1
                                        del self._stems_to_words[stem] #  is this too general?
                                        if newstem not in self._affixes:
                                                self._affixes[newstem] = {}
                                        for affix in newsignature:
                                                self._affixes[newstem][affix] = 1
                                        del self._affixes[stem]
                                        if newstem not in self._stems_to_signatures:
                                                self._stems_to_signatures[newstem]=dict()
                                        self._stems_to_signatures[newstem]=[newsignature]
                                        del self._stems_to_signatures[stem]
	
        #outfile_Rebalancing_Signatures.flush()

        #end of RebalanceSignatureBreaks2
	
        #print "{:50s}{:10d}".format("2e. Rebalance signature breaks.",count)
        if True:
                if FindSuffixesFlag:
                        Affixes = self._suffixes
                else:
                        Affixes = self._prefixes
        #FindSignatureStructure (Lexicon,FindSuffixesFlag, lxalogfile, Affixes, affix_threshold=3) function

		
        # This function assumes that we have the set of stems already in Lexicon.SignatureToStems. It does the rest.

        affix_threshold = 3
        StemList = self._stems_to_words.keys()
        self.StemToSig = {}
        self.StemToAffix = {}# ok
        self.StemToWord = dict()# ok
        self.Signatures = {}
        self.SignatureToStems = {}
        self.WordToSig = {}
        self.StemCounts = {}
 


	#  Signatures with way too many affixes are spurious.
	# If we have already run this function before, we have a set of affixes ranked by frequency,
	# and we can use these now to eliminate low frequency suffixes from signatures with
	# just one affix. 
	# Lexicon.MaximumNumberOfAffixesInASignature
	# Number of affixes we are confident in:

        #print ("2f1. Finding affixes we are confident about.")
        number_of_affixes_per_line = 10



        #if len (Affixes ) ==0:
                #if FindSuffixesFlag:
                        #print ("No suffixes found yet.")
                #else:
                        #print ("No prefixes found yet.")
        if len(Affixes) != 0:
                NumberOfConfidentAffixes = 50
                ConfidentAffixes = dict()
                SortedAffixes = list(Affixes.keys())
                SortedAffixes.sort(key=lambda affix:Affixes[affix], reverse=True)

                #print ("Confidence affixes:")
                count = 0
                for affixno in range(NumberOfConfidentAffixes):
                        ConfidentAffixes[SortedAffixes[affixno]]=1
                        #print  "{:6s} ".format(SortedAffixes[affixno]),
                        count += 1
                        if count % number_of_affixes_per_line == 0:
                                print
                for sig in self._signatures_to_stems:
                        stems = self._signatures_to_stems[sig]
                        newsig = list()
                        if len(stems) == 1:
                                for affix in sig:
                                        if affix in ConfidentAffixes:
                                                newsig.append(affix)
                print

	# Creates Lexicon.StemToWord, Lexicon.StemToAffix.
        #print ("Reanalyzing words.")
        if True:
                #print ("Word number:")
                for i in range(len(self._wordlist.mylist)):
                        #if i % 2500 == 0:
                                #print "{:7,d}".format(i),
                                #sys.stdout.flush()
                        word = self._wordlist[i].Key

                        WordAnalyzedFlag = False
                        for i in range(len(word)-1 , self.MinimumStemLength-1, -1):
                                if FindSuffixesFlag:
                                        stem = word[:i]
                                else:
                                        stem = word[-1*i:]
                                if stem in StemList:
                                        affixlength = len(word)-i
                                        if FindSuffixesFlag:
                                                affix = word[i:]
                                        else:
                                                affix = word[:affixlength]
                                        if len(affix) > self.MaximumAffixLength:
                                                continue
                                        # the next line involves putting a threshold on frequencies of affixes.
                                        if Affixes and affix in Affixes and Affixes[affix] < affix_threshold:
                                                continue
                                        #print stem, suffix
                                        if stem not in self._stems_to_words:
                                                self._stems_to_words[stem] = dict()
                                        self._stems_to_words[stem][word]=1
                                        if stem not in self._affixes:
                                                self._affixes[stem] = dict()
                                        self._affixes[stem][affix] = 1
                                        if stem in self._word_unigram_counter:
                                                self._stems_to_words[stem][word] = 1
                                                self._affixes[stem]["NULL"] = 1
                                        if stem not in self.StemCounts:
                                                self.StemCounts[stem] = 0
                                        self.StemCounts[stem]+= self._word_unigram_counter[word]

        print
        self.LettersInStems =0
        self.TotalLetterCostOfAffixesInSignatures =0

        if (False):
                for stem in self._affixes:
                        if len(_affixes[stem]) > self.MaximumNumberOfAffixesInASignature:
                                for sig in self._signatures_to_stems:
                                        stems = self._signatures_to_stems[sig]
                                        newsig = list()
                                        if len(stems) == 1:
                                                for affix in sig:
                                                        if affix in ConfidentAffixes:
                                                                newsig.append(affix)

        print ("Finding a first set of signatures.")
			 
        # From StemToAffix and StemToWord:
        # Create SignaturesToStems, Suffixes, StemToSig, WordToSig
        StemsToEliminate = list()
        for stem in self._stems_to_words:
                self.LettersInStems += len(stem)
                signature = list(self._affixes[stem])
                #print stem, signature

                signature.sort()
                signature_tuple = tuple(signature)

                if len(signature) == 1:
                        StemsToEliminate.append(stem)
                        continue

                if signature_tuple not in self._signatures_in_stems:
                        self._signatures_in_stems[signature_tuple] = dict()
                        for affix in signature:
                                self.TotalLetterCostOfAffixesInSignatures += len(affix)
                                if affix not in Affixes:
                                        Affixes[affix]=1
                                else:
                                        Affixes[affix] +=1
                self._signatures_to_stems[signature_tuple][stem] = 1

                self._stems_to_signatures[stem] = signature_tuple
                for word in self._stems_to_words[stem]:
                        if word not in self._words_to_signatures:
                                self._words_to_signatures[word] = list()
                        self._words_to_signatures[word].append(signature_tuple)
                        self.LettersInAnalyzedWords += len(word)

        for stem in StemsToEliminate:
                del self._affixes[stem]
                del self._stems_to_words[stem]

        #print "{:50s}{:>10,}".format("2h. Signatures.", len(self._signatures_to_stems))

        print  ("Finished redoing structure")
        #print >>lxalogfile, "Finished redoing structure \n "
        #print "Number of sigs: ", len(Lexicon.SignatureToStems)
        #print >>lxalogfile, "Number of sigs: ", len(self._signatures_to_stems)

        self.ComputeRobustness()

        #print >> lxalogfile, "{:40s}{:10,d}".format("Number of analyzed words", self.NumberOfAnalyzedWords)
        #print >> lxalogfile, "{:40s}{:10,d}".format("Number of unanalyzed words", self.NumberOfUnanalyzedWords)
        #print >> lxalogfile, "{:40s}{:10,d}".format("Letters in stems", self.LettersInStems)
        #print >> lxalogfile, "{:40s}{:10,d}".format("Letters in affixes", self.AffixLettersInSignatures)
        #print >> lxalogfile, "{:40s}{:10,d}".format("Total robustness in signatures", self.TotalRobustnessInSignatures)

        return


    def _make_all_signature_objects(self):
       # self._stems_to_words = signature.make_stems_to_words(
       #     self.wordlist(), self.parameters_['min_stem_length'],
       #     self.parameters_['max_affix_length'], self.parameters_['suffixing'],
       #     self.parameters_['min_sig_count'])

        self.MakeSignatures(1)

        self._signatures_to_stems = signature.make_signatures_to_stems(
            self._stems_to_words, self.parameters_['max_affix_length'],
            self.parameters_['min_sig_count'], self.parameters_['suffixing'])

        self._stems_to_signatures = signature.make_stems_to_signatures(
            self._signatures_to_stems)

        self._words_to_signatures = signature.make_words_to_signatures(
            self._stems_to_words, self._stems_to_signatures)

        self._signatures_to_words = signature.make_signatures_to_words(
            self._words_to_signatures)

        self._words_to_sigtransforms = signature.make_words_to_sigtransforms(
            self._words_to_signatures, self.parameters_['suffixing'])

        self._signatures = set(self._signatures_to_stems.keys())

        self._affixes_to_signatures = signature.make_affixes_to_signatures(
            self._signatures)

        self._words_in_signatures = set(self._words_to_signatures.keys())
        self._affixes = set(self._affixes_to_signatures.keys())
        self._stems = set(self._stems_to_words.keys())

    def run_signature_module(self, verbose=False):
        """
        Run the signature module.
        """
        vprint('Morphological signatures...', verbose=verbose)
        self._make_all_signature_objects()

    # --------------------------------------------------------------------------
    # for the "manifold" module

    def words_to_neighbors(self):
        """
        Return a dict of words to syntactic neighbors.

        :rtype: dict(word: list(str))
        """
        if self._words_to_neighbors is None:
            self._make_all_manifold_objects()
        return self._words_to_neighbors

    def words_to_contexts(self):
        """
        Return a dict of words to contexts with counts.

        :rtype: dict(str: dict(tuple(str): int))
        """
        if self._words_to_contexts is None:
            self._make_all_manifold_objects()
        return self._words_to_contexts

    def contexts_to_words(self):
        """
        Return a dict of contexts to words with counts.

        :rtype: dict(tuple(str): dict(str: int))
        """
        if self._contexts_to_words is None:
            self._make_all_manifold_objects()
        return self._contexts_to_words

    def neighbor_graph(self):
        """
        Return the syntactic word neighborhood graph.

        :rtype: networkx undirected graph
        """
        if self._neighbor_graph is None:
            self._make_all_manifold_objects()
        return self._neighbor_graph

    def _make_all_manifold_objects(self):
        self._words_to_neighbors, self._words_to_contexts, \
        self._contexts_to_words = manifold.run(
            self.word_unigram_counter(),
            self.word_bigram_counter(),
            self.word_trigram_counter(),
            self.parameters_['max_word_types'],
            self.parameters_['n_neighbors'],
            self.parameters_['n_eigenvectors'],
            self.parameters_['min_context_count'])
        self._neighbor_graph = manifold.compute_graph(self._words_to_neighbors)

    def run_manifold_module(self, verbose=False):
        """
        Run the phon module.
        """
        vprint('Syntactic word neighbors...', verbose=verbose)
        if self.corpus_file_object:
            self._make_all_manifold_objects()

    # --------------------------------------------------------------------------
    # for the "phon" module

    def phone_unigram_counter(self):
        """
        Return a dict of phone unigrams with counts.

        :rtype: dict(str: int)
        """
        if self._phone_unigram_counter is None:
            self._make_all_phon_objects()
        return self._phone_unigram_counter

    def phone_bigram_counter(self):
        """
        Return a dict of phone bigrams with counts.

        :rtype: dict(tuple(str): int)
        """
        if self._phone_bigram_counter is None:
            self._make_all_phon_objects()
        return self._phone_bigram_counter

    def phone_trigram_counter(self):
        """
        Return a dict of phone trigrams with counts.

        :rtype: dict(tuple(str): int)
        """
        if self._phone_trigram_counter is None:
            self._make_all_phon_objects()
        return self._phone_trigram_counter

    def phone_dict(self):
        """
        Return a dict of phone unigrams to Phone objects.
        A Phone instance has the methods
        ``spelling()``, ``count()``, ``frequency()``, and ``plog()``.

        :rtype: dict(str: Phone instance)
        """
        if self._phone_dict is None:
            self._make_all_phon_objects()
        return self._phone_dict

    def biphone_dict(self):
        """
        Return a dict of phone bigrams to Biphone objects.
        A Biphone instance has the methods
        ``spelling()``, ``count()``, ``frequency()``, ``MI()``, and
        ``weighted_MI()``.

        :rtype: dict((str, str): Biphone instance)
        """
        if self._phone_dict is None:
            self._make_all_phon_objects()
        return self._biphone_dict

    def word_phonology_dict(self):
        """
        Return a dict of words to Word objects.
        A Word instance has the methods
        ``spelling()``, ``phones()``, ``count()``, ``frequency()``,
        ``unigram_plog()``, ``avg_unigram_plog()``,
        ``bigram_plog()``, and ``avg_bigram_plog()``.

        :rtype: dict(str: Word instance)
        """
        if self._word_dict is None:
            self._make_all_phon_objects()
        return self._word_dict

    def words_to_phones(self):
        """
        Return a dict of words with their phones.

        :rtype: dict(str: list(str))
        """
        return self._words_to_phones

    def _make_all_phon_objects(self):
        word_unigram_counter = self.word_unigram_counter()
        words_to_phones = self.words_to_phones()

        self._phone_unigram_counter, self._phone_bigram_counter, \
        self._phone_trigram_counter = phon.make_word_ngrams(
            word_unigram_counter, words_to_phones)

        self._phone_dict = phon.make_phone_dict(self._phone_unigram_counter)
        self._biphone_dict = phon.make_biphone_dict(self._phone_bigram_counter,
                                                    self._phone_dict)
        self._word_dict = phon.make_word_dict(self.word_unigram_counter(),
                                              self.phone_dict(),
                                              self.biphone_dict(),
                                              self.words_to_phones())

    def run_phon_module(self, verbose=False):
        """
        Run the phon module.
        """
        vprint('Phonology...', verbose=verbose)
        self._make_all_phon_objects()

    # --------------------------------------------------------------------------
    # for the "trie" module

    def broken_words_left_to_right(self):
        """
        Return a dict of words to their left-to-right broken form.

        :rtype: dict(str: list(str))
        """
        #if self._broken_words_left_to_right is None:
            #self._make_all_trie_objects()
        return self._broken_words_left_to_right

    def broken_words_right_to_left(self):
        """
        Return a dict of words to their right-to-left broken form.

        :rtype: dict(str: list(str))
        """
        #if self._broken_words_right_to_left is None:
            #self._make_all_trie_objects()
        return self._broken_words_right_to_left

    def successors(self):
        """
        Return a dict of word (sub)strings to their successors.

        :rtype: dict(str: set(str))
        """
        #if self._successors is None:
            #self._make_all_trie_objects()
        return self._successors

    def predecessors(self):
        """
        Return a dict of word (sub)strings to their predecessors.

        :rtype: dict(str: set(str))
        """
        #if self._predecessors is None:
            #self._make_all_trie_objects()
        return self._predecessors

    """
    def _make_all_trie_objects(self):
        self._broken_words_left_to_right, self._broken_words_right_to_left, \
        self._successors, self._predecessors = trie.run(
            self._wordlist().mylist, self.parameters_['min_stem_length'])
    """
    def run_trie_module(self, verbose=False):
        """
        Run the trie module.
        """
        vprint('Tries...', verbose=verbose)
        #self._make_all_trie_objects()

    def print_wordlist(self, filename): #for a dx1 file

        infile = open(filename)
        filelines = infile.readlines()

        for line in filelines:
	#print line
                pieces = line.split()
                word=pieces[0]
                if (word == "#"):
                        continue
                if word.isalpha() == False:
                        continue
                if len(pieces) > 1:
                        count = int(pieces[1])
                else:
                        count =1
                word = word.lower()
                if (BreakAtHyphensFlag and '-' in word):
                        words = word.split('-')
                        for word in words:
                                if word.isalpha() == False:
                                        continue
                                if word not in self._word_unigram_counter:
                                        self._word_unigram_counter[word]=0
                                self._word_unigram_counter[word]+=count
                                self.TotalLetterCountInWords  += len(word)
                else:
                        if word not in self._word_unigram_counter:
                                self._word_unigram_counterword=0
                        self._word_unigram_counter[word]  = count
                        Lexicon.TotalLetterCountInWords += len(word)
                if len(self._word_unigram_counter) >= numberofwords:
                        break
        self._wordlist.append(Word(word))
        self.sort_wordlist()
        print ("\n1. Finished reading word list.\n")
        self.PrintWordList(outfile_WordList) #what file does it get printed to?

