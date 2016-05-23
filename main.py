import sys
import os.path
import operator
import math
import argparse
import re
 
#default messages used in the system
message = {
	"corpus_not_found": "Corpus file does not exist",
	"filter_not_found": "Filter file does not exist",
	"filters_default": "NOTE: No filters will be used in the program"
}

#configuration of parser to handle command line parameters
parser = argparse.ArgumentParser(description='Corpus processing and collocation extraction')
parser.add_argument('corpus', help='relative or absolute path to the corpus file', type=str, metavar='/path/to/corpus')
parser.add_argument('--output', help='flag to indicate what is the required output, e.g.: --output xi pmi. By default it shows all the computed collocation extranctions', nargs='+', default=["2-gram", "pmi", "xi", "dlog"], type=str, choices=["2-gram", "pmi", "xi", "dlog"])
parser.add_argument('--filters', help='flag to indicate relative or absolute paths to skip gram vectors files, e.g.: --filters /path/to/punctuation.txt /path/to/stoplist.txt', nargs='+', type=str)

args = parser.parse_args()

#return the highest n values from a dictionary obj, in where each key has as value a float number
def the_most (obj, count):
	return sorted(obj.items(), key=operator.itemgetter(1))[len(obj) - count:][::-1]

def log (number, base = 10):
	#fallback in case of the number is 0
	if number == 0:
		number = .1e-300
	return math.log(number, base)

#turn a text into lower case, remove non-necesary spaces and apply the skip gram vectors filters if they are refered
def clean_text (line, skip_gram_vectors):
	#converts all tokens to lowercase
	words = line.lower().strip()

	if len(skip_gram_vectors):
		words = words.split(" ")
		for s in skip_gram_vectors:
			#add to words list just those that are different than space and different to any stopword or punctuation mark
			words = [w for w in words if w and w != s]
		
		return words
	else:
		#remove more than one space between tokens
		return " ".join(words.split()).split(" ")

#handling relative or absolute paths
if args.corpus[0] == ".":
	args.corpus = os.path.dirname(os.path.realpath(__file__)) + "/" + args.corpus

#exeption if the corpus file does not exist
if not os.path.isfile(args.corpus):
	print message["corpus_not_found"]
	sys.exit(0)

skip_gram_vectors = []
if args.filters:
	for _file in args.filters:
		#handling relative or absolute paths for filter files
		if _file == ".":
			_file = os.path.dirname(os.path.realpath(__file__)) + "/" + _file

		#exeption if the filter file does not exist
		if not os.path.isfile(_file):
			print message["filter_not_found"]
			sys.exit(0)

		with open(_file) as f:
			lines = f.readlines()
			for i in lines:
				#removing the EOL element of each element of the file
				skip_gram_vectors.append(i.strip())
		f.closed
else:
	print message["filters_default"]

with open(args.corpus) as f:
	print "Processing the corpus..."
	#"clean" every line according to skip_gram_vectors
	content = clean_text(f.read().replace("\n", ""), skip_gram_vectors)
	#2-gram collocations computing
	_2_grams = {}
	_map = {}
	
	N = len(content)
	print "Words in corpus: " + str(N)
	for i in range(N):
		word = content[i]
		
		if i < N - 1:
			#computes all 2-gram collocations of the corpus and their frequencies
			current_2_gram = word + " " + content[i + 1]

			if current_2_gram not in _2_grams:
				_2_grams[current_2_gram]  = 1
			else:
				_2_grams[current_2_gram] += 1

		#adding individual words of the corpus to _map and computing their frequencies

		if word not in _map:
			_map[word]  = 1
		else:
			_map[word] += 1
	
	print " RESULTS\n========="

	if "2-gram" in args.output:
		#outputs 20 most frequent collocations and their corresponding frequencies
		most_frequent = the_most(_2_grams, 20)

		print "Most Frequent 2-grams Collocations:\n-----------------------------------"
		print most_frequent

	if "pmi" in args.output:
		pointwise = {}

	if "xi" in args.output:
		chi_square = {}

	if "dlog" in args.output:
		dunning = {}

		#defining function L
		def L(k, n, x):
			return pow(x, k) * pow(1 - x, n - k)

	for i in _2_grams:
		parts = i.split(" ")
		c12 = float(_2_grams[i])
		c1  = float(_map[parts[0]])
		c2  = float(_map[parts[1]])

		#computing collocations as per point-wise mutual information
		if "pmi" in args.output:
			P_xy = c12 / N
			P_x  = c1 / N
			P_y  = c2 / N
			pointwise[i] = log(P_xy / (P_x * P_y), 2)

		#computing collocations as per Pearson's Chi-squared test
		if "xi" in args.output:
			O11 = c12
			O12 = c1 - c12
			O21 = c2 - c12
			O22 = N - O11 - O12 - O21
			chi_square[i] = N * pow(O11*O22 - O12*O21, 2) / ((O11 + O12)*(O11 + O21)*(O12 + O22)*(O21 + O22))

		#computing collocations as per Dunning's log-likelihood test
		if "dlog" in args.output:
			p  = c2 / N
			p1 = c12 / c1
			p2 = (c2 - c12) / (N - c1)
			
			#Calculating Dunning's log-likelihood and then multiplying each result by -2 to see the statistical significance
			dunning[i] = -2 * (
				log( L(c12, c1, p) ) + 
				log( L(c2 - c12, N - c1, p) ) - 
				log( L(c12, c1, p1) ) - 
				log( L(c2 - c12, N - c1, p2) )
			)

	if "pmi" in args.output:
		highest_scoring_pmi = the_most(pointwise, 20)
		print "\nHighest Scoring Point-wise Mutual Information:\n----------------------------------------------"
		print highest_scoring_pmi

	if "xi" in args.output:
		highest_scoring_chi_square = the_most(chi_square, 20)
		print "\nHighest Scoring Pearson's Chi-squared test:\n-------------------------------------------"
		print highest_scoring_chi_square

	if "dlog" in args.output:
		highest_scoring_dunning = the_most(dunning, 20)
		print "\nHighest Scoring Dunning's log-likelihood test:\n----------------------------------------------"
		print highest_scoring_dunning

f.closed