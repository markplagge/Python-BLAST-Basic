import itertools as itr
import numpy as np
import argparse
import os
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from scipy.spatial.distance import hamming
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2


class hmatch():
	def __init__(self, qs, ds, score, qstart, qend, dstart, dend):
		self.qs = ""
		self.ds = ""
		for i in qs:
			self.qs += i
		for i in ds:
			self.ds += i
		self.score = score
		self.qstart = qstart
		self.qend = qend
		self.dstart = dstart
		self.dend = dend

	def __str__(self):

		outs = "score= " + str(self.score) + " Q:s, e= " + \
			   str(self.qstart) + " " + str(self.qend) + \
			   " D: s, e=" + str(self.dstart) + " " + str(self.dend) + "\n"
		match = ""
		for c, v in zip(list(self.qs), list(self.ds)):
			if c == v:
				match += "|"
			else:
				match += " "
		outs += self.qs + "\n"
		outs += match + "\n"
		outs += self.ds + "\n"

		return outs

	def __gt__(self, other):
		return self.score > other.score

	def __hash__(self, **kwargs):
		hs = self.qs + self.ds + str(self.score)
		return hash(hs)

	def __eq__(self, other):
		tv = str(self) == str(other)
		tv = tv or self.__hash__() == other.__hash__()
		return tv

	def __lt__(self, other):

		return self.score < other.score


class resultTable():
	def __init__(self, queryS, database, match, mismatch, gap, l, k, t, v, m, e, g,
				 validStrs="CGAT", searchForStrings=True):
		self.q = queryS
		self.d = database

		self.match = match
		self.mismatch = mismatch
		self.gap = gap
		self.l = l
		self.K = k
		self.T = t
		self.V = v
		self.M = m
		self.E = e
		self.G = g
		self.alphabet = validStrs
		if searchForStrings:
			self.downloadData()
		else:
			#given a string
			self.q = Seq(queryS)
			self.d = Seq(database)
		self.t = np.floor((self.T - self.l*self.match)/(self.miss - self.match))

		self.words = {}
		#Generate neighbors
		#self.alphabet = []

		#for p in itr.product(validStrs, repeat=l):
		#	self.alphabet.append(p)
	def __addWord(self,word,type,position):
		if word not in self.words.keys():
			self.words[word] = []
		self.words[word].append(type,position)

	def getW(self,word):
		try:
			return self.words[word]
		except:
			return False
	def calcWordScore(self,op,new):
		#t = np.floor((self.T - l * self.match) / (self.mismatch * self.match))
		score = 0
		for i in range(0, len(new)):
			score += self.S(op[i],new[i])
		return score

	def S(self,i,j):
		if i == j:
			return self.match
		return self.miss



	def findNeighbors(self, word, position):
		for i in itr.combinations(range(len(word)), self.l):
			cc = [[char] for char in word]
			for i in cc:
				orig_char = word[i]
				cc[i] = [l for l in self.alphabet]
			for poss in itr.product(*cc):
				yield ''.join(poss)

	def addWordWithNeighbors(self, word, position):
		pass

	def speedHamIT(self,testWord):
		ct = 0
		for i,j in zip(testWord,self.cword):
			if i != j:
				ct +=1
			if ct > self.t:
				return False
		return True



	def walkQuery(self):
		for i in range(0,len(self.q) - self.l):
			word = list(self.q[i:(i+self.l)])
			it = itr.product(self.alphabet, repeate=self.l)
			self.__addWord(word,"Q",i)
			self.cword = word

			rit = filter(self.speedHamIT, it)
			for res in rit:
				self.__addWord("".join(res),"Q",i)













	def downloadData(self):
		self.d = self.downloadData(self.d)
		self.q = self.downloadData(self.q)



	def downloadDataS(dataID):
		fn = str(dataID) + ".gb"
		Entrez.email = "plaggm@rpi.edu"
		if not os.path.isfile(fn):
			kc = Entrez.efetch("nuccore", id=dataID, rettype="gb",
							   retmode="text")
			kp = SeqIO.parse(kc, "gb")
			with open(fn, "w") as f:
				SeqIO.write(kp, f, "gb")
			kc = kp

		# with open("hum.gb","r") as f:
		# kc = SeqIO.parse(f,"gb")
		kc = SeqIO.read(fn, "gb")
		return kc.seq


query = "AGCTTTTTGGTTAATTCCCACTT"
db = "CATAAGCTTCTGACTCTTACCTCCCTCTCTCCTACTCCTGCTCGCATCTGCTATAGTGGAGGC"

match = 1
mismatch = -1
gap = -1
l = 5
K = 10
T = 2
V = 5
M = 4
E = 2
G = 20


def downloadData(dataID):
