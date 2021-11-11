import numpy as np
import copy 

class genome:
	n=0
	__array = np.zeros((n, 4))
	def __init__(self,length):
		self.array = np.zeros((length,4))
		self.n = length


	def setg(self,i,value):
		self.array[i][0]=value	


	def setc(self,i,value):
		self.array[i][1]=value

	def setd(self,i,value):
		self.array[i][2]=value

	def setx(self,i,value):
		self.array[i][3]=value

	def getg(self,i):
		return self.array[i][0]

	def getc(self,i):
		return self.array[i][1]

	def getd(self,i):
		return self.array[i][2]

	def getx(self,i):
		return self.array[i][3]

	def getn(self):
		return self.n

	def get_row(self,i):
		return self.array[i]

	def print_genome(self):
		print(self.array)
 

	def set_array(self,temp):
		self.array = copy.deepcopy(temp)
		self.n = len(self.array)

	def get_array(self):
		return self.array 