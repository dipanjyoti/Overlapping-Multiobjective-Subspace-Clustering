from genome import genome
import numpy as np
import random
from pandas import read_table
from sklearn import preprocessing as pp
import math
from datetime import datetime
import time
import copy 
import pygmo as pg
import pandas as pd
import itertools

NOT_ZERO_DIVISION_SECURITY = 1e-10
start_time = time.clock()

filename= 'breast.arff' #N=10,c=10,T=1000,n=60
filenametrue= 'breast.true' #N=10,c=10,T=1000,n=60

genome_size=100 #Numbewr of tuples
sample_size = 150
ss=sample_size

N=10 # Number of Solutions
c=34 #Maximum Clusters
T=2000 #Number of iterations

threshold_fuzzy_Value = 0.1


def No_of_Clusters(M,member):
	Model=M
	membershipList=member
	clusters=0
	for m in range(len(Model)):
		if np.count_nonzero(Model[m]!=0) and np.count_nonzero(membershipList[m]!=0):
			clusters=clusters+1
	return clusters

def FeatureSet(Model,m):
	featureSet=[]
	for i in range (dimension):
		if Model[m][i]!=0:
			featureSet.append(i)
	return set(featureSet)

def dimension_non_redundancy(M,member): 
	Model=M
	membershipList=member
	featureSet=0
	Total_Cluster= No_of_Clusters(Model,membershipList)
	for i in range(len(Model)-1):
		if np.count_nonzero(Model[i]!=0) and np.count_nonzero(membershipList[i]!=0):
			featureSet_i=FeatureSet(Model,i)
			for j in range(i+1, len(Model)):
				if np.count_nonzero(Model[j]!=0) and np.count_nonzero(membershipList[j]!=0):
					featureSet_j=FeatureSet(Model,j)
					#featureSet+=(len(featureSet_i.intersection(featureSet_j)))/(dimension+0.0)  # Modified # confuse/ doubt
					featureSet+=len(featureSet_i.intersection(featureSet_j))
	return ((featureSet * 2)/((0.0+NOT_ZERO_DIVISION_SECURITY)+Total_Cluster*(Total_Cluster-1)))

def Feature_Per_Cluster(Model, membershipList):
	elements=0
	Total_Cluster=No_of_Clusters(Model,membershipList)
	for i in range(len(Model)):
		elements=elements+np.count_nonzero(Model[i])
	FPC=(elements+0.0)/Total_Cluster
	FPC=abs((dimension/2)-FPC)
	return FPC  #/(dimension+0.0)


def Calculate_PSM(M,member):
	Model=M
	membershipList=member
	Model=np.asarray(Model)
	membershipList=np.asarray(membershipList)

	DNR=dimension_non_redundancy(Model,membershipList)
	#DC=dimension_coverage(Model,membershipList)
	FPC=Feature_Per_Cluster(Model, membershipList)
	PSM=DNR+FPC
	return PSM

def membership_non_redundency(mem_array):
    temp = np.sum(mem_array,axis=1)
    n = len(mem_array) - len([i for i in range(len(temp)) if temp[i] == 0])
    temp1 = np.sum(mem_array,axis=0)
    m_n_r=sum(temp1)-ss
    
    # m_n_r = 0
    # for i in range(len(mem_array)):
    #     for j in range(i+1,len(mem_array)):
    #         temp1=mem_array[i]
    #         temp2=mem_array[j]
    #         new=temp1+temp2
    #         for k in range(len(new)):
    #             if(new[k]==2):
    #                 m_n_r=m_n_r+1
    # print m_n_r
    if(n==1):
        return m_n_r
    else:
        return (m_n_r*2)/((n)*(n-1)+0.0)
            
            
#feature_non_redundency
def feature_non_redundency(c):
    temp_array = np.delete(c, 0, 1)

    C = len(temp_array)
    
    if(C==1):             ## error detection
        print temp_array
        print 'halt'
    
    D = d
    for i in range(C):
        for j in range(D):
            if(temp_array[i][j]== None):
                temp_array[i][j] = 0
            else:
                temp_array[i][j] = 1
    #print temp_array

    n = len(temp_array) 
    f_n_r = 0
    for i in range(len(temp_array)):
        for j in range(i+1,len(temp_array)):
            temp1=temp_array[i]
            #print 'a',temp1
            temp2=temp_array[j]
            #print 'b',temp2
            new=temp1+temp2
            #print 'c',new
            for k in range(len(new)):
                if(new[k]==2):
                    f_n_r=f_n_r+1
    if(C==1):
        d_n_r=f_n_r
    else:
        d_n_r=(f_n_r*2)/((C)*(C-1)+0.0)
    dim_cov_array= np.sum(temp_array,axis=0)
    #print 'a',dim_cov_array
    dim_cov = len(dim_cov_array) - len([i for i in range(len(dim_cov_array)) if dim_cov_array[i] == 0])
    #print 'b',dim_cov
    return (d_n_r *D)/(dim_cov+0.0)



#function to calculate distance between a data and its corrosponding cluster center
def Calculate_compactness(M,S, member):
	#print member
	membershipList=member
	Model=M
	Model=np.asarray(Model)
	membershipList=np.asarray(membershipList)
	sampleData=S
	sample_size=len(sampleData)
	SDmax=len(member)
	compact=0.0
	flag=0
	NoCluster=0
	for m in range(SDmax):
		E_c_dis=0.0
		if np.count_nonzero(Model[m]!=0) and np.count_nonzero(membershipList[m]!=0):
			points=np.count_nonzero(membershipList[m])
			NoCluster=NoCluster+1
			flag=1
			for s in range(len(sampleData)):
				if (membershipList[m][s]!=0):
					DisCal=0
					for d in range(dimension):
						if float(Model[m][d])==0.0:
							DisCal=DisCal+ abs(float(sampleData[s][d])-0.0)
						else:
							DisCal=DisCal+ abs(float(sampleData[s][d])-float(Model[m][d]))
					E_c_dis=E_c_dis+(membershipList[m][s]*DisCal)
		if (flag==1):
			compact=compact+(E_c_dis/(points+0.0))
			flag=0 
	#print 'ss',compact/(NoCluster+0.0)
	return compact/(NoCluster+0.0)


def totalCulster(mem_array):
	count=0
	for i in range (len(mem_array)):
		if (np.count_nonzero(mem_array[i])!=0):
			count=count+1

	#print 'count',count
	return count


def init_genome(T,c,d,x_max):
 	for i in range(0,T.getn()):
		r = random.random() 
		if r>0.33:
			g=0
			T.setg(i,g)
		else:
			g=1
			T.setg(i,g)

	for i in range(0,T.getn()):
		C = random.randint(1,c)
		T.setc(i,C)

	for i in range(0,T.getn()):	
		D = random.randint(1,d)
		T.setd(i,D)
		x=x_max[D-1]
		j = random.randint(-1000,1000)
		T.setx(i,j*(x/1000))		
		


def cluster_center(T):
	k=0
	old=[]
	centers=[]
	for i in range(0,T.getn()):
		G=T.getg(i)
		if(G!=0):
			if(T.getc(i) not in old):
				centers.append([])
				centers[k].append(T.getc(i))
				for j in range(d):
					centers[k].append([])
					centers[k][j+1] = None
				D=int(T.getd(i))
				centers[k][D]=T.getx(i)
				k=k+1
				old.append(T.getc(i))
			else:
				indx =old.index(T.getc(i))
				D=int(T.getd(i))
				if(centers[indx][D]==None):
					centers[indx][D]=T.getx(i)
				else:
					centers[indx][D]=(centers[indx][D]+T.getx(i))/2
	centers = np.array(centers)
	return centers


def MembershipDegree(c,s,fuzzy):
	C = len(c)
	S = len(s)
	mem_array = np.zeros((C,S),dtype=int)
	#for i in range(0,C):
	#	mem_array[i][0]=c[i][0]
    

	for j in range(0,S):
		sample = np.array(s[j])
		center = np.zeros(d)
		for l in range(0,d):
			if(c[0][l+1] != None):
				center[l]=c[0][l+1]
		minm = np.sum(np.absolute(np.subtract(center,sample)))
		ind = 0
		for k in range(0,C):
			sample = np.array(s[j])
			center = np.zeros(d)
			for l in range(0,d):
				if(c[k][l+1] != None):
					center[l]=c[k][l+1]
			

			min_temp = np.sum(np.absolute(np.subtract(center,sample)))
			if(min_temp<minm):
				minm=min_temp
				ind=k
		#print 'ind' ,ind
		mem_array[ind][j]=1
        
        
		if(fuzzy==1):        
			threshold=minm+threshold_fuzzy_Value*minm        #fuzzy clustering
			for i in range(0,C):
				sample = np.array(s[j])
				center = np.zeros(d)
				for l in range(0,d):
					if(c[i][l+1] != None):
						center[l]=c[i][l+1]
				min_temp = np.sum(np.absolute(np.subtract(center,sample)))
				if(min_temp <= threshold):
					mem_array[i][j]=1
                    
	return mem_array


def convert_model(c):
	temp_array = np.delete(c, 0, 1)
	C = len(temp_array)
	D = d
	for i in range(C):
		for j in range(D):
			if(temp_array[i][j]== None):
				temp_array[i][j] = 0
	return temp_array


# def mutation_point_subt(T):
# 	i = random.randint(0,T.getn()-1) 
# 	k = random.randint(1,4)
# 	#print 'i',i
# 	#print 'k',k

# 	if(k==1):
# 		val= random.randint(0,1)
# 		T.setg(i,val)
# 	elif(k==2):
# 		val = random.randint(1,c)
# 		T.setc(i,val)
# 	elif(k==3):
# 		val = random.randint(1,d)
# 		T.setd(i,val)
# 	else:
# 		D = int(T.getd(i))
# 		x = x_max[D-1]
# 		j = random.randint(-1000,1000)
# 		T.setx(i,j*(x/1000))

def mutation_point_subt(T, no_of_pt_mut):

	for m in range(10*no_of_pt_mut):

		i=random.randint(1,T.getn()-1)
		j=1

		if(j==1): 
			k = random.randint(1,4)

			if(k==1):
				val= random.randint(0,1)
				T.setg(i,val)

			elif(k==2):
				val = random.randint(1,c)
				T.setc(i,val)

			elif(k==3):
				val = random.randint(1,d)
				T.setd(i,val)

			else:
				D = int(T.getd(i))
				x = x_max[D-1]
				j = random.randint(-1000,1000)
				T.setx(i,j*(x/1000))
	#print T.get_array()
	return T

def merge_operator_M(Y1,Y2,K):
	if(K==1):
		Y2[0]=Y1[0]
		return Y2
	elif(K==2):
		Y2[0]=Y1[0]
		Y2[1]=Y1[1]
		return Y2
	elif(K==3):
		Y1[3]=Y2[3]
		return Y1
	else:
		return Y1


def bionomial(L):
	mu=0.0005
	s= np.random.binomial(L, mu)
	return s

def large_deletions(T1):
	array = T1.get_array()

	while(1):
		i =  random.randint(0,T1.getn()-1)
		j =  random.randint(0,T1.getn()-1)
		if(i<=j):
			if((j-i)<=0.25*T1.getn()):
				break
		else:
			if((i-j)>=0.75*T1.getn()):
				break

	k = random.randint(1,4)
	M=[]

	if(i<=j):
		row_b=[]
		row_a=[]
		for x in range(0,i-1):
			row_b.append(x)
		for y in range(j+1,len(array) ):
			row_a.append(y)
		start = array[row_b, :]
		end = array[row_a, :]
		#start , end = T1.delete(i,j)
		temp = merge_operator_M(T1.get_row(i),T1.get_row(j),k)
		M.append([])
		for z in range(4):
			M[0].append(temp[z])
		#print 'start',start
		#print 'end', end
		#print 'M' , M
		#print 'len',len(start[0]),len(end[0]),len(M)
		final_array = np.concatenate((start,M,end), axis=0)
		T1.set_array(final_array)
		return T1
	else:
		row_m=[]
		for z in range(j+1,i-1):
			row_m.append(z)
		mid = array[row_m, :]
		#mid = T1.delete(i,j)
		temp = merge_operator_M(T1.get_row(i),T1.get_row(j),k)
		M.append([])
		for z in range(4):
			M[0].append(temp[z])
		#print 'mid', mid 
		#print 'M' , M
		#print 'len',len(mid[0]),len(M)
		final_array = np.concatenate((mid,M), axis=0)
		T1.set_array(final_array)
		return T1


# In[13]:


def large_duplications(T1):
	array = T1.get_array()
	while(1):
		i =  random.randint(0,T1.getn()-1)
		j =  random.randint(0,T1.getn()-1)
		if(i<=j):
			if((j-i)<=0.25*T1.getn()):
				break
		else:
			if((i-j)>=0.75*T1.getn()):
				break

	p =  random.randint(0,T1.getn()-1)
	k = random.randint(1,4)
	#print ' i j p k are', i ,j,p,k
	a=[]
	b=[]
	row_s=[]
	row_e=[]
	for x in range(0,p-1):
		row_s.append(x)
	for z in range(p+1,T1.getn()):
		row_e.append(z)
	start = array[row_s, :]
	end = array[row_e, :]
	
	if(i<=j):
		row_m=[]
		for y in range(i+1,j-1):
			row_m.append(y)
		mid = array[row_m, :]
		
		temp_a = merge_operator_M(T1.get_row(p),T1.get_row(i),k)
		temp_b = merge_operator_M(T1.get_row(j),T1.get_row(p),k)
		a.append([])
		for m in range(4):
			a[0].append(temp_a[m])
		b.append([])
		for n in range(4):
			b[0].append(temp_b[n])
		final_array = np.concatenate((start,a,mid,b,end), axis=0)
		T1.set_array(final_array)
		return T1

	else:

		temp_a = merge_operator_M(T1.get_row(p),T1.get_row(j),k)
		temp_b = merge_operator_M(T1.get_row(i),T1.get_row(p),k)
		a.append([])
		for m in range(4):
			a[0].append(temp_a[m])
		b.append([])
		for n in range(4):
			b[0].append(temp_b[n])

		row_m_1=[]
		row_m_2=[]
		for y in range(i+1,T1.getn()):
			row_m_1.append(y)

		for w in range(0,j-1):
			row_m_2.append(w)

		mid_1 = array[row_m_1, :]
		mid_2 = array[row_m_2, :]

		final_array = np.concatenate((start,a,mid_1,mid_2,b,end), axis=0)
		T1.set_array(final_array)
		return T1		


def large_translocations(T1):
	array = T1.get_array()
	while(1):
		i =  random.randint(0,T1.getn()-1)
		j =  random.randint(0,T1.getn()-1)
		if(i<=j):
			if((j-i)<=0.25*T1.getn()):
				break
		else:
			if((i-j)>=0.75*T1.getn()):
				break

	if (i<=j):	
		while(1):
			p =  random.randint(0,T1.getn()-1)
			if (i<=p<=j):
				pass
			else:
				break
	else:
		while(1):
			p =  random.randint(0,T1.getn()-1)
			if (j<p<i):
				break
					
	k = random.randint(1,4)
	#print ' i j p k are', i ,j,p,k


	row_s=[]
	row_e=[]
	row_ms=[]
	row_me=[]
	a=[]
	b=[]
	
	if(i<=j):
		if(0<=p<=i):
			for x in range(0,p-1):
				row_s.append(x)
			for z in range(p+1,i):
				row_e.append(z)

			for u in range(i,j):
				row_ms.append(u)
			for v in range(j,T1.getn()):
				row_me.append(v)

			start = array[row_s, :]
			end = array[row_e, :]
			mid1 = array[row_ms, :]
			mid2 = array[row_me, :]

			#print len(array), i, j, T1.get_row(i),T1.get_row(j), k	
			temp_a = merge_operator_M(T1.get_row(p),T1.get_row(i),k)
			temp_b = merge_operator_M(T1.get_row(p),T1.get_row(j),k)

			a.append([])
			for m in range(4):
				a[0].append(temp_a[m])
			b.append([])
			for n in range(4):
				b[0].append(temp_b[n])
			final_array = np.concatenate((start, a, mid1, end, b, mid2), axis=0)
			T1.set_array(final_array)
			#print  'c11', T1.getn()
			return T1

		elif(j<=p<=T1.getn()):
			for x in range(0,i):
				row_s.append(x)
			for z in range(i,j):
				row_e.append(z)

			for u in range(j+1,p-1):
				row_ms.append(u)
			for v in range(p,T1.getn()):
				row_me.append(v)

			start = array[row_s, :]
			end = array[row_e, :]
			mid1 = array[row_ms, :]
			mid2 = array[row_me, :]


			temp_a = merge_operator_M(T1.get_row(i),T1.get_row(p),k)
			temp_b = merge_operator_M(T1.get_row(j),T1.get_row(p),k)

			a.append([])
			for m in range(4):
				a[0].append(temp_a[m])
			b.append([])
			for n in range(4):
				b[0].append(temp_b[n])
			final_array = np.concatenate((start, a, mid1, end, b, mid2), axis=0)
			T1.set_array(final_array)
			#print  'c12', T1.getn()
			return T1


	else:
		
		for x in range(p,i):
			row_s.append(x)
		for z in range(j,p):
			row_e.append(z)

		start = array[row_s, :]
		end = array[row_e, :]	



		temp_a = merge_operator_M(T1.get_row(p),T1.get_row(j),k)
		temp_b = merge_operator_M(T1.get_row(i),T1.get_row(p),k)
		a.append([])
		for m in range(4):
			a[0].append(temp_a[m])
		b.append([])
		for n in range(4):
			b[0].append(temp_b[n])

		row_m_1=[]
		row_m_2=[]
		for y in range(i+1,T1.getn()):
			row_m_1.append(y)

		for w in range(0,j-1):
			row_m_2.append(w)

		mid_1 = array[row_m_1, :]
		mid_2 = array[row_m_2, :]

		final_array = np.concatenate((end,a,mid_1,mid_2,b,start), axis=0)
		T1.set_array(final_array)
		return T1		




def population(fuzzy):
#	print 'population'
	for i in range(N):
		centers_array.append(cluster_center(genome_array[i]))
		member_array.append(MembershipDegree(centers_array[i],scaled,fuzzy))
		model_array.append(convert_model(centers_array[i]))
		#print 'ca',len(centers_array[i]) ###for checking

		mnd_array.append([])
		fnd_array.append([])
		compactness_array.append([])

		# fnd = feature_non_redundency(centers_array[i])
		# fnd_array[i].append(fnd)

		mnd_array[i].append(membership_non_redundency(member_array[i]))
		fnd_array[i].append(Calculate_PSM(model_array[i],scaled))
		#fnd_array[i].append(feature_non_redundency(centers_array[i]))
		compactness_array[i].append(Calculate_compactness(model_array[i],scaled,member_array[i]))

	

def population_after_mutation(fuzzy):
# #	print 'population after mut'
# 	for j in range(N):
# 		#print j
# 		T1_1 = copy.deepcopy(genome_array[j])
# 		mut_prob = random.random()

# 		if(mut_prob<=lar_del):
# 			large_deletions(T1_1)
# 		if(lar_dup_low<mut_prob<=lar_dup_Up):
# 			large_duplications(T1_1)

# 		no_of_pt_mut = int(T1_1.getn()*pt_mut)
		
# 		for k in range(no_of_pt_mut):
# 			mutation_point_subt(T1_1)

# 		genome_array_1.append(T1_1)


	for j in range(N):
		#print j
		T1_1 = copy.deepcopy(genome_array[j])
		no_of_pt_mut = bionomial(T1_1.getn())


		for i in range(max(1,no_of_pt_mut)): 
			num = random.randint(0,5)
			#print 'num',num

			if (num==0):
				#print 'z0',T1_1.getn()
				T1=large_deletions(T1_1)
				#print 'z0',T1.getn()
				T2=large_duplications(T1)
				#print 'z0',T2.getn()
				T3=large_translocations(T2)
				#print 'z0',T3.getn()
				T1_1=T3
				#print 'z',T3.getn()

			elif (num==1):
				#print 'z1',T1_1.getn()
				T1=large_deletions(T1_1)
				#print 'z1',T1.getn()
				T2=large_translocations(T1)
				#print 'z1',T2.getn()
				T3=large_duplications(T2)
				#print 'z1',T3.getn()
				T1_1=T3
				#print 'z',T3.getn()

			elif (num==2):
				#print 'z2',T1_1.getn()
				T1=large_duplications(T1_1)
				#print 'z2',T1.getn()
				T2=large_deletions(T1)
				#print 'z2',T2.getn()
				T3=large_translocations(T2)
				#print 'z2',T3.getn()
				T1_1=T3
				#print 'z',T3.getn()

			elif (num==3):
				#print 'z3',T1_1.getn()
				T1=large_duplications(T1_1)
				#print 'z3',T1.getn()
				T2=large_translocations(T1)
				#print 'z3',T2.getn()
				T3=large_deletions(T2)
				#print 'z3',T3.getn()
				T1_1=T3
				#print 'z',T3.getn()
			

			elif (num==4):
				#print 'z4',T1_1.getn()
				T1=large_translocations(T1_1)
				#print 'z4',T1.getn()
				T2=large_deletions(T1)
				#print 'z4',T2.getn()
				T3=large_duplications(T2)
				#print 'z4',T3.getn()
				T1_1=T3
				#print 'z',T3.getn()
			

			else:
				#print 'z5',T1_1.getn()
				T1=large_translocations(T1_1)
				#print 'z5',T1.getn()
				T2=large_duplications(T1)
				#print 'z5',T2.getn()
				T3=large_deletions(T2)
				#print 'z5',T3.getn()
				T1_1=T3
				#print 'z',T3.getn()


		no_of_pt_mut = bionomial(T3.getn()) #int(T1_1.getn()*0.3)
		T4=mutation_point_subt(T3, max(1,no_of_pt_mut))
		genome_array_1.append(T4)

		if t==T-1:
			print 'Genome_Length', T4.getn()


	for l in range(N):
		centers_array_1.append(cluster_center(genome_array_1[l]))
		#print len(centers_array_1[l])
		if(len(centers_array_1[l])==0):
			genome_array_1[l].print_genome()
			print 'halt'
		member_array_1.append(MembershipDegree(centers_array_1[l],scaled,fuzzy))
		model_array_1.append(convert_model(centers_array_1[l]))

		mnd_array_1.append([])
		fnd_array_1.append([])
		compactness_array_1.append([])
        
		#fnd = feature_non_redundency(centers_array_1[l])
		mnd_array_1[l].append(membership_non_redundency(member_array_1[l]))
		fnd_array_1[l].append(Calculate_PSM(model_array_1[l],scaled))
		#fnd_array_1[l].append(feature_non_redundency(centers_array_1[l]))
		compactness_array_1[l].append(Calculate_compactness(model_array_1[l],scaled,member_array_1[l]))

def non_dominating(P):
	population1=[]
	population1=P

	F_one_set=[]
	for i in range(2*N): #population_size=N
		F_one_set.append([])
	ndf, dl, dc, ndr = pg.fast_non_dominated_sorting(points=population1) #ndf, dl, dc,
	# for i in range(2*N):
	# 	F_one_set[ndr[i]].append(population1[i])
	# F_one_set2 = [x for x in F_one_set if x != []]
	# return F_one_set2
	return ndf


#Functions for Evaluation Metrics

def _validity_cluster_checking(found_clusters_effective,
                               threshold_cluster_validity=0.0):
    return found_clusters_effective >= threshold_cluster_validity

def _mapped(contingency_table):
    mapped_clusters = (contingency_table.T * 1. / contingency_table.sum(1)).T
    return mapped_clusters == mapped_clusters.max(0)

def compute_only_f1(contingency_table,
                    valid_clusters,
                    mapped_clusters):
    num = mapped_clusters * contingency_table
    num = num.loc[:, valid_clusters]
    num = num.sum(1)
    denum_recall = contingency_table.sum(1)
    rec = num * 1. / (denum_recall + NOT_ZERO_DIVISION_SECURITY)
    denum_precision = mapped_clusters * contingency_table.sum(0) * valid_clusters
    denum_precision = denum_precision.sum(1)
    precis = num * 1. / (denum_precision + NOT_ZERO_DIVISION_SECURITY)
    denum = rec + precis
    num = 2 * rec * precis
    return sum(num * 1.0 / (denum + NOT_ZERO_DIVISION_SECURITY)) * 1. / (len(num) + NOT_ZERO_DIVISION_SECURITY)

def f1(contingency_table,
       threshold_cluster_validity=0.0):
    mapped_clusters = _mapped(contingency_table)
    found_clusters_effective = contingency_table.sum(0)
    valid_clusters = _validity_cluster_checking(found_clusters_effective, threshold_cluster_validity)
    return compute_only_f1(contingency_table, valid_clusters, mapped_clusters)

def compute_only_entropy(contingency_table,
                         valid_clusters):
    contingency_table = contingency_table.loc[:, valid_clusters]
    found_clusters_effective = contingency_table.sum(0)
    p_h_in_c = contingency_table * 1. / (found_clusters_effective + NOT_ZERO_DIVISION_SECURITY)
    
    log_p_h_in_c = np.log(p_h_in_c)

    pre_ec = -1. * p_h_in_c * log_p_h_in_c

    pre_ec = pre_ec.fillna(0)
    #print pre_ec
    ec = pre_ec.sum(0)
    num = (ec * found_clusters_effective).sum()
    denum = found_clusters_effective.sum() * np.log(len(contingency_table.index))
    return 1. - (num * 1. / denum)

def entropy(contingency_table,
            threshold_cluster_validity=0.0):
    #contingency_table = pd.crosstab(cluster_hidden, cluster_found)
    valid_clusters = _validity_cluster_checking(contingency_table.sum(0), threshold_cluster_validity)
    return compute_only_entropy(contingency_table, valid_clusters)

def compute_only_accuracy(contingency_table, valid_clusters, found_clusters_effective):
    best_matching_hidden_cluster = contingency_table == contingency_table.max(0)
    best_matching_hidden_cluster_weight = 1. / best_matching_hidden_cluster.sum(0)
    correctly_predicted_objects = contingency_table * best_matching_hidden_cluster * best_matching_hidden_cluster_weight
    correctly_predicted_objects *= valid_clusters
    return sum(correctly_predicted_objects.sum(0)) * 1. / (sum(found_clusters_effective)+NOT_ZERO_DIVISION_SECURITY)

def accuracy(contingency_table, threshold_cluster_validity=0.0):
    #contingency_table = pd.crosstab(cluster_hidden, cluster_found)
    found_clusters_effective = contingency_table.sum(0)
    valid_clusters = _validity_cluster_checking(found_clusters_effective, threshold_cluster_validity)
    return compute_only_accuracy(contingency_table, valid_clusters, found_clusters_effective)

# Start of Functions Computes a max weight perfect matching in a bipartite graph
def improveLabels(val):
    """ change the labels, and maintain minSlack.
    """
    for u in S:
        lu[u] -= val
    for v in V:
        if v in T:
            lv[v] += val
        else:
            minSlack[v][0] -= val

def improveMatching(v):
    """ apply the alternating path from v to the root in the tree.
    """
    u = T[v]
    if u in Mu:
        improveMatching(Mu[u])
    Mu[u] = v
    Mv[v] = u

def slack(u,v): return lu[u]+lv[v]-w[u][v]

def augment():
    """ augment the matching, possibly improving the lablels on the way.
    """
    while True:
        # select edge (u,v) with u in S, v not in T and min slack
        ((val, u), v) = min([(minSlack[v], v) for v in V if v not in T])
        assert u in S
        if val>0:
            improveLabels(val)
        # now we are sure that (u,v) is saturated
        assert slack(u,v)==0
        T[v] = u                            # add (u,v) to the tree
        if v in Mv:
            u1 = Mv[v]                      # matched edge,
            assert not u1 in S
            S[u1] = True                    # ... add endpoint to tree
            for v in V:                     # maintain minSlack
                if not v in T and minSlack[v][0] > slack(u1,v):
                    minSlack[v] = [slack(u1,v), u1]
        else:
            improveMatching(v)              # v is a free vertex
            return

def maxWeightMatching(weights):
    """ given w, the weight matrix of a complete bipartite graph,
        returns the mappings Mu : U->V ,Mv : V->U encoding the matching
        as well as the value of it.
    """
    global U,V,S,T,Mu,Mv,lu,lv, minSlack, w
    w  = weights
    n  = len(w)
    U  = V = range(n)
    lu = [ max([w[u][v] for v in V]) for u in U]  # start with trivial labels
    lv = [ 0                         for v in V]
    Mu = {}                                       # start with empty matching
    Mv = {}
    while len(Mu)<n:
        free = [u for u in V if u not in Mu]      # choose free vertex u0
        u0 = free[0]
        S = {u0: True}                            # grow tree from u0 on
        T = {}
        minSlack = [[slack(u0,v), u0] for v in V]
        augment()
    #                                    val. of matching is total edge weight
    val = sum(lu)+sum(lv)
    return (Mu, Mv, val)

def compute_CE_RNIA(Total_Cluster, Model_Selected, PredCluster, trueCluster):

	dimensionSet=[]
	PredLevelModel=[]
	dimensionSetIndex=[]

	for i in range (Total_Cluster):
		dimensionSetIndex.append([])

	for s in range (len(mydata_list_of_list)):
		for m in range (len(Model_Selected)):
			if membershipList[m][s]==1:
				PredLevelModel.append(m)
				break

	for i in range (Total_Cluster):
		X=PredCluster[i][0]
		M=PredLevelModel[X]
		dimensionSet.append(Model_Selected[M])

	for i in range (Total_Cluster):
		for d in range (dimension):
			if dimensionSet[i][d]!=0:
				dimensionSetIndex[i].append(d)

	size=max(len(trueCluster),len(PredCluster))

	confusionMatrix=[]

	for i in range(size):
		confusionMatrix.append([])
		for j in range(size):
			confusionMatrix[i].append(0)

	for i in range (len(PredCluster)):
		for j in range (len(trueCluster)):
			same=len(set(PredCluster[i]) & set(trueCluster[j]))
			common= len(set(dimensionSetIndex[i]) & set(trueSubspace[j]))
			confusionMatrix[i][j]=common*same

	confusionMatrix = [[confusionMatrix[j][i] for j in range(len(confusionMatrix))] for i in range(len(confusionMatrix[0]))]
	D_max= (list(maxWeightMatching(confusionMatrix))[-1])
	return D_max, confusionMatrix, dimensionSetIndex



# def calculate_Union(PredCluster, PredSubspace): #Only for non-overlapping but not for overlapping

# 	Pred_All_U=[]
# 	True_All_U=[]
# 	for i in range(len(PredCluster)):
# 		for k in itertools.product(PredCluster[i],PredSubspace[i]):
# 			Pred_All_U.append(k)
# 	Pred_All_U=set(Pred_All_U)

# 	for i in range(len(trueCluster)):
# 		for k in itertools.product(trueCluster[i],trueSubspace[i]):
# 			True_All_U.append(k)
# 	True_All_U=set(True_All_U)

#    	All_U= Pred_All_U.union(True_All_U)  #set(j).union(set(k))
#    	U=len(All_U)
#    	return U


def calculate_Union(PredCluster, PredSubspace): # More Accurate for both overlapping and non overlapping

	Pred_All_U=[]
	True_All_U=[]

	for i in range(len(PredCluster)):
		for k in itertools.product(PredCluster[i],PredSubspace[i]):
			Pred_All_U.append(k)

	for i in range(len(trueCluster)):
		for k in itertools.product(trueCluster[i],trueSubspace[i]):
			True_All_U.append(k)

	for i in range (len(True_All_U)):
		if True_All_U[i] not in Pred_All_U:
			Pred_All_U.append(True_All_U[i])

   	U=len(Pred_All_U)
   	return U


def No_of_Clusters(M,member):
	Model=M
	membershipList=member
	clusters=0
	for m in range(len(Model)):
		if np.count_nonzero(Model[m]!=0) and np.count_nonzero(membershipList[m]!=0):
			clusters=clusters+1
	return clusters

def evaluation(Model_Selected,membershipList):
	membershipListNonZero=[]
	PredLevel=[]
	PredCluster=[]
	Model_Selected=np.asarray(Model_Selected)
	membershipList=np.asarray(membershipList)
	Total_Cluster= totalCulster(membershipList)
	Model_Selected=Model_Selected.tolist()

	for i in range (len(membershipList)):
		if np.count_nonzero(membershipList[i]!=0):
			membershipListNonZero.append(membershipList[i])

	# for s in range (len(mydata_list_of_list)):
	# 	for m in range (len(membershipListNonZero)): 
	# 		if membershipListNonZero[m][s]==1:
	# 			PredLevel.append(m)
	# 			break

	for i in range(Total_Cluster):
		PredCluster.append([])
	
	for m in range (len(membershipListNonZero)): 
		for i in range(len(mydata_list_of_list)):
			if(membershipListNonZero[m][i]==1):
				PredCluster[m].append(i)
	#print PredCluster

	#exit()
	length=[]
	for i in range(len(PredCluster)):
		length.append(len(PredCluster[i]))

	contingency_table=[]

	for i in range(len(trueCluster)):
		contingency_table.append([])

	for i in range (len(trueCluster)):
		truesame=[]
		for j in range (len(PredCluster)):
			truesame.append(len(set(trueCluster[i]) & set(PredCluster[j])))
			contingency_table[i]=truesame
	contingency_table1 = pd.DataFrame(contingency_table)
	contingency_table=contingency_table1
	D_max, confusionMatrix, PredSubspace=compute_CE_RNIA(Total_Cluster, Model_Selected, PredCluster, trueCluster)

	F_measure=f1(contingency_table, threshold_cluster_validity=0.0)
	Entropy = entropy(contingency_table, threshold_cluster_validity=0.0)
	Accuracy=accuracy(contingency_table,threshold_cluster_validity=0.0)
	I=0
	for i in range(Total_Cluster):
		I=I+sum(confusionMatrix[i])
	#U=calculate_Union(PredCluster, PredSubspace)  #For non-overlap Clustering
	U=calculate_Union(PredCluster, PredSubspace) #For overlap  clustering
	U1=dataCount*dimension # U value will be different for syntetic dataset
	
	CE=(float)(D_max)/(U+0.0)
	RNIA=(float)(I)/(U+0.0)

	Avgdim=0.0
	for i in range (Total_Cluster):
		Avgdim=Avgdim+len(PredSubspace[i])
	Avgdim=Avgdim/(Total_Cluster+0.0)
	return F_measure, Accuracy, CE, RNIA, Entropy, Avgdim, Total_Cluster


#***************************************************************************  Main Program Starts  **********************************************************************************




##########-mutation-parameters-#############
# lar_del=0.2
# lar_dup_low=0.2
# lar_dup_Up=0.5
# pt_mut=0.2

##############################
trueSubspace=[]
trueSubspace1=[]
trueCluster=[]
scaled=[]

dataframe = read_table(filename, sep=',',header=None)
data = np.array(dataframe)
n_sample = len(data)      #length of sample 
dataCount=n_sample
d1 = len(data[0])
data = data[:,0:d1-1]
d=len(data[0])
dimension=d

SS = pp.StandardScaler(copy=True, with_mean=True, with_std=True)   
scaled_All = SS.fit_transform(data) 

# geometric_center = np.sum(scaled_All,axis =0)
# print"geometric_center", geometric_center

x_max = np.amax(scaled_All, axis=0)

old=[n_sample]
for i in range(0,sample_size):
	while True:
		k=(random.randint(0,n_sample-1))
		#print 'k=',k
		if k not in old:
			old.append(k)
			break
	scaled.append(scaled_All[k])


readfile=open(filenametrue,'r') 
for line in readfile.readlines():
	line = line.strip()
	true_list=[]
	for word in line.split(' '):
		true_list.append(word)
	true_list = [int(x) for x in true_list]
	trueSubspace1.append(true_list[:d])
	trueCluster.append(true_list[dimension+1:]) 
trueSubspace1=np.array(trueSubspace1)

for i in range (len(trueSubspace1)):
	index=np.where(trueSubspace1[i] == 1)[0]
	trueSubspace.append(index)
trueSubspace = np.array(trueSubspace).tolist()
readfile.close()



#####genome_init#####
genome_array =[]
for i in range(N):
	T1_initial = genome(int(genome_size+20))
	init_genome(T1_initial,c,d,x_max)
	genome_array.append(T1_initial)

#print T1_initial.get_array()


#####generations#####

for t in range(T):


	if(t%100==0):
		print 'Iteration',t
	
	randSample=(random.randint(0,sample_size-1))
	randData=(random.randint(0,n_sample-1))
	scaled[randSample]=scaled_All[randData]

	centers_array =[]
	model_array = []
	member_array = []
	mnd_array = []
	fnd_array = []
	compactness_array = []

	population(1)
	obj_funs_array = np.concatenate((mnd_array,fnd_array,compactness_array), axis=1)
	#print obj_funs_array

	genome_array_1 =[]
	centers_array_1 =[]
	model_array_1 = []
	member_array_1 = []
	mnd_array_1 = []
	fnd_array_1 = []
	compactness_array_1 = []

	population_after_mutation(1)

	obj_funs_array_1 = np.concatenate((mnd_array_1,fnd_array_1,compactness_array_1), axis=1)
	
	Total_obj_funs_array = np.concatenate((obj_funs_array,obj_funs_array_1), axis=0)
	# nds = pg.sort_population_mo(Total_obj_funs_array)
	# #print Total_XB_PBM_array
	#print Total_obj_funs_array


	F=non_dominating(Total_obj_funs_array)
	#print F
	#GenerateSolution(F)

	new_genome = []

	K=N
	indexList=[]
	for i in range(len(F)):
		if(len(F[i])<=K):
			K=K-len(F[i])
			for j in range(len(F[i])):
				indexList.append(F[i][j])
			if (K==0):
				break

		else:
			for j in range(K):
				indexList.append(F[i][j])
			K=0
			break

	for i in range(len(indexList)):
		if(indexList[i]>=N):
			T1_1_1 = copy.deepcopy(genome_array_1[int(indexList[i]-N)])
			new_genome.append(T1_1_1)
					
		else:
			T1_1_1 = copy.deepcopy(genome_array[indexList[i]])
			new_genome.append(T1_1_1)

	# new_genome = []
	# for i in range(N):
	# 	if(nds[i]>=N):
	# 		#print i 
	# 		#print nds[i]
	# 		T1_1_1 = copy.deepcopy(genome_array_1[int(nds[i]-N)])
	# 		new_genome.append(T1_1_1)
	# 	else:
	# 		T1_1_1 = copy.deepcopy(genome_array[nds[i]])
	# 		new_genome.append(T1_1_1)

	for j in range(N):
		T1_dup = copy.deepcopy(new_genome[j])
		genome_array[j]=T1_dup


Total_Time=(time.clock() - start_time)


#####final_generation####

centers_array =[]
model_array = []
member_array = []
mnd_array = []
fnd_array = []
compactness_array = []
#ss=dataCount
population(1)

Model_dict=model_array
mydata_list_of_list=scaled_All

SDmax=c
NoOfIteration=T
print 'Iteration=', NoOfIteration
#Subspace Clusters Evaluation

F_measure_All=[]
Accuracy_All=[]
CE_All=[]
RNIA_All=[]
Entropy_All=[]
Avgdim_All=[]
Total_Cluster_All=[]
#print 'T',T
for i in range(len(Model_dict)):
	membershipList= MembershipDegree(centers_array[i],mydata_list_of_list,1)
	#member=totalCulster(membershipList)
	#print 'len(Model_dict[i])', len(membershipList)
	
	F_measure, Accuracy, CE, RNIA, Entropy, Avgdim, Total_Cluster= evaluation(Model_dict[i],membershipList)

	F_measure_All.append(F_measure)
	Accuracy_All.append(Accuracy)
	CE_All.append(CE)
	RNIA_All.append(RNIA)
	Entropy_All.append(Entropy)
	Avgdim_All.append(Avgdim)
	Total_Cluster_All.append(Total_Cluster)



print 
print 'sample_size=',sample_size
print 'SDmax=', SDmax

print 
print 'F_measure: MAX, MIN, ALL', max(F_measure_All), min(F_measure_All), F_measure_All
print 
print 'Accuracy: MAX, MIN, ALL', max(Accuracy_All), min(Accuracy_All), Accuracy_All
print 
print 'CE: MAX, MIN, ALL', max(CE_All), min(CE_All), CE_All
print 
print 'RNIA: MAX, MIN, ALL', max(RNIA_All), min(RNIA_All), RNIA_All
print 
print 'Entropy: MAX, MIN, ALL', max(Entropy_All), min(Entropy_All), Entropy_All
print 
print 'Coverage=', 1
print 
print 'NumClusters: MAX, MIN, ALL', max(Total_Cluster_All), min(Total_Cluster_All), Total_Cluster_All
print 
print 'AvgDim: MAX, MIN, ALL', max(Avgdim_All), min(Avgdim_All), Avgdim_All
print 
print 'run_Time=',Total_Time
print 

print round(max(F_measure_All), 2), round(min(F_measure_All), 2), round(max(Accuracy_All), 2), round(min(Accuracy_All), 2), round(max(CE_All), 2), round(min(CE_All), 2), round(max(RNIA_All), 2), round(min(RNIA_All), 2), round(max(Entropy_All), 2), round(min(Entropy_All), 2), 1.00, 1.00, round(max(Total_Cluster_All), 2), round(min(Total_Cluster_All), 2), round(max(Avgdim_All), 2), round(min(Avgdim_All), 2), Total_Time, Total_Time
