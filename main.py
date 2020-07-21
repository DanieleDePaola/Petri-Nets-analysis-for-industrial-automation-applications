import numpy as np
from scipy.linalg import lu
import scipy.optimize
from numpy.linalg import svd

def nullspace(A, atol=1e-13, rtol=0):
    A = np.atleast_2d(A)
    u, s, vh = svd(A)
    tol = max(atol, rtol * s[0])
    nnz = (s >= tol).sum()
    ns = vh[nnz:].conj().T
    return ns



# Input and output matrix acquisition

print("Insert the number of places: ");
#p=input();  
p=int(6);
print("Insert the number of transitions: ");
#t=input();
t=int(5);  

print("Send me c if you want to send me the c matrix, io if you want to put the input and output matrices")


choice=input(); 
while(choice!="c" and choice!="io"):
	print("Wrong choice, put c or io")
	choice = input(); 

#if(choice == "io"):
	
	#Output Matrix acquisition
	print("Let's create the Output Marix!");
	print(); 
	OutputMatrix=[]; 

	for i in range(int(p)):
		print("Values in line", i); 
		OutputLine= [];
		for j in range(int(t)): 
			print("Insert number: ")
			c=int(input());
			OutputLine.append(c); 
		
		OutputMatrix.append(OutputLine); 


	#Input matrix acquisition
	print("Let's create the Input Marix!");
	print(); 
	InputMatrix = []; 		
	for i in range(int(p)):
		print("Values in line", i); 
		InputLine= [];
		for j in range(int(t)): 
			print("Insert number: ")
			c=int(input());
			InputLine.append(c); 
		
		InputMatrix.append(InputLine); 


	# Now we can create the Cmatrix, the C matrix is given by  OutputMatrix - InputMatrix	

	print(); 
	print("OutputMatrix is", OutputMatrix); 
	print(); 
	print("Input Matrix is", InputMatrix); 

	print(); 
	print("The CMatrix is given  by OutputMatrix - InputMatrix"); 
	print()

	outArray = np.squeeze(np.asarray(OutputMatrix)); 
	inpArray = np.squeeze(np.asarray(InputMatrix)); 

	print("numpyOut", outArray); 

	CMatrix=np.subtract(outArray,inpArray);  
	
#elif(choice == "c"): 

	#C matrix acquisition
	print("Let's create the C Matrix!");
	print(); 
	CMatrix = []; 		
	for i in range(int(p)):
		print("Values in line", i); 
		CLine= [];
		for j in range(int(t)): 
			print("Insert number: ")
			c=int(input());
			CLine.append(c); 
		
		CMatrix.append(CLine); 

CMatrix=[[-1,-1,0,0,2],[-1,-1,1,1,0],[1,0,-1,0,0],[0,1,0,-1,0],[0,0,1,0,-1],[0,0,0,1,-1]]

CMatrix = np.squeeze(np.asarray(CMatrix)); 
print("The CMatrix is \n", CMatrix);

#Fsm or Marked Graph

#FSM if I have only one 1 and only one -1 in each column, it means that each transition exits from only one place and enters only one place
one = False; 
minusOne = False; 
fsm= True; 

for j in range(int(t)): 
	if(fsm==False): break 

	for i in range(int(p)):
		if(CMatrix[i][j]==1):
			 if(one): fsm=False; 
			 one=True; 

		if(CMatrix[i][j]==-1): 
			if(minusOne): fsm = False; 
			minusOne=True;
		if(CMatrix[i][j]>1 or CMatrix[i][j]<-1):fsm=False	
		one=False; 
		minusOne=False; 	

#TODO if(minusOne==False or one==False): fsm=False; #in case I found only one or minus one in on ecolumn
print()

#Marked Graph, the same procedure used above, but for rows. This property says that the if one place enters different transitions, all these transitions are not output for other places

one=False; 
minusOne=False; 
markedGraph=True; 
for i in range(int(p)): 
	if(markedGraph==False): break 

	for j in range(int(t)):	
		if(CMatrix[i][j]==1):
			 if(one): markedGraph=False; 
			 one=True; 

		if(CMatrix[i][j]==-1): 
			if(minusOne): markedGraph = False; 
			minusOne=True;
		if(CMatrix[i][j]>1 or CMatrix[i][j]<-1):markedGraph=False		
	one=False; 
	minusOne=False; 		

#TODO if(minusOne==False or one==False): markedGraph=False; #in case I found only one or minus one in one column

print("is Fsm: ", fsm); 
print("is Marked Graph:", markedGraph); 
print()

#Finding out T invariants solving a linear system in p equations. I'll try to use here the true power of numpy in matrix manipulation. T-invariants are studied to understand if the net is reversible or could not be reversible. 

print("T-invariants analysis")
b=np.zeros(int(p)); 


#CMatrix = np.squeeze(np.asarray(CMatrix)); 
#CMatrix=CMatrix.transpose(); 
#print("CMatrix transposed first time\n",CMatrix); 

#CMatrix = np.vstack([CMatrix, b]);
#print("CMatrix with zeros\n",CMatrix); 


#CMatrix=CMatrix.transpose(); 

#CMatrix=scipy.optimize._remove_redundancy._remove_redundancy(CMatrix, np.zeros_like(CMatrix[:, 0]))[0]
 
print("CMatrix system \n",CMatrix); #TODO Wrong transformations, work on it 

pl, u = lu(CMatrix, permute_l=True)
u=scipy.optimize._remove_redundancy._remove_redundancy(u, np.zeros_like(u[:, 0]))[0] #Remove equal rows


#print("Cmatrix equal lines optimized\n", u)

u=u[~(u==0).all(1)]#remove all zeros rows
print("Cmatrix all zeroes eliminated:\n", u)

kernelRank=np.linalg.matrix_rank(u, tol=None, hermitian=False)


if((t-kernelRank)>0): 
	print ("\nThe system has inf ^", (t-kernelRank))
	print("solutions, it means at least a T-invariant exixts, so the P-net COULD NOT be reversible, we should explore the marking graph to be sure about reversibility\n")
	#With more time I can implement a BFS here to explore the graph


#P-Invariant analysis, from theory, P-invariants are the generators of the Kernel space of for the transposed C matrix, if the sum of the n p-invariants is a p-invariant, the net is conservative (strictly conservative for p-inv=[1,1,1,1,1,1]). If conservative the net is also limited. The tokens in one place will never be more than k, with k <inf. 

CMatrix=CMatrix.transpose(); 
pl, pInvOptimization = lu(CMatrix, permute_l=True)
pInvOptimization=scipy.optimize._remove_redundancy._remove_redundancy(CMatrix, np.zeros_like(CMatrix[:, 0]))[0]
pInvOptimization=pInvOptimization[~(pInvOptimization==0).all(1)]
print("Cmatrix for the ")


ns=nullspace(CMatrix)
ns=scipy.optimize._remove_redundancy._remove_redundancy(ns, np.zeros_like(ns[:, 0]))[0]

print("Result for p inv:\n ", ns)


#print("basis columns:\n", basis_columns )








  

#CPinv= np.linalg.pinv(CMatrix); 
#TInvariants = CPinv.dot(b)


#print("The T invariants are: ", TInvariants); 


	