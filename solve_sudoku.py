# -*- coding: utf-8 -*-

#####
# Mengmiao Wu
###

import pdb #Use pdb.set_trace() to make a break in your code.
import numpy as np
import Queue

###################
# Solve sudoku #
###################

#####
# reduce: Fonction used by AC3 to reduce the domain of Xi using constraints from Xj.
#
# Xi: Variable (tuple (Y,X)) with reduced domain, if possible.
#
# Xj: Variable (tuple (Y,X)) with reduced domain, if possible.
#
# csp: Object of class CSP containing all information relative to constraint 
#      satisfaction for Sudoku.
#
# return: A tuple containing a boolean indicating if there were changes to the domain, and the csp.
### 

def reduce(Xi, Xj, csp):
	for x in csp.domains[Xi]:
		not_conflict=[]
		for y in csp.domains[Xj]:
			if x==y:
				not_conflict.append('no')
			else:
				not_conflict.append('yes')
		if 'yes' not in not_conflict:
			csp.domains[Xi].remove(x)
			return True,csp
	'''
		if len(csp.domains[Xj])==1:
			if x in csp.domains[Xj]:
				csp.domains[Xi].remove(x) 
				
	'''
	#print False
	return False,csp



#####
# AC3: Function used to reduce the domain of variables using AC3 algorithm.
#
# csp: Object of class CSP containing all information relative to constraint 
#      satisfaction for Sudoku.
#
# return: A tuple containing the optimized csp and a boolean variable indicating if there are no violated constraints.
### 

def AC3(csp):

	queue=Queue.Queue()
	for Xi in csp.variables:
		for Xj in csp.constraints[Xi]:
			queue.put((Xi,Xj))

	while not queue.empty():
		(x,y)=queue.get()
		#print (x,y)
		ifr,csp=reduce(x,y,csp)
		if ifr:
			#print "here1"
			if len(csp.domains[x])==0:
				return csp,False
			if len(csp.domains[x])==1:
				#print "here3"'''
				for n in csp.constraints[x]:
					#if ''.join(csp.domains[x]) in csp.domains[n]:
					queue.put((n,x))
	#print "here"
	return csp,True


#####
# is_compatible: Fonction verifying the correctness of an assignment.
#
# X: Tuple containing the position in y and in x of the cell concerned by the assignment.
#
# v: String representing the value (between [1-9]) affected by the assignment.
#
# assignment: dict mapping cells (tuple (Y,X)) to values.
#
# csp: Object of class CSP containing all information relative to constraint 
#      satisfaction for Sudoku.
#
# return: A boolean indicating if the assignment of the value v in cell X is legal.
### 

def is_compatible(X,v, assignment, csp):
	for c in csp.constraints[X]:
		
		#print ''.join(csp.domains[c])

		if ''.join(csp.domains[c])==v:
			return False
		
		if c in assignment:
			if assignment[c]==v: 
				return False
	return True


def select_next(assignments,csp):
	incr=1
	while incr != 10:
		for var in csp.variables:
			if var not in assignments:
				if len(csp.domains[var])==incr:
					#print "incr:"
					#print incr
					return var
		incr += 1

#####
# backtrack : Function used to find the missing assignments on the Sudoku grid using Backtracking Search.
#
#
# assignment: dict mapping cells (tuple (Y,X)) to values.
#
# csp: Object of class CSP containing all information relative to constraint 
#      satisfaction for Sudoku.
#
# retour: The dictionary of assignments (cell => value)
### 

def backtrack(assignments, csp):
	if len(assignments)==81:
		return assignments
	
	var=select_next(assignments,csp)	 #minimum-remaining-value
	#print var
	#print csp.domains[var]
	
	for value in csp.domains[var]:
		assignments[var]=value
		#print var
		#print value
		#print csp.domains[var]

		csp2=csp.copy()
		csp2.domains[var]=[value]

		if is_compatible(var,value,assignments,csp):
			csp2,respect=AC3(csp2)
			if respect:
				result=backtrack(assignments,csp2)				
				if result != None:
					#print "return result"
					return result
				#print "am i here"
		#print "deleting"
		#print assignments[var]
		del assignments[var]
		#print assignments
	#print "return assignments"
	return None

#####
# backtracking_search : Main function for backgracking
#
# csp: Object of class CSP containing all information relative to constraint 
#      satisfaction for Sudoku.
# The member variables are:
#      'variables'   : list of cases (tuple (Y,X)) 
#      'domaines'    : dict mapping a cell to a list of possible values
#      'contraintes' : dict mapping a cell to a list of cells who's value must be different from the first cell
#
# return: The dictionary of assignments (cell => value)
### 

def backtracking_search(csp):
	#print csp.domains[(0,1)]
	#reduce((0,1),(0,2),csp)
	#print csp.domains[(0,1)]
	#print csp.variables
	#print csp.domains
	#x='5'
	#csp.domains[(7,3)].remove(x)
	#print csp.domains
	#print csp.constraints
	return backtrack({}, csp.copy())
