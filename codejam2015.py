import numpy as np

def nonlin(x,deriv=False):
	if(deriv==True):
	    return x*(1-x)

	return 1/(1+np.exp(-x))
    
def neural_net(training_data_feature_matrix, training_data_output):
	X = training_data_feature_matrix               
	y = training_data_output
	np.random.seed(1)

	# randomly initialize our weights with mean 0
	#syn0 = 2*np.random.random((3,4)) - 1
	#syn1 = 2*np.random.random((4,1)) - 1
	syn0 = 2*np.random.random((265,265)) - 1
	syn1 = 2*np.random.random((265,1)) - 1

	# print X[0][0]
	# print X[0][1]
	# print type(X[0][0])
	# print type(syn0[0][0])
	# print 'syn0'
	# print syn0.shape

	#for j in xrange(1):
	for j in range(60000):

		# Feed forward through layers 0, 1, and 2
	    l0 = X
	    l1 = nonlin(np.dot(l0,syn0))
	    l2 = nonlin(np.dot(l1,syn1))

	    # how much did we miss the target value?
	    l2_error = y - l2
	    
	    if (j% 10000) == 0:
	        print ("Error:" + str(np.mean(np.abs(l2_error))))
	        
	    # in what direction is the target value?
	    # were we really sure? if so, don't change too much.
	    l2_delta = l2_error*nonlin(l2,deriv=True)

	    # how much did each l1 value contribute to the l2 error (according to the weights)?
	    l1_error = l2_delta.dot(syn1.T)
	    
	    # in what direction is the target l1?
	    # were we really sure? if so, don't change too much.
	    l1_delta = l1_error * nonlin(l1,deriv=True)

	    syn1 += l1.T.dot(l2_delta)*0.000001
	    syn0 += l0.T.dot(l1_delta)*0.000001

	#print syn0
	#print syn1
	return syn0, syn1
	

def modify(xtrain):
    chemolist=[]
    for i in range(len(xtrain)):
        for j in range(1,11):
            if xtrain[i][j] == 'YES' or xtrain[i][j]=='Yes':
                xtrain[i][j]=1
            elif xtrain[i][j]=='NO' or xtrain[i][j]=='No':
                xtrain[i][j]=-1
            elif xtrain[i][j] =='POS':
                xtrain[i][j]=1
            elif xtrain[i][j]=='NEG':
                xtrain[i][j]=-1
            elif xtrain[i][j]=='ND' or xtrain[i][j]=='NotDone':
                xtrain[i][j]=0
            elif xtrain[i][j]=='F':
                xtrain[i][j]=1
            elif xtrain[i][j]=='M':
                xtrain[i][j]=2
        if xtrain[i][11] in chemolist:
            xtrain[i][11]=chemolist.index(xtrain[i][11])
        elif xtrain[i][11] not in chemolist:
            chemolist.append(xtrain[i][11])
            xtrain[i][11]=chemolist.index(xtrain[i][11])
            
    return xtrain,chemolist

def modify_output(ytrain):
	for i in range(len(ytrain)):
		if ytrain[i][0] == 'COMPLETE_REMISSION':
			ytrain[i][0] = 0
		else:
			ytrain[i][0] = 1
	return ytrain



def filter_data():
	with open("trainingData.txt") as f:
		mydata=[[value for value in line.split()]for line in f]

	feature_names = mydata[2]
	data = mydata[2:]

	xtrain=[]
	ytrain=[]
	for i in range(len(data)):
		xtrain.append(data[i][:266])
		ytrain.append(data[i][266:])
		(xtrain_mod, chemolist) = modify(xtrain)

		ytrain_mod = modify_output(ytrain)


		patients = []
		for i in range(len(xtrain_mod)):
			patient_data = xtrain_mod[i]
			for j in range(1, len(patient_data)):
				if patient_data[j] == 'NA':
					patient_data[j] = 0.0
				else:
					patient_data[j] = float(patient_data[j])
			patients.append(patient_data[1:])
		print (patients[0])



		remission_output = []

		for i in range(len(ytrain_mod)):
			remission_output.append([ytrain_mod[i][0]])

		x_data = np.array(patients)
		# Remission data
		y_data = np.array(remission_output)
		print (x_data)
		print (y_data)
		return x_data, y_data



def evaluate(weights0, weights1, new_datax):

	#Forward propagation
	l0 = new_datax
	l1 = nonlin(np.dot(l0,weights0))
	l2 = nonlin(np.dot(l1,weights1)) 

	return l2


# feats = np.array([[0,0,1],
#             [0,1,1],
#             [1,0,1],
#             [1,1,1]])

# res = np.array([[0],
# 			[1],
# 			[1],
# 			[0]])

#neural_net(feats, res)

(features, results) = filter_data()
#print len(features) 166
#print results
(weights0, weights1) = neural_net(features[0:160], results[0:160])

guess = evaluate(weights0, weights1, features[161:len(features)])

print (guess)
print (results[161:len(results)])
