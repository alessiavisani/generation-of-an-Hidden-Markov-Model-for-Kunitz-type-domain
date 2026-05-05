#!/Users/aleonmac/miniconda3/envs/lab_bioinformatics_env/bin/python
import sys
import numpy as np

def get_preds(fname): #extracts the information from the prediciton file 
    preds=[]
    fh=open(fname)

    for line in fh:
        v=line.rstrip().split() #remove the last element and obtain the different positions 
        preds.append([v[0], float(v[1]), int(v[2])]) # e value in the column 2 and real value in the column 3 
    return preds



def get_confusion_matrix(preds, threshold=0.001):
    # i need to go line by line and assignening the cases in the correct position
    confusion_matrix=np.zeros((2,2)) 
    n=len(preds)
    for k in range(n):
        j=0
        i=int(preds[k][2]) #the ones considered to be the truth
        
        if len(preds[k]) < 3: continue
        if float(preds[k][1]) <= threshold: 
            j=1 #if the preds in position i 1 è minore del threshold j=1 # e value in position one i will get the prediction
        else: j=0
        confusion_matrix[i,j]=confusion_matrix[i,j]+1
    return confusion_matrix

def get_accuracy(confusion_matrix):
    total=np.sum(confusion_matrix)
    return(confusion_matrix[0,0]+confusion_matrix[1,1])/total


def get_mcc(cm): #matthew correlation coeffcient
    true_pos=cm[1,1]
    true_neg=cm[0,0]
    false_neg=cm[1,0]
    false_pos=cm[0,1]

    numerator = (true_pos *true_neg)-(false_pos*false_neg)
    d=(true_pos + false_pos) * (true_pos + false_neg) * (true_neg +false_pos) * (true_neg +false_neg)
    mcc=numerator/np.sqrt(d)
    return mcc



if __name__ == "__main__":
    fname=sys.argv[1]
    threshold=float(sys.argv[2])
    preds=get_preds(fname)
    

    cm= get_confusion_matrix(preds,threshold)
    

    q2=get_accuracy(cm)
    

    mcc=get_mcc(cm)
    print("TH:",threshold , "Q2:",q2,"MCC",mcc)
    
    #python performance.py blast_preds.txt 0.001

    #for i in `seq 1 15`; do python performance.py kunitz_set_1.txt 1e-$i; done > kunitz_set_1.performance
