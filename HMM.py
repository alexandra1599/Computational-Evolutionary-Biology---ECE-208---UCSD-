from math import *

INF = float("inf")
def log_sum_exp(numlist):
    try:
        minx = min([x for x in numlist if x != -INF])
    except:
        return -INF    
    s = sum(exp(x-minx) for x in numlist)
    return minx + log(s) if s > 0 else -INF

class HMM:
    def __init__(self):
        #the alphabet of the HMM model. This is a list of characters.
        self.alphabet = []
        # emission probability of the I states, one dictionary for each I in the model
        self.eI = [] 
        # emission probability of the M states, one dictionary for each M in the model
        # the first M state is called 'B'; it never emits anything so the associated dictionary is always empty
        self.eM = [{}] 
        # transition probability, one dictionary for each set of states (D,M,I) in the model
        self.t = [] 
        # number of matching states in the model, excluding B and E
        self.nstate = 0
    
    def load(self,hmmfile):
    # only load the first model in the given hmmfile if there are two or more models
        with open(hmmfile,'r') as fin:
            for line in fin:
                stream = line.strip().split()
                if stream[0] == "LENG":
                    self.nstate = int(stream[1])
                if stream[0] == "HMM": 
                    # read alphabet
                    self.alphabet = stream[1:]
                    # read transition order
                    stream = fin.readline().strip().split()
                    trans_order = [(y[0]+y[3]).upper() for y in stream]
                    # read the next line, if it is the COMPO line then ignore it and read one more line
                    stream = fin.readline().strip().split()
                    if stream[0] == "COMPO":
                        stream = fin.readline().strip().split()
                    # now the stream should be at the I0 state; read the emission of the I0 state 
                    e = {}
                    for (x,y) in zip(self.alphabet,stream):
                        e[x] = -float(y)
                    self.eI.append(e)    
                    # now the stream should be at the B state; read in the transition probability
                    stream = fin.readline().strip().split()
                    tB = {'MM':-INF,'MD':-INF,'MI':-INF,'IM':-INF,'II':-INF,'ID':-INF,'DM':-INF,'DI':-INF,'DD':-INF}
                    for x,y in zip(trans_order,stream):
                        tB[x] = -INF if y == '*' else -float(y)
                    self.t.append(tB)
                    break    
            
            for i in range(1,self.nstate+1):
                
                stream = fin.readline().strip().split() # this one is the emission of the M state
                if float(stream[0]) != i:
                    print("Warning: incosistent state indexing in hmm file; expecting state "+ str(i) + "; getting state " + stream[0])
                a = {}
                for x,y in zip(self.alphabet,stream[1:]):
                    a[x] = -INF if y == "*" else -float(y) 
                self.eM.append(a)    
                
                stream = fin.readline().strip().split() 
                emissionp = {}
                for x,y in zip(self.alphabet,stream):
                    emissionp[x] = -INF if y == "*" else -float(y) 
                self.eI.append(emissionp)                    
                
                stream = fin.readline().strip().split() # this one is the transition prob
                transp = {'MM':-INF,'MD':-INF,'MI':-INF,'IM':-INF,'II':-INF,'ID':-INF,'DM':-INF,'DI':-INF,'DD':-INF}
                for x,y in zip(trans_order,stream):
                    transp[x] = -INF if y == '*' else -float(y)
                self.t.append(transp)
    
    def compute_llh(self,query):
        # compute likelihood of an aligned query to the HMM
        # return -inf if the query is not properly aligned
        j = 0
        prev_state = 'M'
        llh = 0

        for c in query:
            if c == '-': # deletion
                curr_state = 'D'
            elif c >= 'a' and c <= 'z': # insertion
                curr_state = 'I'
            elif c >= 'A' and c <= 'Z': # match
                curr_state = 'M'
            else: # invalid symbol
                return -INF
            
            transition = prev_state + curr_state

            # penalize the transition
            if transition in self.t[j]:
                llh += self.t[j][transition]
            else: # invalid transition
                return -INF 
            
            # update the state index
            j += (curr_state != 'I') # move to the next state unless this is a move towards an 'I'
            if j > self.nstate: # reach end of the HMM chain but not end of the sequence
                return -INF

            # penalize the emission
            if curr_state == 'M':
                llh += self.eM[j][c]
            elif curr_state == 'I': 
                llh += self.eI[j][c.upper()]
            
            # update state
            prev_state = curr_state

        if j != self.nstate: # does not reach the E state at the end
            return float("-inf")
        
        # towards the 'E' state
        transition = prev_state + 'M'
        llh += self.t[j][prev_state+'M']    
        return llh      
  
    def Viterbi(self,query):
        # implement the Viterbi algorithm to find the most probable path of a query sequence and its likelihood
        # return the log-likelihood and the path. If the likelihood is 0, return "$" as the path 
        n = len(query)
        m = self.nstate
        
        VM_alex = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        VI_alex = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        VD_alex = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
                
        Vscore = -INF
        aln = "$"
        # initialization
        VM_alex[0][1] = self.t[0]['MM'] + self.eM[1][query[0]]
        VD_alex[0][1] = self.t[0]['MD']                        
        VI_alex[0][0] = self.t[0]['MI'] + self.eI[0][query[0]] 
        
        # YOUR CODE HERE
        
        pathVM_alex = [[list() for j in range(m+2)] for i in range(n+1)]
        pathVI_alex = [[list() for j in range(m+2)] for i in range(n+1)]
        pathVD_alex = [[list() for j in range(m+2)] for i in range(n+1)]
        
        #initialize VD first
        for j in range(2,m+2):
            VD_alex[0][j] = VD_alex[0][j-1]+self.t[j-1]['DD']
            pathVD_alex[0][j].append(2)
        
        #initialize VM
        for j in range(2,m+1):
            VM_alex[0][j] = self.eM[j][query[0]] + VD_alex[0][j-1]+self.t[j-1]['DM']
            pathVM_alex[0][j].append(2)
        
        #intitialize VI
        #first row
        for j in range(1,m+1):
            VI_alex[0][j] = self.eI[j][query[0]] + VD_alex[0][j] + self.t[j]['DI']
            pathVI_alex[0][j].append(2)
            
        for i in range(1,n):
            VI_alex[i][0] = self.eI[0][query[i]] + VI_alex[i-1][0]+ self.t[0]['II']
            pathVI_alex[i][0].append(1)
        
        #j is  model state; j is horozontal
        #i is the point in the subsequence that we are at; i is vertical
        a=0
        for j in range(1,m+1):
            for i in range(1,n):
                
                maxVM_alex = [VM_alex[i-1][j-1] + self.t[j-1]['MM'], VI_alex[i-1][j-1] + self.t[j-1]['IM'], VD_alex[i][j-1] + self.t[j-1]['DM'] ]
                VM_alex[i][j] = self.eM[j][query[i]]+ max(maxVM_alex)
                indic = maxVM_alex.index(max(maxVM_alex))
                if indic == 0: 
                    tempVM_alex = pathVM_alex[i-1][j-1][:]
                if indic == 1: 
                    tempVM_alex = pathVI_alex[i-1][j-1][:]
                if indic == 2:
                    tempVM_alex = pathVD_alex[i][j-1][:]
                tempVM_alex.append(indic)
                
                
                maxVI_alex = [VM_alex[i-1][j] + self.t[j]['MI'], VI_alex[i-1][j] + self.t[j]['II'], VD_alex[i][j] + self.t[j]['DI']]
                VI_alex[i][j] = self.eI[j][query[i]]+ max(maxVI_alex)
                indic = maxVI_alex.index(max(maxVI_alex))
                if indic == 0: 
                    tempVI_alex = pathVM_alex[i-1][j][:]
                if indic == 1: 
                    tempVI_alex = pathVI_alex[i-1][j][:]
                if indic == 2:
                    tempVI_alex = pathVD_alex[i][j][:]
                tempVI_alex.append(indic)
                
                maxVD_alex = [VM_alex[i-1][j-1] + self.t[j-1]['MD'], VI_alex[i-1][j-1] + self.t[j-1]['ID'], VD_alex[i][j-1] + self.t[j-1]['DD']]
                VD_alex[i][j] = max(maxVD_alex)
                indic = maxVD_alex.index(max(maxVD_alex))
                if indic == 0: 
                    tempVD_alex = pathVM_alex[i-1][j-1][:]
                if indic == 1: 
                    tempVD_alex = pathVI_alex[i-1][j-1][:]
                if indic == 2:
                    tempVD_alex = pathVD_alex[i][j-1][:]
                tempVD_alex.append(indic)
                
                pathVM_alex[i][j] = tempVM_alex[:]
                pathVI_alex[i][j] = tempVI_alex[:]
                pathVD_alex[i][j] = tempVD_alex[:]
                
                
                
        #last row of VD
        i = n
        for j in range(0, m+1):
            maxVD_alex = [VM_alex[i-1][j-1] + self.t[j-1]['MD'], VI_alex[i-1][j-1] + self.t[j-1]['ID'], VD_alex[i][j-1] + self.t[j-1]['DD']]
            VD_alex[i][j] = max(maxVD_alex)
            indic = maxVD_alex.index(max(maxVD_alex))
            if indic == 0: 
                pathVD_alex[i][j] = pathVM_alex[i-1][j-1][:]
            if indic == 1: 
                pathVD_alex[i][j] = pathVI_alex[i-1][j-1][:]
            if indic == 2:
                pathVD_alex[i][j] = pathVD_alex[i][j-1][:]
            pathVD_alex[i][j].append(indic)
        
        i=n
        j=m+1
        maxVM_alex = [VM_alex[i-1][j-1] + self.t[j-1]['MM'], VI_alex[i-1][j-1] + self.t[j-1]['IM'], VD_alex[i][j-1] + self.t[j-1]['DM'] ]
        VM_alex[i][j] = max(maxVM_alex)
        
        indic = maxVM_alex.index(max(maxVM_alex))
        if indic == 0: 
            pathVM_alex[i][j] = pathVM_alex[i-1][j-1][:]
        if indic == 1: 
            pathVM_alex[i][j] = pathVI_alex[i-1][j-1][:]
        if indic == 2:
            pathVM_alex[i][j] = pathVD_alex[i][j-1][:]
        pathVM_alex[i][j].append(indic)
        
        if VM_alex[i][j] != -INF:
            aln = ""
            l = 0
            for k in pathVM_alex[i][j]:
                if k == 0:
                    aln = aln + (query[l].upper())
                    l = l+1
                elif k == 1:
                    aln = aln + (query[l].lower())
                    l = l+1
                elif k == 2:
                    aln = aln + '-'
        
        Vscore = VM_alex[n][m+1]
        return Vscore,aln        


    def Forward(self,query):
        # implement the Forward algorithm to compute the marginal probability of a query sequence
        n = len(query)
        m = self.nstate
        
        FM_alex = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        FI_alex = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        FD_alex = [[float("-inf") for j in range(m+2)] for i in range(n+1)]
        
        # initialization
        FM_alex[0][1] = self.t[0]['MM'] + self.eM[1][query[0]]  
        FD_alex[0][1] = self.t[0]['MD']           
        FI_alex[0][0] = self.t[0]['MI'] + self.eI[0][query[0]] 

        # YOUR CODE HERE
        
        #initialize FD first
        for j in range(2,m+2):
            FD_alex[0][j] = FD_alex[0][j-1]+self.t[j-1]['DD']
            
        
        #initialize FM
        for j in range(2,m+1):
            FM_alex[0][j] = self.eM[j][query[0]] + FD_alex[0][j-1]+self.t[j-1]['DM']
        
        
        #intitialize FI
        #first row
        for j in range(1,m+1):
            FI_alex[0][j] = self.eI[j][query[0]] + FD_alex[0][j] + self.t[j]['DI']
           
            
        for i in range(1,n):
            FI_alex[i][0] = self.eI[0][query[i]] + FI_alex[i-1][0]+ self.t[0]['II']
        
        
               
        #j is  model state; j is horozontal
        #i is the point in the subsequence that we are at; i is vertical
        a=0
        for j in range(1,m+1):
            for i in range(1,n):
                
                maxFM_alex = [FM_alex[i-1][j-1] + self.t[j-1]['MM'], FI_alex[i-1][j-1] + self.t[j-1]['IM'], FD_alex[i][j-1] + self.t[j-1]['DM'] ]
                FM_alex[i][j] = self.eM[j][query[i]]+ log_sum_exp(maxFM_alex)
                
                
                maxFI_alex = [FM_alex[i-1][j] + self.t[j]['MI'], FI_alex[i-1][j] + self.t[j]['II'], FD_alex[i][j] + self.t[j]['DI']]
                FI_alex[i][j] = self.eI[j][query[i]]+ log_sum_exp(maxFI_alex)
               
                
                maxFD_alex = [FM_alex[i-1][j-1] + self.t[j-1]['MD'], FI_alex[i-1][j-1] + self.t[j-1]['ID'], FD_alex[i][j-1] + self.t[j-1]['DD']]
                FD_alex[i][j] = log_sum_exp(maxFD_alex)
                
                #last row of FD
        i = n
        for j in range(0, m+1):
            maxFD_alex = [FM_alex[i-1][j-1] + self.t[j-1]['MD'], FI_alex[i-1][j-1] + self.t[j-1]['ID'], FD_alex[i][j-1] + self.t[j-1]['DD']]
            FD_alex[i][j] = log_sum_exp(maxFD_alex)
           
        
        i=n
        j=m+1
        maxFM_alex = [FM_alex[i-1][j-1] + self.t[j-1]['MM'], FI_alex[i-1][j-1] + self.t[j-1]['IM'], FD_alex[i][j-1] + self.t[j-1]['DM'] ]
        FM_alex[i][j] = log_sum_exp(maxFM_alex)
       
        Fscore = FM_alex[n][m+1]        
        return Fscore
