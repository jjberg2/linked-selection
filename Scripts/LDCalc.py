# sys.argv[1] gives the ms output file prefix
# sys.argv[2] gives the number of discrete loci along the chromosome
# sys.argv[3] gives the number of ms output files to loop over
import sys,pdb,numpy as np,bisect;
def get_pos_line(file):
    pos = file.readline()
    pos = pos.strip()
    pos = pos.split(":")[1]
    pos = pos.split(" ")[1:]
    return pos

num_loc = int ( sys.argv[2] )
my_bins = [ float(i/num_loc) for i in range ( 0,num_loc+1) ]
my_lens = []
for filenum in range(1,int(sys.argv[3])): 
    with open (sys.argv[1] + str(filenum) , "r") as ldf:
        pos=[]
        for _ in range(5):
            next(ldf)
        pos = get_pos_line(ldf)
        my_lens.append(str(len(pos))) 
               
 
my_max = max(list(map(int, my_lens)))
mycov = np.zeros([my_max,my_max,int(sys.argv[3])])
for filenum in range(1,int(sys.argv[3])): 
    with open (sys.argv[1] + str(filenum) , "r") as ldf:
            pos=[]
            for _ in range(5):
                next(ldf)
            pos = get_pos_line(ldf)
            chromread=ldf.readlines()
            mychrom=[]
            for idx,line in enumerate(chromread):
                line=line.strip()
                mychrom.append(list(line))
            x=np.array(mychrom).T
            mycov[0:len(pos),0:len(pos),filenum-1] = np.cov(x)

pdb.set_trace()
for idx,i in enumerate(my_bins):
    idx_list.append([np.searchsorted(pos,i,side="right"),np.searchsorted(pos,i+float(1/num_loc),side="left")])
            
            
            
