import sys
parentdir = __file__.split('/tests')[0]
sys.path.extend([parentdir+"/src", parentdir + "/src/vcsolver"]) 
print(sys.path)