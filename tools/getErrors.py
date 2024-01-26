import os

nFile, nError = 0, 0  
for root, directories, filenames in os.walk('.'):
    for filename in filenames:
        if filename[-4:] == '.err' :
            fName = os.path.join(root,filename)
            nFile += 1 
            #print("file={0:s}".format(fName))
            iLast = -1 
            for i, line in enumerate(open(fName,'r').readlines()) :
                lineLower = line.lower() 
                if '====' in lineLower : continue
                if 'warning' in lineLower : continue
                if 'info in' in lineLower : continue
                if 'mass expected' in lineLower : continue
                if 'rm: no match.' in lineLower : continue
                if 'cloning into' in lineLower : continue
                if 'in member function' in lineLower : continue
                if '%msg' in lineLower : continue
                if 'nanoaodtools' in lineLower : continue
                if 'double val' in lineLower : continue
                if '^~~' in line : continue
                if 'checking out files' in lineLower : continue 
                if i > iLast+1 : 
                    print('***\nError in file{0:s}'.format(fName)) 
                    nError += 1 
                print("    [{0:4d}]={1:s}".format(i,line.strip()))
                iLast = i
        
print("Checked {0:d} files.  {1:d} errors found.".format(nFile,nError))

 
