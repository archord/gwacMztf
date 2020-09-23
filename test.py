# -*- coding: utf-8 -*-
import os
import h5py as h5
import numpy as np

def readFile(fullpath):
    
    print(fullpath)
    f=h5.File(fullpath,'r')
    for k in f.keys():
        if k.find('_Ind')>0:
            continue
        tdata=np.array(f[k])
        if tdata.shape[0]==54:
            td=tdata[:,(tdata[48,:]==2) or (tdata[48,:]==3)]
        else:
            print("%s column error: column number is %d"%(fullpath, tdata.shape[0]))
        
    f.close()

def readAll(tpath):
    
    rstData = np.array([])
    
    files = os.listdir(tpath)
    for fname in files:
        print("read file %s"%(fname))
        fullpath = "%s/%s"%(tpath, fname)
        f=h5.File(fullpath,'r')
        for k in f.keys():
            if k.find('_Ind')>0:
                continue
            tdata=np.array(f[k])
            if tdata.shape[0]==54:
                td=tdata[:,((tdata[47,:]>=2) & (tdata[47,:]<=3))]
                if td.shape[1]>0:
                    #print(td[47])
                    if len(rstData.shape)==2:
                        rstData = np.concatenate((rstData, td), axis=1)
                    else:
                        rstData = td
            else:
                print("%s column error: column number is %d"%(fullpath, tdata.shape[0]))
            
        f.close()
    
    if len(rstData.shape)==2:
        print("total get %d records, and save to allOut.txt"%(rstData.shape[1]))
        np.savetxt("allOut.txt", rstData.transpose())

if __name__ == '__main__':
    
    fname='ZTF/ztfDR1var'
    readAll(fname)