# -*- coding: utf-8 -*-
import os
import h5py as h5
import numpy as np
import math
import pandas as pd
import traceback


def getGreatCircleDistance(ra1, dec1, ra2, dec2):
    rst = 57.295779513 * math.acos(math.sin(0.017453293 * dec1) * math.sin(0.017453293 * dec2)
            + math.cos(0.017453293 * dec1) * math.cos(0.017453293 * dec2) * math.cos(0.017453293 * (math.fabs(ra1 - ra2))));
    return rst   
        
def tmatch():
    
    gwacData='/data/work/wj/crossmatch200721'
    
    dirs = os.listdir(gwacData)
    dirs.sort()
    
    allDatas = np.array([])
    for iii, fname in enumerate(dirs):

        try:
            if fname.find('Field_')==0:
                print("%d, %s"%(iii, fname))
                fullPath = "%s/%s"%(gwacData, fname)
                tds = np.loadtxt(fullPath, dtype='str', delimiter=',')
                if allDatas.shape[0]==0:
                    allDatas=tds[1:]
                else:
                    allDatas=np.concatenate((allDatas, tds[1:]), axis=0)
                
                #print(allDatas.shape)
        except Exception as e:
            print(str(e))
            tstr = traceback.format_exc()
            print(tstr)
        
        #if iii>=0:
        #    break
    
    tdf = pd.DataFrame(allDatas)
    tdfSort = tdf.sort_values([14, 15])
    
    destDir = '/data/work/wj/crossmatch200721Merge'
    if not os.path.exists(destDir):
        os.system("mkdir -p %s"%(destDir))
    savePath = '%s/gwac_%d.csv'%(destDir, tdfSort.shape[0])
    print(savePath)
    tdfSort.to_csv(savePath, index=False)
    
    ztfData='/data/work/wj/ztfData'
    fname = ''
    kname = ''
    
    mchRows = []
    for index, td in tdfSort.iterrows():
        tid = int(float(td[13]))
        tfname = td[14]
        tkname = td[15]
        
        if fname!=tfname:
            fname=tfname
            tfullpath = "%s/%s"%(ztfData, fname)
            f=h5.File(tfullpath,'r')

        if kname!=tkname:
            kname=tkname
            tdata=np.array(f[kname])
            mag = tdata[8,:]
            tdata=tdata[:,mag<18]
            
        mchRow = tdata[:,tid]
        mchRows.append(mchRow)
    
    destDir = '/data/work/wj/crossmatch200721Merge'
    GPS1_df = pd.DataFrame(mchRows)
    savePath = '%s/ztf_%d.csv'%(destDir, GPS1_df.shape[0])
    print(savePath)
    GPS1_df.to_csv(savePath, index=False)
    

if __name__ == '__main__':
    
    tmatch()