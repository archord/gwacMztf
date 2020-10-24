# -*- coding: utf-8 -*-
import os
import h5py as h5
import numpy as np
import math
import pandas as pd
from astropy_healpix import HEALPix
from astropy import units as u
import traceback


def getGreatCircleDistance(ra1, dec1, ra2, dec2):
    rst = 57.295779513 * math.acos(math.sin(0.017453293 * dec1) * math.sin(0.017453293 * dec2)
            + math.cos(0.017453293 * dec1) * math.cos(0.017453293 * dec2) * math.cos(0.017453293 * (math.fabs(ra1 - ra2))));
    return rst   
 

def addIdx(tidx, hp, ra, dec, origIdx):
    
    hpix = hp.lonlat_to_healpix(ra * u.deg, dec * u.deg)

    j=0
    for tpix in hpix:
        tidx[tpix].append((ra[j],dec[j],j, origIdx))
        j=j+1


def buildIdx(hp, tpath):
    
    arc2deg = 180/np.pi
    tidx = []
    origIdx=[]
    
    for i in range(hp.npix):
        tidx.append([])
    
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
                origIdx.append((fname, k))
                ra = tdata[0,:]
                dec = tdata[1,:]
                mag = tdata[8,:]
                ra=ra[mag<18]*arc2deg
                dec=dec[mag<18]*arc2deg
                addIdx(tidx, hp, ra, dec, len(origIdx)-1)
            else:
                print("%s column error: column number is %d"%(fullpath, tdata.shape[0]))
            
        f.close()
        #break
    
    return tidx, origIdx

def domatch1(hp, tidxs, origIdx, tpath, fname):
    
    maxDist = 30.0/60/60
    fullPath = "%s/%s"%(tpath, fname)
        
    print(fname)
    
    tdata00 = np.loadtxt(fullPath, dtype='str')
    
    tmatched = []
    for td in tdata00:
        
        try:
            #print(td)
            ra1=float(td[1])
            dec1=float(td[2])
            hpixs = hp.cone_search_lonlat(ra1 * u.deg, dec1 * u.deg, radius=maxDist * u.deg)
            #print(hpixs)
            
            for ti in hpixs:
                tposs = tidxs[ti]
                for tpos in tposs:
                    ra2=tpos[0]
                    dec2=tpos[1]
                    #print(tpos)
                    tdist = getGreatCircleDistance(ra1, dec1, ra2, dec2)
                    if tdist<=maxDist:
                        rst0 = np.append(td, tpos[0:3])
                        rst0 = np.append(rst0, origIdx[tpos[3]])
                        tmatched.append(rst0) #ra[j],dec[j],j, fname, k
                        
            #break
        except Exception as e:
            print("domatch error")
            print(str(e))
            tstr = traceback.format_exc()
            print(tstr)
    
    print("local tmatched %d"%(len(tmatched)))
    if len(tmatched)>0:
        destDir = '/data/work/wj/crossmatch200925'
        if not os.path.exists(destDir):
            os.system("mkdir -p %s"%(destDir))
        savePath = '%s/%s_%04d.csv'%(destDir,fname.split('.')[0],len(tmatched))
        print(savePath)
        GPS1_df = pd.DataFrame(tmatched)
        GPS1_df.to_csv(savePath, index=False)
        
def tmatch():
    
    
    ztfData='/data/work/wj/ztfData'
    gwacData='/data/work/wj/dat'
    gwacDataFiles='/data/work/wj/filteredFileList.txt'
    tidx00='/data/work/wj/ztfIdx2.npz'
    
    hp = HEALPix(nside=256)
    print(hp.npix)
    if os.path.exists(tidx00):
        print("load ztfIdx from cache")
        aa = np.load(tidx00, allow_pickle=True)
        tidxs=aa['tidxs']
        origIdx=aa['origIdx']
    else:
        tidxs, origIdx = buildIdx(hp, ztfData)
        print("save ztfIdx to cache")
        np.savez_compressed(tidx00, tidxs=tidxs, origIdx=origIdx)
      
    dirs = np.loadtxt(gwacDataFiles, dtype='str')
    
    iii=0
    for tdir in dirs:

        try:
            if tdir.find('Field_')==0:
                print(iii)
                domatch1(hp, tidxs, origIdx, gwacData, tdir)

        except Exception as e:
            print(str(e))
            tstr = traceback.format_exc()
            print(tstr)
        
        #if iii>=0:
        #    break
        iii = iii + 1

if __name__ == '__main__':
    
    tmatch()