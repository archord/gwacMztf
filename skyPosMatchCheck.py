# -*- coding: utf-8 -*-
import os
import h5py as h5
import numpy as np
import math
import pandas as pd
from astropy_healpix import HEALPix
from astropy import units as u
import traceback
import matplotlib.pyplot as plt


def getGreatCircleDistance(ra1, dec1, ra2, dec2):
    rst = 57.295779513 * math.acos(math.sin(0.017453293 * dec1) * math.sin(0.017453293 * dec2)
            + math.cos(0.017453293 * dec1) * math.cos(0.017453293 * dec2) * math.cos(0.017453293 * (math.fabs(ra1 - ra2))));
    return rst   


def buildIdx(hp, tdata):
    
    tidx = []
    
    for i in range(hp.npix):
        tidx.append([])
    
    ra = tdata[:,1].astype(np.float32)
    dec = tdata[:,2].astype(np.float32)    
    hpix = hp.lonlat_to_healpix(ra * u.deg, dec * u.deg)

    j=0
    for tpix in hpix:
        tidx[tpix].append((ra[j],dec[j],j))
        j=j+1
    
    return tidx


def crossMatchStatistics():
    
    
    gwacData='/data/work/wj/dat'      
    dirs = os.listdir(gwacData)
    dirs.sort()
    
    tNames = []
    tDatas = []
    obsNum = 0
    for iii, tdir in enumerate(dirs):

        try:
            if tdir.find('21450595')>0:
                obsNum = obsNum + 1
                fullPath = "%s/%s"%(gwacData, tdir)
                tdata00 = np.loadtxt(fullPath, dtype='str')
                if tdata00.shape[0]>10000:
                    tNames.append(tdir)
                    tDatas.append(tdata00)
                    print("fName is %s, num is %d"%(tdir, tdata00.shape[0]))

        except Exception as e:
            print(str(e))
            tstr = traceback.format_exc()
            print(tstr)

    print("total file is %d, large 10000 is %s"%(obsNum, len(tNames)))
    
    maxDist = 30/60/60
    hp = HEALPix(nside=256)
    
    for j in range(len(tNames)):
        print("\n\n*****************")
        print("start match %d file: %s, has %d stars"%(j+1, tNames[j], tDatas[j].shape[0]))
        templateData = tDatas[j]
        tidxs = buildIdx(hp, templateData)
        mchList = []
        for i in range(len(tNames)):
            tdata = tDatas[i]
            if i!=j:
                tnum00 =0
                for td in tdata:
                    ra1=float(td[1])
                    dec1=float(td[2])
                    hpixs = hp.cone_search_lonlat(ra1 * u.deg, dec1 * u.deg, radius=maxDist * u.deg)
                    
                    isMatch = False
                    for ti in hpixs:
                        tposs = tidxs[ti]
                        for tpos in tposs:
                            ra2=tpos[0]
                            dec2=tpos[1]
                            try:
                                tdist = getGreatCircleDistance(ra1, dec1, ra2, dec2)
                                if tdist<=maxDist:
                                    isMatch = True
                                    tnum00 = tnum00+1
                                    if tnum00<5:
                                        print("ra1=%f,dec1=%f,ra2=%f,dec2=%f,dist=%f"%(ra1, dec1, ra2, dec2, tdist))
                                    break
                            
                            except Exception as e:
                                print("domatch error")
                        if isMatch:
                            break
                
                print("%d, fName is %s, has %d stars, match num is %d"%(i+1, tNames[i], tdata.shape[0], tnum00))
                if tnum00>1000:
                    mchList.append([i, tnum00])
        print(mchList)

#sky-ccd sky-coornidate statistic and hist plot
def tmatch():
    
    
    gwacData='/data/work/wj/dat'      
    dirs = os.listdir(gwacData)
    dirs.sort()
    
    tNames = []
    tDatas = []
    for iii, tdir in enumerate(dirs):

        try:
            if tdir.find('21450595')>0:
                fullPath = "%s/%s"%(gwacData, tdir)
                tdata00 = np.loadtxt(fullPath, dtype='str')
                if tdata00.shape[0]>10000:
                    tNames.append(tdir)
                    tDatas.append(tdata00)
                    print("fName is %s, num is %d"%(tdir, tdata00.shape[0]))

        except Exception as e:
            print(str(e))
            tstr = traceback.format_exc()
            print(tstr)

    skyMapName = {}
    skyMapData = {}
    for i in range(len(tNames)):
        tkey = tNames[i][8:18]
        if tkey not in skyMapName:
            skyMapName[tkey] = []
            skyMapData[tkey] = []
        skyMapName[tkey].append(tNames[i])
        skyMapData[tkey].append(tDatas[i])

    allSkyDist = []
    for tkey in skyMapName:
        
        print("\n\n*****************")
        print("match sky %s"%(tkey))
        
        tNames = skyMapName[tkey]
        tDatas = skyMapData[tkey]
        
        maxNumIdx = -1
        maxNum = -1
        for i in range(len(tNames)):
            if tDatas[i].shape[0]>maxNum:
                maxNum = tDatas[i].shape[0]
                maxNumIdx = i
                
        print("maxNum is %d, fName is %s"%(maxNum, tNames[maxNumIdx]))
        
        maxDist = 30.0/60/60
        templateData = tDatas[maxNumIdx]
        hp = HEALPix(nside=256)
        tidxs = buildIdx(hp, templateData)
        
        allDists1 = []
        for i in range(len(tNames)):
            tdata = tDatas[i]
            allDists2 = []
            if i!=maxNumIdx and tdata.shape[0]>10000:
                for td in tdata:
                    #print(td)
                    ra1=float(td[1])
                    dec1=float(td[2])
                    hpixs = hp.cone_search_lonlat(ra1 * u.deg, dec1 * u.deg, radius=maxDist * u.deg)
                    #print(hpixs)
                    
                    minDist = 100000
                    for ti in hpixs:
                        tposs = tidxs[ti] #tidx[tpix].append((ra[j],dec[j],j))
                        for tpos in tposs:
                            ra2=tpos[0]
                            dec2=tpos[1]
                            
                            try:
                                tdist = getGreatCircleDistance(ra1, dec1, ra2, dec2)
                                if tdist<=minDist:
                                    minDist = tdist
                            
                            except Exception as e:
                                print("domatch error")
                                #print(str(e))
                                #tstr = traceback.format_exc()
                                #print(tstr)
                
                    if minDist <= maxDist:
                        allDists2.append(minDist)
                print("fName is %s, has %d stars, match num is %d"%(tNames[i], tdata.shape[0], len(allDists2)))
            allDists1.append(allDists2)
        allSkyDist.append(allDists1)
        
        for j, tdists in enumerate(allSkyDist):
    
            maxNum = 0
            maxIdx = -1
            for i, tds in enumerate(tdists):
                if len(tds)>maxNum:
                    maxNum = len(tds)
                    maxIdx = i
            
            if maxIdx>-1:
                print((maxIdx, maxNum))
                x = np.array(tdists[maxIdx])*3600
                bins=np.arange(0,30,1)
                plt.hist(x,bins,color='fuchsia',alpha=0.5)
                plt.xlabel('distance')
                plt.ylabel('count')
                plt.grid()
                plt.xlim(0,30)
                plt.show()

def templateChoseByStarMchNum():
    
    
    gwacData='/data/work/wj/dat'      
    dirs = os.listdir(gwacData)
    dirs.sort()
    
    tNames = []
    tDatas = []
    for iii, tdir in enumerate(dirs):

        try:
            fullPath = "%s/%s"%(gwacData, tdir)
            tdata00 = np.loadtxt(fullPath, dtype='str')
            if tdir.find('Field_')==0 and tdata00.shape[0]>10000:
                tNames.append(tdir)
                tDatas.append(tdata00.shape[0])
                print("fName is %s, num is %d"%(tdir, tdata00.shape[0]))

        except Exception as e:
            print(str(e))
            tstr = traceback.format_exc()
            print(tstr)

    skyMapName = {}
    skyMapData = {}
    for i in range(len(tNames)):
        tkey = tNames[i][8:18]
        if tkey not in skyMapName:
            skyMapName[tkey] = []
            skyMapData[tkey] = []
        skyMapName[tkey].append(tNames[i])
        skyMapData[tkey].append(tDatas[i])

    allSkyDist = []
    for tkey in skyMapName:
        
        print("\n\n*****************")
        print("match sky %s"%(tkey))
        
        tNames = skyMapName[tkey]
        tDatas = skyMapData[tkey]
        
        maxNumIdx = -1
        maxNum = -1
        for i in range(len(tNames)):
            if tDatas[i].shape[0]>maxNum:
                maxNum = tDatas[i].shape[0]
                maxNumIdx = i
                
        print("maxNum is %d, fName is %s"%(maxNum, tNames[maxNumIdx]))
        
        maxDist = 30.0/60/60
        templateData = tDatas[maxNumIdx]
        hp = HEALPix(nside=256)
        tidxs = buildIdx(hp, templateData)
        
        allDists1 = []
        for i in range(len(tNames)):
            tdata = tDatas[i]
            allDists2 = []
            if i!=maxNumIdx and tdata.shape[0]>10000:
                for td in tdata:
                    #print(td)
                    ra1=float(td[1])
                    dec1=float(td[2])
                    hpixs = hp.cone_search_lonlat(ra1 * u.deg, dec1 * u.deg, radius=maxDist * u.deg)
                    #print(hpixs)
                    
                    minDist = 100000
                    for ti in hpixs:
                        tposs = tidxs[ti] #tidx[tpix].append((ra[j],dec[j],j))
                        for tpos in tposs:
                            ra2=tpos[0]
                            dec2=tpos[1]
                            
                            try:
                                tdist = getGreatCircleDistance(ra1, dec1, ra2, dec2)
                                if tdist<=minDist:
                                    minDist = tdist
                            
                            except Exception as e:
                                print("domatch error")
                                #print(str(e))
                                #tstr = traceback.format_exc()
                                #print(tstr)
                
                    if minDist <= maxDist:
                        allDists2.append(minDist)
                print("fName is %s, has %d stars, match num is %d"%(tNames[i], tdata.shape[0], len(allDists2)))
            allDists1.append(allDists2)
        allSkyDist.append(allDists1)
        
        for j, tdists in enumerate(allSkyDist):
    
            maxNum = 0
            maxIdx = -1
            for i, tds in enumerate(tdists):
                if len(tds)>maxNum:
                    maxNum = len(tds)
                    maxIdx = i
            
            if maxIdx>-1:
                print((maxIdx, maxNum))
                x = np.array(tdists[maxIdx])*3600
                bins=np.arange(0,30,1)
                plt.hist(x,bins,color='fuchsia',alpha=0.5)
                plt.xlabel('distance')
                plt.ylabel('count')
                plt.grid()
                plt.xlim(0,30)
                plt.show()
                
if __name__ == '__main__':
    
    tmatch()