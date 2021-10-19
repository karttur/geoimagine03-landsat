'''
Created on 7 Oct 2018

@author: thomasgumbricht
'''

#import geoimagine.gis.mj_gis_v80 as mj_gis
import zipper.explode as zipper
import support.karttur_dt as mj_dt
#from geoimagine.kartturmain import RasterProcess
#import geoimagine.sentinel.gml_transform as gml_transform
#from geoimagine.kartturmain import Composition, LayerCommon, RasterLayer
from support import ConvertLandsatScenesToStr,EarthSunDist
from mask import SingleBandMasking
from gdalutilities import GDALstuff
from support import Masker, LandsatConfidence, LandsatMasker

#from sentinelsat.sentinel import SentinelAPI
import os
#import xml.etree.ElementTree as ET
from shutil import rmtree, move, copyfileobj
import urllib.request
import shutil
import subprocess
import landsatxplore
import numpy as np
import math
import gc
from geoimagine.landsat.atmcorr import AtmCorr


class LandsatComposition:
    '''
    class for landsat compositions - duplicate of SentinelComposition
    '''
    def __init__(self, compD):  
        for key in compD:
            if '_' in compD[key]:
                exitstr = 'the "%s" parameter can not contain underscore (_): %s ' %(key, compD[key])
                exit(exitstr) 
            setattr(self, key, compD[key])
        if not hasattr(self, 'folder'):
            exitstr = 'All landsat compositions must contain a folder'
            exit(exitstr)
    
class LandsatScene(LayerCommon):
    '''Class for landsat scenes - duplicate of SentinelTile'''
    def __init__(self, composition, locusD, datumD, filepath, FN): 
        """The constructor expects an instance of the composition class."""
        LayerCommon.__init__(self)

        self.comp = composition
        
        self.locus = locusD['locus']
        self.locuspath = locusD['path']

        self.path = filepath
        self.FN = FN

        self.datum = lambda: None
        for key, value in datumD.items():
            setattr(self.datum, key, value)
        if self.datum.acqdate:
            self._SetDOY()
            self._SetAcqdateDOY()
        self._SetPath()
        
    def _SetPath(self):
        """Sets the complete path to sentinel tiles"""
        
        self.FP = os.path.join('/Volumes',self.path.volume, self.comp.system, self.comp.source, self.comp.division, self.comp.folder, self.locuspath, self.datum.acqdatestr)
        self.FPN = os.path.join(self.FP,self.FN)
        if ' ' in self.FPN:
            exitstr = 'EXITING landsat FPN contains space %s' %(self.FPN)
            exit(exitstr)
            
    def _GetMetaDataFromProduct(self,lsatprodid):
        '''Get metadata from landsat product
        '''
        self.lsatprodid = lsatprodid
        
        #splite the product id to get all the data
        sceneCompL = lsatprodid.split('_')
        self.source = sceneCompL[0]
        self.product = sceneCompL[1]
        
        self.acqdatestr = sceneCompL[3] 
        self.acqdate = mj_dt.yyyymmddDate(self.acqdatestr)
        self.acqdoy = mj_dt.DateToDOY(self.acqdate)
        
        self.wrspath = int(sceneCompL[2][0:3])
        self.wrsrow = int(sceneCompL[2][3:6])

        self.collnr = int(sceneCompL[5])
        self.collcat = sceneCompL[6]
        
    def _SetMetaDataFromDict(self, metaD):
        for key in metaD:
            setattr(self, key, metaD[key])

    def _InsertSceneMeta(self):
        #Get all the meta for this granule
        '''
        paramL = ['product','']
            metaRec = self.session._GetGranuleMeta(granuleid)

                paramD = ['source','product', 'proclevel', 'cloudcover', 'sensopmode', 's2datatakeid', 'procbase', 'platformid', 'platformname', 'instrument']
                self.metaquery = dict(zip(paramD,metaRec))
            else:
                exit('No metadata for granule in _ExplodeSentinel2Granules')
                
            granuleRec = self.session._GetGranuleTile(granuleid)
            if granuleRec:                
                paramD = ['orbitid', 'acqdate', 'acqtime', 'sunazimuth', 'sunelevation', 'doy', 'source', 'product', 'folder', 'filetype', 'filename', 'downloaded', 'organized', 'exploded', 'deleted', 'declouded', 'maskstatus', 'metacheck', 'tgnote']
                self.tilequery = dict(zip(paramD,granuleRec))
                popL = []
                for key in self.tilequery:
                    if self.tilequery[key] == None:
                        popL.append(key)
                for key in popL:
                    self.tilequery.pop(key)
            else:
                exit('No metadata for granule in _ExplodeSentinel2Granules') 
                
        #Insert the meta and tile data in the postgres db
        self.session._InstertTileMeta(self.metaquery)
        '''
            
    def _InsertScene(self,session):
        '''
        '''
        sceneQD = {'lsatprodid':self.lsatprodid,'source':self.source,
                   'product':self.product, 'acqdate':self.acqdate,'acqdatestr':self.acqdatestr,'acqdoy':self.acqdoy,
                   'collnr':self.collnr,'collcat':self.collcat,
                   'wrspath':self.wrspath,'wrsrow':self.wrsrow,
                   'proccat':self.proccat}
        session._InsertScene(sceneQD)
        
    def _UpdateSceneStatus(self,session,statusD):
        statusD['lsatprodid'] = self.lsatprodid
        session._UpdateSceneStatus(statusD)

    def _Explode(self,explodeD,ziptype,filetype):
        """Explode layers and files from downloaded and organized landsat files.
        """
        #Checks and/or extracts the extractL of compressed archives, i.e. scenes and files downloaded from internet          
        if ziptype == 'zip':
            return self.UnZip(explodeD)
        elif ziptype in ['gzip']:
            return self.GunZip(explodeD)
        elif ziptype in ['tar.gz']:
            if filetype[0:7] in ['_LE7HDF','_LE5HDF' ]:
                return self.UnTarGzHDF(explodeD)
            else:
                return self.UnTarGz(explodeD)
        elif self.ziptype in ['tar']:
            return self.UnTar(explodeD)
            
    def UnTarGz(self,explodeD):
        """Explodes layers from tar.gz files.
        """
        import tarfile

        #print ('opening %(tar)s...' %{'tar':self.FN})
        #for l in explodeD:
        #    print ('explode',l, explodeD[l])
        #with tarfile.open(self.FPN, 'r:*') as tarextractL:
        #    for member in tarextractL.getmembers():
        #        print ('memberfull',member.name)
        #        filename = os.path.basename(member.name)
        #        print ('memberfn',filename)
        #    for fn in tarextractL.getnames():
        #        print ('names',fn)
        #BALLE
        try:
            with tarfile.open(self.FPN, 'r:*') as tarextractL:
                for member in tarextractL.getmembers():
                    filepath = member.name
                    filename = os.path.basename(member.name)
                    if not filename or 'gap_mask' in filepath or filepath.endswith('gz'):
                        continue  
                    print ('filename',filename)
                    layer = self._CreateLayerOut(filename,explodeD)
    
                    if layer:
                        print ('extracting', layer.FPN)
                        source = tarextractL.extractfile(member)
                        
                        target = open(layer.FPN, "wb")        
                        shutil.copyfileobj(source, target)
        except:
            
            return False
        return True
            
                    
          
    def UnTarGzHDFTEST(self,explodeL):
        """Explode layers from tar.gz files containing hdf layers.
        """
        import tarfile
        import time   
        for e in explodeL:
            print ('    explodeL', e.FN,e.pattern)
            if e.Exists():
                continue
            if 'lndcal' in e.pattern:     
                tempFPN = os.path.join(self.FP,'lndcal.hdf')
            elif 'lndsr' in e.pattern:
                tempFPN = os.path.join(self.FP,'lndsr.hdf')
            else:
                continue
            #copy the file to memory and extract the hdf straight from memory? 
            cmd = '/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/gdal_translate '
            cmd = '%(cmd)s HDF4_EOS:EOS_GRID:"%(hdf)s"%(band)s %(tar)s' %{'cmd':cmd,'hdf':tempFPN,'band':e.hdfGrid, 'tar':e.FPN}
            print (cmd)
            os.system(cmd)
        print ('opening HDF %(tar)s...' %{'tar':self.FN})
        with tarfile.open(self.FPN, 'r:*') as tarexplodeL:
            for member in tarexplodeL.getmembers(): 
                filename = os.path.basename(member.name)
                if not member.isfile():
                    continue 
                print ('filename in zip', filename)
                if os.path.splitext(filename)[1].lower() == '.hdf':
                    explodeLayerL = self.ListHDFLayersOut(filename,explodeL) 
                    if len(explodeLayerL) > 0:       
                        #create a temporary extraction of the hdf file and extract from the temporary file
                        #There are some issues with HDF implementation in GDAL/OGR for python
                        if 'lndcal' in filename:     
                            tempFPN = os.path.join(self.FP,'lndcal.hdf')
                        elif 'lndsr' in filename:
                            tempFPN = os.path.join(self.FP,'lndsr.hdf')
                        else:
                            print ('unknown hdf file', filename)
                        if not os.path.isfile(tempFPN):
                            source = tarexplodeL.extractfile(member)
                            target = file(tempFPN, "wb")
                            shutil.copyfileobj(source, target) 
                            prevfilesize = 0 
                            while True:
                                time.sleep(1.0) # seconds
                                filesize = os.path.getsize(tempFPN)
                                print ('filesize',filesize)
                                if prevfilesize == filesize:
                                    break
                                prevfilesize = filesize
                        for e in explodeLayerL:
                            if e.Exists():
                                continue
                            #copy the file to memory and extract the hdf straight from memory? 
                            cmd = '/Library/Frameworks/GDAL.framework/Versions/1.11/Programs/gdal_translate '
                            cmd = '%(cmd)s HDF4_EOS:EOS_GRID:"%(hdf)s"%(band)s %(tar)s' %{'cmd':cmd,'hdf':tempFPN,'band':e.hdfGrid, 'tar':e.FPN}
                            print (cmd)
                            os.system(cmd)
                            BALLE
                else: #non hdf extraction
                    extract = self.CreateLayerOut(filename,explodeL)
                    if extract:
                        source = tarexplodeL.extractfile(member)
                        target = file(extract, "wb")
                        shutil.copyfileobj(source, target)
                    
    def _CreateLayerOut(self,aFN,explodeD):
        """Creates the target layer within the extractfile method - no public access
        """
        #Check which file this is
        layerok = False
        for l in explodeD:
            if explodeD[l]['params']['pattern'] in aFN:
                layerok = l
                break  
        if layerok:
            #print (explodeD[l]['layer'].FPN)
            if explodeD[l]['layer']._Exists():
                #session._InsertSceneBand(explodeD[l]['layer'])
                return False
            else:
                #print ('exploding',explodeD[l]['layer'].FN)

                return explodeD[l]['layer']        
        else:
            pass
            #print 'WARNING %s not extracted' %(aFN)
         
    def ListHDFLayersOutTEST(self,aFN,extractL):
        """Internal list process for the extraction of hdf files.
        """
        import fnmatch
        extractLayerL = []
        for extractLayer in extractL:
            if fnmatch.fnmatch(aFN,extractLayer.pattern):
                extractLayerL.append(extractLayer)
        if len(extractLayerL) == 0:
            pass
        return extractLayerL
    
class ProcessLandsat():
    'class for all Landsat specific processes'   
    def __init__(self, process, session, verbose):
        self.session = session
        self.verbose = verbose
        self.process = process
        print (self.process.proc.processid) 
        
        #direct to subprocess
        if self.process.proc.processid == 'managebulkmetaurl':
            query = {}
            for item in process.proc.paramsD:
                if item in ['today','creator']:
                    continue
                query[item]= process.proc.paramsD[item]

            self.session._ManageBulkMetaUrl(query)
        elif self.process.proc.processid == 'downloadbulkmeta':
            self._DownLoadBulkMeta()
          
        elif self.process.proc.processid == 'landsatmetadb':
            self._LandsatMetaDB()  
        elif self.process.proc.processid == 'insertbulkmeta':
            self._InsertBulkMetaCsv()
        elif self.process.proc.processid == 'searchLandsatSingleScene':
            self._searchLandsatSingleScene()
        elif self.process.proc.processid == 'explodeLandsatSingleScene':
            self._ExplodeLandsatSingleSceneIni()
        elif self.process.proc.processid == 'landsatscenetemplate':
            self._LandsatSceneTemplate()
        elif self.process.proc.processid == 'landsatbandtemplate':
            self._LandsatLayerTemplate()
        elif self.process.proc.processid == 'landsatsupporttemplate':
            self._LandsatLayerTemplate()
        elif self.process.proc.processid == 'landsatmetatemplate':
            self._LandsatLayerTemplate()
        elif self.process.proc.processid == 'getmetaLandsatSingleScene':
            self._SingleSceneIni('getmeta')
        elif self.process.proc.processid == 'TOARFILandsatSingleScene':
            #self._TOARFILandsatSingleSceneIni()
            self._SingleSceneIni('dntotoa','reflectance')
        elif self.process.proc.processid == 'EmissivityKelvinSingleScene':
            #self._EmissivityKelvinSingleSceneIni()
            self._SingleSceneIni('dntokelvin','emissivity')
        elif self.process.proc.processid == 'DNtoSRFILandsatSingleScene':
            self._SingleSceneIni('dntosrfi','reflectance')
        elif self.process.proc.processid == 'LandsatTMDosSingleScene':
            self._SingleProduct('dostosrfi','reflectance')
        elif self.process.proc.processid == 'SRFIfromDOSLandsatSingleScene':
            self._SingleSceneIni('srfifromdos','reflectance')
        elif self.process.proc.processid == 'dbupdateLandsatSingleScene':
            self._SingleSceneIni('dbupdate',None)
        elif self.process.proc.processid == 'RGBLandsatSingleScene':
            self._SingleSceneIni('rgb','reflectance')
        elif self.process.proc.processid == 'MetaMaskLandsatSingleScene':
            self._SingleSceneIni('metamask','bqa')
            
        else:
            print (self.process.proc.processid)
            SNULLEBULLE
                       
    def _DownLoadBulkMeta(self):  
        query = {}     
        for item in self.process.proc.paramsD:
            if item in ['today','creator','filetype','daysago']:
                continue
            query[item]= self.process.proc.paramsD[item]
        paramL = ['csvurl', 'xmlurl']
        urlL = self.session._SelectBulkMetaUrl(query, paramL)

        if self.process.params.filetype == 'xml':
            srcUrl = urlL[1]
            dstFN = os.path.split(urlL[1])[1]
        else:
            srcUrl = urlL[0]
            dstFN = os.path.split(urlL[0])[1]
        dstFP = os.path.join ('/Volumes',self.process.dstpath.volume,'landsat-bulk-meta')
        print (self.process.dstpath.volume)
        print (dstFP)
        if not os.path.exists(dstFP):
            os.makedirs(dstFP)
        dstFPN = os.path.join(dstFP,dstFN)
        print ('srcUrl',srcUrl)
        print ('dstFPN',dstFPN)
        if os.path.isfile(dstFPN):

            #Get the date of the existing file
            t = os.path.getmtime(dstFPN) 
            filedate = mj_dt.DateFromTmTime(t) 
            today = mj_dt.Today()

            daysOld = mj_dt.DateDiff(today,filedate)

            if daysOld == 0: #never update form the same day
                printstr = '%s already downloaded earlier today' %(dstFPN)
                print (printstr)

            elif daysOld < self.process.params.daysago:
                printstr = 'Recent version of %s already exists' %(dstFPN)
                print (printstr)
            else:
                #move and rename the existing file
                retiredFN = '%s_%s' %(mj_dt.DateToStrDate(filedate),dstFPN)
                retiredFP = os.path.join(self.tarpath.mainpath,'retired')
                if not os.path.exists(retiredFP):
                    os.makedirs(retiredFP)
                retiredFPN = os.path.join(retiredFP,retiredFN)
                shutil.move(dstFPN,retiredFPN)
                printstr = 'downloading bulk meta file %s, will take a while...' %(dstFPN)  
                print (printstr)     
                urllib.urlretrieve(metFileUrl, tarFPN)
        else: #file does not exist
            printstr = 'downloading bulk meta file %s, will take a while...' %(dstFPN)  
            print (printstr)
        
            # Download the file from `url` and save it locally under `file_name`:
            with urllib.request.urlopen(srcUrl) as response, open(dstFPN, 'wb') as out_file:
                copyfileobj(response, out_file)     
        #Explode the downloaded gunzip
        srcFPN = self._Explode(dstFPN)
        if self.process.params.filetype == 'xml':
            pass
        else:
            self._InsertBulkMetaCsv(srcFPN)
        
    def _Explode(self,zipFPN):
        import gzip
        dstFPN = os.path.splitext(zipFPN)[0]

        if os.path.isfile(dstFPN):
            t = os.path.getmtime(dstFPN) 
            filedate = mj_dt.DateFromTmTime(t) 
            today = mj_dt.Today()

            daysOld = mj_dt.DateDiff(today,filedate)

            if daysOld == 0: #never update form the same day
                printstr = '%s already downloaded earlier today' %(dstFPN)
                print (printstr)
                return dstFPN

            elif daysOld < self.process.params.daysago:
                printstr = 'Recent version of %s already exists' %(dstFPN)
                print (printstr)
                return dstFPN
            else:
                os.remove(dstFPN)
        with gzip.open(zipFPN, 'r') as f_in, open(dstFPN, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
        return dstFPN
    
    def _LandsatMetaDB(self):
        for item in self.process.proc.metadef.paramsD: 
            self.session._ManageMetaDb(self.process,self.process.proc.metadef.paramsD[item])
        #Manually run the sql
        
    def _CsvMetaHeader(self,headerL,paramTT):
        testL = headerL[:]
        checkok = True
        self.seqL = []
        self.seqD = {}
        self.attrHeaderD = {}
        lacking = 0
        #Construct the headerL for writing the COPY command to postgresql
        for paramT in paramTT:
            #Create the dfferent csv file names
            if paramT[1] not in self.csvFNL:
                self.attrHeaderD[paramT[1]] = ['sceneid']  
                #self.attrHeaderD[paramT[1]].append('sceneid')
                if paramT[1] == 'main':
                    self.attrHeaderD['main'].append('satelnr')
                    self.attrHeaderD['main'].append('sensorid')

                self.csvFNL.append(paramT[1])  
            self.attrHeaderD[paramT[1]].append(paramT[2])
            
        for paramT in paramTT:
            print ('    paramT',paramT)
            #check the parameters
            if paramT[0].lower() != 'false':                        
                if paramT[0] not in headerL:
                    lacking -= 1
                    if paramT[3].lower() == 'y':
                        printstr = 'FATAL Can not find required parameter %s' %(paramT[0])
                        print (printstr)
                        checkok = False
                    else:
                        printstr = 'Can not find  parameter %s setting it to default' %(paramT[0])
                        print (printstr,paramT)
                        self.seqL.append(lacking)                                               
                        self.seqD[lacking] = paramT
                else:
                    #param in headerL
                    seq = headerL.index(paramT[0])
                    self.seqL.append(seq)
                    index = testL.index(paramT[0])
                    testL.pop(index)    
                     
        print ('checkok',checkok)
        print ('unused items in metadata',testL)

        return checkok
                        
    def _InsertBulkMetaCsv(self,srcFPN):
        #import fileinput
        import csv
        
        srcFP,srcFN = os.path.split(srcFPN)
        copyfile = os.path.splitext(srcFN)[0]

        query = {'sensor':'all'}
        paramTT = self.session._SelectMetaTranslate('csv',query)
        if self.process.params.collection != 'LandsatARD':
            query = {'sensor':'collections'}
            paramTTx = self.session._SelectMetaTranslate('csv',query)
            paramTT.extend(paramTTx)

        csvFNL = []
        for paramT in paramTT:
            #Create the dfferent csv file names
            if paramT[1] not in csvFNL:
                csvFNL.append(paramT[1]) 
        
        print ('csvFNL',csvFNL)
        
        popItemL =[]
        for item in csvFNL:
            FN = '%s-%s.csv' %(copyfile,item)
            FPN = os.path.join(srcFP,FN) 
            print ('checking',item)
            #Check if dst csv files exists, and, if so, the dates
            if os.path.isfile(FPN):
                t = os.path.getmtime(FPN) 
                filedate = mj_dt.DateFromTmTime(t) 
                today = mj_dt.Today()
    
                daysOld = mj_dt.DateDiff(today,filedate)
    
                if daysOld == 0: #never update form the same day
                    printstr = '    %s already processed earlier today' %(FPN)
                    print (printstr)
                    popItemL.append(item)

                elif daysOld < self.process.params.daysago:
                    printstr = '    Recent version of %s already processed' %(FPN)
                    print (printstr)
                    popItemL.append(item)

                else:

                    os.remove(FPN)
            else:
                printstr = '    %s does not exists' %(FPN)
                print (printstr)
           
        for item in popItemL:
            index = csvFNL.index(item)
            csvFNL.pop(index)     
        if not self.process.overwrite and len(csvFNL) == 0:
            return
 
        
        
        self.csvFNL = []
        with open(srcFPN) as f:
            reader = csv.reader(f)
            #Get the header
            headerL = next(reader)
            sceneIdindex = headerL.index('sceneID')
            checkok = self._CsvMetaHeader(headerL,paramTT)

            if not checkok:
                SNULLEBULLE
            #Add acqdoy to main
            self.attrHeaderD['main'].append('acqdoy')
            #append the copy soruce file to all meat tables - to allow bulk deleting
            for key in self.attrHeaderD:          
                self.attrHeaderD[key].append('copyfile')
                
            #Set the csv dst filenames
            metaFPND = {}
            csvFD = {}
            for key in self.csvFNL:
                FN = '%s-%s.csv' %(copyfile,key)
                FPN = os.path.join(srcFP,FN)

                metaFPND[key] = FPN

            #If new files are to be created
            if len(csvFNL) > 0:

                #Create the csv files to use for COPY to db
                for key in self.csvFNL:
                    FN = '%s.csv' %(key)
                    #open the files for writing   
                    F = open(metaFPND[key],'w')
                    csvFD[key] = F             
                    #write the headerL
                    wr = csv.writer(csvFD[key])
                    wr.writerow(self.attrHeaderD[key])
                x = 0
                sceneidL = []
                for row in reader:
                    #Get the sceneid
                    sceneId = row[sceneIdindex]
                    if sceneId in sceneidL:
                        continue
                    else:
                        sceneidL.append(sceneId)
                    #Extract satellite and sensor codes ffom sceneid
                    satelnr = int(sceneId[2:3])
                    sensorId = sceneId[1:2]
                    #Translate acqdate to doy
                    
                    #paramD = dict([ ('main',[sceneId,satelnr,sensorId]), ('geo',[sceneId]), ('sub',[sceneId]), ('url',[sceneId]), ('coll',[sceneId]) )])
                    paramD = dict([ ('main',[sceneId,satelnr,sensorId]), ('geo',[sceneId]), ('sub',[sceneId]), ('url',[sceneId]), ('coll',[sceneId])])
    
                    for i,s in enumerate(self.seqL):
                        #Simplify and sort parameters
                        if s < 0:
                            #no data column in bulk file, use template data
                            paramD[self.seqD[s][1]].append(self.seqD[s][4])     
                        elif row[s]:  
                            paramD[paramTT[i][1]].append(row[s])
                        else:
                            #Null data column in bulk file, use default
                            paramD[paramTT[i][1]].append(paramTT[i][4])
                    #add doy from acddate
                    acqD = mj_dt.yyyy_mm_dd_Str_ToDate(paramD['main'][4])
                    acqdoy = mj_dt.DateToDOY(acqD)
                    paramD['main'].append(acqdoy)
    
                    for key in paramD:
                        paramD[key].append(copyfile)
                        wr = csv.writer(csvFD[key])
                        wr.writerow(paramD[key])        
                    #ConnLandsat.InsertBulkMetaItems(paramD,paramD['sceneId'][1],sensorid,sensorDbD[sensorid])                
                    x += 1
                    if x/100000 == int(x/100000):
                        print ('    ',x)
                for key in csvFD:
                    csvFD[key].close()
        #If this point is reached, insert the data
        self.session._CopyBulkMeta(metaFPND,self.attrHeaderD,copyfile)

    def _searchLandsatSingleScene(self):
        '''seach and download a single scene position
        '''
        print ('    Searching Landsat scene',self.process.params.wrspath,self.process.params.wrsrow)
        #create a temp folder to which the download will be directed, only when the download is complete will the data be moved in place
        self.tempFP = os.path.join('/Volumes',self.process.dstpath.volume, 'landsat', 'temp')
        if not os.path.exists(self.tempFP):
            os.makedirs(self.tempFP)   
        if self.process.params.asscript:
            shFN = '%(sens)s-%(p)d-%(r)d.sh' %{'sens':self.process.params.sensorid,'p':self.process.params.wrspath,'r':self.process.params.wrsrow}
            shFP = os.path.join(self.tempFP, 'script')
            if not os.path.exists(shFP):
                os.makedirs(shFP)
            shFPN = os.path.join(shFP,shFN)
            self.scriptF = open(shFPN,'w')
            
        paramL = ['m.acqdate','m.path','m.row','s.imgqual1','m.datatypel1']

        #Search the data to download
        statusD = {}

        #if self.process.params.downloaded:
        #    statusD['l.downloaded'] = {'val':'Y', 'op':'='}
        #else:
        #statusD['l.downloaded'] = {'val':'N', 'op':'='}
        
        statusD['l.downloaded'] = {'val':'N', 'op':'='}
        if self.process.params.imgqualmin:
            statusD['s.imgqual1'] = {'val':self.process.params.imgqualmin, 'op':'>='}
            
        if self.process.params.collnr:
            statusD['c.collnr'] = {'val':self.process.params.collnr, 'op':'='}
            
        if self.process.params.geormse:
            statusD['g.geormse'] = {'val':self.process.params.geormse, 'op':'<'}
            
        if self.process.params.satelnr:
            statusD['m.satelnr'] = {'val':self.process.params.satelnr, 'op':'='}
            
        if len(self.process.params.datatypel1) > 4:
            statusD['m.datatypel1'] = {'val':self.process.params.datatypel1, 'op':'='}
        
        if len(self.process.params.sensorid) > 0:
            statusD['m.sensorid'] = {'val':self.process.params.sensorid, 'op':'='}
            
        if self.process.params.maxcloudcovland < 100:
            statusD['s.cloudcovland'] = {'val':self.process.params.maxcloudcovland, 'op':'<'}
            paramL.append('s.cloudcovland')
        elif self.process.params.maxcloudcov < 100:
            statusD['s.cloudcov'] = {'val':self.process.params.maxcloudcov, 'op':'<'}
            paramL.append('s.cloudcov')
        else:
            paramL.append('s.cloudcov')
            
        paramL.append('c.lsatprodid')
        paramL.append('m.sceneid')
        statusD['m.path'] = {'val':self.process.params.wrspath, 'op':'='}
        statusD['m.row'] = {'val':self.process.params.wrsrow, 'op':'='}
        scenes = self.session._SelectUSGSLandsatScenes(self.process.srcperiod, statusD, paramL)
        pL = [p.split('.')[1] for p in paramL]
        sceneD = {}
        for scene in scenes:
            sceneD[scene[len(scene)-1]] = dict(zip(pL,scene))
            #All scenes that can be retrieved like this are process category DN
            sceneD[scene[len(scene)-1]]['proccat'] = 'DN'
        if self.process.params.download:
            self._DownloadLandsatScenes(sceneD)
        else:
            self._SeachLandsatScenes(sceneD)

    def _SeachLandsatScenes(self,sceneD):
        #Loop over the scenes and download, organize and explode as requested
        dstFP = os.path.join('/Volumes', self.process.dstpath.volume, 'landsat-bulk-download')

        for scene in sceneD:
            print ('scene',scene,sceneD[scene])
            dstFN = '%(id)s.tar.gz' %{'id':sceneD[scene]['lsatprodid']}
            self.dstFPN = os.path.join(dstFP,dstFN)
            print ('    searching ...',self.dstFPN)
            #Check if scene is in bulk download folder
            if os.path.exists(self.dstFPN):
                print ('    ....not found')
                self._OrganizeDNScene(sceneD[scene])
                continue
            #Check if scene is organized
            lsatScene = self._ConstructScene(sceneD[scene],self.process.dstpath)
            
            if os.path.exists(lsatScene.FPN):
                print ('    ....Already organized')
                lsatScene._GetMetaDataFromProduct(sceneD[scene]['lsatprodid'])
                lsatScene.proccat = 'DN'
                lsatScene._InsertScene(self.session)
                lsatScene._UpdateSceneStatus(self.session,{'column':'downloaded', 'status': 'Y'})

                lsatScene._UpdateSceneStatus(self.session,{'column':'organized', 'status': 'Y'})
                continue
            else:
                print ('    ....Not downloaded')
                
            
    def _DownloadLandsatScenes(self,sceneD):
        '''Download landsat scenes
        '''
        from landsatxplore.earthexplorer import EarthExplorer
        import netrc
        
        #Create the destination folder
        dstFP = os.path.join('/Volumes', self.process.dstpath.volume, 'landsat-bulk-download')
        if not os.path.exists(dstFP):
            os.makedirs(dstFP)

        #Set the connection to Earth Explorer
        HOST = 'earthexplorer'
        secrets = netrc.netrc()
        username, account, password = secrets.authenticators( HOST )
        ee = EarthExplorer(username, password)
        #dstFP = '/Volumes/africa/landsat-bulk-download'
        
        #Loop over the scenes and download, organize and explode as requested
        for scene in sceneD:
            print ('scene',scene,sceneD[scene])
            dstFN = '%(id)s.tar.gz' %{'id':sceneD[scene]['lsatprodid']}
            self.dstFPN = os.path.join(dstFP,dstFN)
            print ('    downloading ...',self.dstFPN)
            #Check if scene is in bulk download folder
            if os.path.exists(self.dstFPN):
                print ('    ....Already downloaded')
                self._OrganizeDNScene(sceneD[scene])
                continue
            #Check if scene is organized
            lsatScene = self._ConstructScene(sceneD[scene],self.process.dstpath)
            
            if os.path.exists(lsatScene.FPN):
                print ('    ....Already organized')
                lsatScene._GetMetaDataFromProduct(sceneD[scene]['lsatprodid'])
                lsatScene.proccat = 'DN'
                lsatScene._InsertScene(self.session)
                lsatScene._UpdateSceneStatus(self.session,{'column':'downloaded', 'status': 'Y'})

                lsatScene._UpdateSceneStatus(self.session,{'column':'organized', 'status': 'Y'})
                continue

            try:
                result = ee.download(scene_id=scene, output_dir=dstFP)
                print ('    ....starting')
                print ('result',result)
                self._OrganizeDNScene(sceneD[scene])
            except:
                print ('    ...could not download, continuing')
                pass


    def _DownloadLandsatESPAbulkOrder(self):
        '''
        '''
        if self.process.params.htmlfn == 'all':
            FP = os.path.join('/Volumes',self.process.srcpath.volume,'earthexplorerorder') 
            FL = [ f for f in os.listdir(FP) if os.path.isfile(os.path.join(FP,f)) ]
        else:
            FPN = os.path.join('/Volumes',self.process.srcpath.volume,'earthexplorerorder',self.process.params.htmlfn)
            if not os.path.exists(FPN):
                exitstr = 'ESPA bulk download file does not exist: %s'  %(FPN)
                exit(exitstr)
            else:
                FL = [FPN]
            
    
        for FPN in FL:
            if FPN[0] == '.':
                continue
            ext = os.path.splitext(FPN)[1]
            if ext.lower() == '.html':
                #open and read the html file
                #FPN = os.path.join(eartexplFP,F)
                for line in open(FPN):
                    if '>Download</a>' in line:
                        line = line.rstrip('\n')
                        url = line.replace('<a href="','') 
                        url = url.replace('">Download</a>','')
                        srcFN = os.path.split(url)[1]
                        print ('srcFN',srcFN )
                        SNULLE
                        #print url
                        tarFPN = os.path.join(tarFP,srcFN)
                        if not os.path.exists(tarFPN):
                            CurlCmd = 'curl -o %(local)s %(url)s' %{'local':tarFPN,'url':url}
                            os.system( CurlCmd )
                        else:
                            print ('already done',srcFN)
         
    def _OrganizeDNScene(self,scene):
        '''
        '''
        lsatScene = self._ConstructScene(scene,self.process.srcpath)
        lsatScene.proccat = 'DN'
        srcFPN = self.dstFPN
        if os.path.isfile(srcFPN):
            lsatScene._GetMetaDataFromProduct(scene['lsatprodid'])
            lsatScene._InsertScene(self.session)
            lsatScene._UpdateSceneStatus(self.session,{'column':'downloaded', 'status': 'Y'})

        if not os.path.exists(lsatScene.FP):
            os.makedirs(lsatScene.FP)
        if not os.path.exists(lsatScene.FPN):
            printstr = 'Moving scene: %(src)s \n    to %(dst)s' %{'src': srcFPN, 'dst':lsatScene.FPN}
            print (printstr)
            move(srcFPN,lsatScene.FPN)
            lsatScene._UpdateSceneStatus(self.session,{'column':'organized', 'status': 'Y'})

    def _ConstructScene(self,scene,scenePath):
        '''
        '''
        sceneCompL = scene['lsatprodid'].split('_')
        source = sceneCompL[0]
        product = sceneCompL[1]
        
        acqdatestr = sceneCompL[3] 
        acqdate = mj_dt.yyyymmddDate(acqdatestr)
        
        wrspath = int(sceneCompL[2][0:3])
        wrsrow = int(sceneCompL[2][3:6])
        
        division = 'scenes'
        system = 'landsat'
        folder = 'original'
        
        compD = {'source':source,'product':product,'folder':folder,'system':system,'division':division}
        comp = LandsatComposition(compD)
        datumD = {'acqdatestr': mj_dt.DateToStrDate(acqdate), 'acqdate':acqdate}

        #Set the filename
        FN = '%(t)s.tar.gz' %{'t':scene['lsatprodid']}
        #Set the locus
        if 'wrspath' in scene:
            wrsD = ConvertLandsatScenesToStr((scene['wrspath'],scene['wrsrow']))
        else:
            wrsD = ConvertLandsatScenesToStr((scene['path'],scene['row']))
        
        loc = wrsD['prstr']
        #Set the locuspath'pstr','rstr','prstr'
        locusPath = os.path.join(wrsD['pstr'],wrsD['rstr'])
        #Construct the locus dictionary
        locusD = {'locus':loc, 'path':locusPath}
        #Construct the landsat scene            
        return  LandsatScene(comp, locusD, datumD, scenePath, FN)
        
    def _ExplodeLandsatRegion(self): 
        pass

    
    def _ExplodeLandsatSingleSceneIni(self):
        '''
        '''
        wrsD = ConvertLandsatScenesToStr((self.process.params.wrspath,self.process.params.wrsrow))
        self._ExplodeLandsatSingleScene(wrsD)
          
    def _ExplodeLandsatSingleScene(self,wrsD):
            '''
            '''           
            #Search for unexploded scenes for this position
            statusD = {}

            statusD['organized'] = {'val':'Y', 'op':'='}
            statusD['exploded'] = {'val':'N', 'op':'='}
            statusD['redundant'] = {'val':'N', 'op':'='}
            statusD['deleted'] = {'val':'N', 'op':'='}

            statusD['wrspath'] = {'val':wrsD['p'], 'op':'='}
            statusD['wrsrow'] = {'val':wrsD['r'], 'op':'='}
            
            statusD['proccat'] = {'val':self.process.params.proccat, 'op':'='}

            if len(self.process.params.source) == 4:
                statusD['source'] = {'val':self.process.params.source, 'op':'='}
            if len(self.process.params.product) == 4:
                statusD['product'] = {'val':self.process.params.product, 'op':'='}
            if len(self.process.params.collcat) == 2:
                statusD['collcat'] = {'val':self.process.params.collcat, 'op':'='}
            if self.process.params.collnr:
                statusD['collnr'] = {'val':self.process.params.collnr, 'op':'='}
            paramL = ['lsatprodid','source','product','proccat','wrspath','wrsrow','acqdate','proccat']
                    #Search the data to download
            scenesLL = self.session._SelectLocalLandsatScenes(self.process.srcperiod, statusD, paramL)
            for sceneL in scenesLL:
                sceneD = dict(zip(paramL,sceneL))
                lsatScene = self._ConstructScene(sceneD,self.process.srcpath)
                lsatScene._GetMetaDataFromProduct(sceneD['lsatprodid'])
                lsatScene.proccat = sceneD['proccat']
                self._ExplodeScene(lsatScene)
            
    def _ExplodeScene(self,lsatScene):
        '''
        '''
        paramL = ['source', 'product', 'folder', 'band', 'prefix', 'suffix',
                  'dataunit', 'celltype','cellnull', 'scalefac', 'offsetadd', 'measure', 'fileext', 'pattern', 'hdfgrid','filecat']

        queryD = {'source':lsatScene.source,'product':lsatScene.product,
                  'collcat':lsatScene.collcat,'collnr':lsatScene.collnr,
                  'proccat':lsatScene.proccat,'required':'Y'}

        layercomps = self.session._SelectSceneTemplate(queryD,paramL)

        nrExpected = len(layercomps)
        
        if len(layercomps) == 0:
            exitStr = 'No templates found for file',lsatScene.FN
            print (exitStr)
            exit(exitStr)
        layerL = []
        for layercomp in layercomps: 
            layerL.append( dict(zip(paramL,layercomp)) )   
        nrExploded = self._SearchExtractLayers(lsatScene,layerL)
        print (nrExploded)
        print (self.explodeD)
        print (lsatScene.FPN)

        #Here is the explosion
        if len(self.explodeD) > 0:
            if not os.path.isfile(lsatScene.FPN):
                warnstr = 'Can not find file %s' %(lsatScene)
                print (warnstr)
                lsatScene._UpdateSceneStatus(self.session,{'column':'downloaded', 'status': 'N'})
                lsatScene._UpdateSceneStatus(self.session,{'column':'organized', 'status': 'N'})

                return
                #exit(warnstr)
            ewaulr = lsatScene._Explode(self.explodeD,self.process.srcpath.hdrfiletype,'')
            if not ewaulr:
                lsatScene._UpdateSceneStatus(self.session,{'column':'downloaded', 'status': 'N'})

                lsatScene._UpdateSceneStatus(self.session,{'column':'organized', 'status': 'N'})
                warnstr = 'Error when expanding Landsat gz: %s' %(lsatScene.FPN)
                print (warnstr)
                os.remove(lsatScene.FPN)
                return
            #Loop over the explodeD, to fill the db
            for l in self.explodeD:
                layer = self.explodeD[l]['layer']
                bandD = self.explodeD[l]['bandD']

                if os.path.isfile(layer.FPN):
                    nrExploded += 1
                    self.session._InsertBand(bandD)
                    if layer.comp.filecat in ['B','S']:
                        self.session._InsertLayer(layer,self.process.overwrite,self.process.delete)
          
        print ('nrExploded',lsatScene.FN,nrExploded,nrExpected)              
        if nrExploded >= nrExpected:
            lsatScene._UpdateSceneStatus(self.session,{'column':'exploded', 'status': 'Y'})
                 
        '''
        if self.process.params.asscript:
            self.scriptF.close()
            printstr =  'To actually explode the MODIS tiles, you have to execute the shell script file:\n    %(fpn)s' %{'fpn':shFPN}
            print (printstr)
        
            #Accessmeta on the fly
            #procAttrParamD['archive'] = archive
            if 'xml' in metaD:
                metaLayer = metaD['xml']
            else:
                exit('Can not get proper metafile in explode')
            self.AccessSingleMetaFile(LOut.comp,archive,metaLayer,explodeD)
        '''

    def _LandsatSceneTemplate(self):
        queryD = {}
        for item in self.process.proc.paramsD:
            if item not in ['creator', 'today']:
                queryD[item] = self.process.proc.paramsD[item]
        self.session._InsertSceneTemplate(queryD, self.process.overwrite, self.process.delete)
    
    def _SearchExtractLayers(self,lsatScene,compDL):
        '''This is (almost) a copy from modis.modis
        '''
        nrExploded = 0
        self.explodeD = {}
        bandparamL = ['lsatprodid','folder','band','prefix','suffix']

        for compD in compDL:
            
            bandvalueL = [lsatScene.lsatprodid,compD['folder'],compD['band'],compD['prefix'],compD['suffix']] 

            bandD = dict(zip(bandparamL, bandvalueL))

            comp = Composition(compD, self.process.system.dstsystem, self.process.system.dstdivision)
            
            datumD = {'acqdatestr': lsatScene.datum.acqdatestr, 'acqdate':lsatScene.datum.acqdate}
            #Set the locus
            locusPath = lsatScene.locuspath
            #Construct the locus dictionary
            locusD = {'locus':lsatScene.locus, 'wrspath': int(lsatScene.locus[1:4]), 'wrsrow':int(lsatScene.locus[5:8]), 'wrspr':lsatScene.locus,'path':locusPath}
            filepath = lambda: None
            filepath.volume = self.process.dstpath.volume; filepath.hdrfiletype = compD['fileext']
            #Create a standard raster layer
            layer = RasterLayer(comp, locusD, datumD, filepath)
            if not layer._Exists() or self.process.overwrite:
                self.explodeD[compD['band']] = {'layer':layer,'params':compD,'bandD':bandD}
            elif layer._Exists():
                nrExploded += 1
                self.session._InsertBand(bandD)
                if compD['filecat'] in ['B','S']:
                    self.session._InsertLayer(layer,self.process.overwrite,self.process.delete)
        return nrExploded

    def _LandsatLayerTemplate(self):
        queryD = {}
        for item in self.process.proc.paramsD:
            if item not in ['creator', 'today']:
                queryD[item] = self.process.proc.paramsD[item]
        self.session._InsertLayerTemplate(queryD, self.process.overwrite, self.process.delete)
        
    def _SingleProduct(self,lp,bandcat):
        '''
        '''
        wrsD = ConvertLandsatScenesToStr((self.process.params.wrspath,self.process.params.wrsrow))
        statusD ={}
        statusD['lsatprodid'] = {'val':self.process.params.lsatprodid, 'op':'='}
        paramL = ['lsatprodid','source','product','proccat','wrspath','wrsrow','acqdate','proccat','collcat','collnr']
                #Search the data to download
        scenesLL = self.session._SelectLocalLandsatScenes(self.process.srcperiod, statusD, paramL)

        for sceneL in scenesLL:
            sceneD = dict(zip(paramL,sceneL))

            self._GetSceneLayers(sceneD,lp,bandcat)

        
    def _SingleSceneIni(self,lp,bandcat=None):
        '''
        '''
        wrsD = ConvertLandsatScenesToStr((self.process.params.wrspath,self.process.params.wrsrow))
        self._GetSingleScene(wrsD,lp,bandcat)
                      
    def _GetSingleScene(self,wrsD,lp,bandcat):
        '''
        '''               
        #Search for scenes to process for this position
        statusD = {}

        if lp == 'dbupdate':
            pass
        elif lp == 'getmeta':
            statusD['exploded'] = {'val':'Y', 'op':'='}
            statusD['reflcalmeta'] = {'val':'N', 'op':'='}
        elif lp == 'dntosrfi':
            statusD['exploded'] = {'val':'Y', 'op':'='}
            statusD['reflcalmeta'] = {'val':'Y', 'op':'='}
            statusD['atmcal'] = {'val':'N', 'op':'='}
        elif lp == 'srfifromdos':
            statusD['exploded'] = {'val':'Y', 'op':'='}
            statusD['reflcalmeta'] = {'val':'Y', 'op':'='}
            statusD['atmcal'] = {'val':'Y', 'op':'='}
            statusD['srfi'] = {'val':'N', 'op':'='}
        elif lp == 'metamask':
            statusD['exploded'] = {'val':'Y', 'op':'='}
            statusD['metamask'] = {'val':'N', 'op':'='}
        elif lp == 'rgb':
            if self.process.params.theme == 'dn':
                statusD['rgbdn'] = {'val':'N', 'op':'='}
                statusD['exploded'] = {'val':'Y', 'op':'='}
            elif self.process.params.theme == 'toa':
                statusD['rgbtoa'] = {'val':'N', 'op':'='}
                statusD['toa'] = {'val':'Y', 'op':'='}
            elif self.process.params.theme == 'srfi':
                statusD['rgbsrfi'] = {'val':'N', 'op':'='}
                statusD['srfi'] = {'val':'Y', 'op':'='}
            elif self.process.params.theme == 'tercor':
                statusD['rgbtc'] = {'val':'N', 'op':'='}
                statusD['tercor'] = {'val':'Y', 'op':'='}
            elif self.process.params.theme == 'deveg':
                statusD['rgbdeveg'] = {'val':'N', 'op':'='}
                statusD['deveg'] = {'val':'Y', 'op':'='}
            else:
                SNULLEBULLE
        else:
            print ('lp',lp)
            SNULLE
                    
        statusD['wrspath'] = {'val':wrsD['p'], 'op':'='}
        statusD['wrsrow'] = {'val':wrsD['r'], 'op':'='}
        
        statusD['proccat'] = {'val':self.process.params.proccat, 'op':'='}

        if len(self.process.params.source) == 4:
            statusD['source'] = {'val':self.process.params.source, 'op':'='}
        if len(self.process.params.product) == 4:
            statusD['product'] = {'val':self.process.params.product, 'op':'='}
        if len(self.process.params.collcat) == 2:
            statusD['collcat'] = {'val':self.process.params.collcat, 'op':'='}
        if self.process.params.collnr:
            statusD['collnr'] = {'val':self.process.params.collnr, 'op':'='}
            
        paramL = ['lsatprodid','source','product','proccat','wrspath','wrsrow','acqdate','proccat','collcat','collnr']
        '''
        if lp == 'dbupdate':
            paramL = ['lsatprodid','source','product','proccat','wrspath','wrsrow','acqdate','proccat','collcat','collnr',
                      'downloaded','organized','exploded','reflcalmeta','emiscalmeta','atmcal','toarfi','srfi',
                      'emissivity','redundant','deleted','metastatus','maskstatus','tgnote']
        else:
            paramL = ['lsatprodid','source','product','proccat','wrspath','wrsrow','acqdate','proccat','collcat','collnr']
        '''
        scenesLL = self.session._SelectLocalLandsatScenes(self.process.srcperiod, statusD, paramL, True)
        #print (scenesLL)

        for sceneL in scenesLL:
            self.sceneD = dict(zip(paramL,sceneL))

            self.lsatScene = self._ConstructScene(self.sceneD,self.process.srcpath)
            self.lsatScene.lsatprodid = self.sceneD['lsatprodid']
            if lp == 'dbupdate':
                self._DbUpdate(self.sceneD)
            elif lp == 'getmeta':
                self._MetaDataScene(self.sceneD)
            else:
                self._GetSceneLayers(self.sceneD,lp,bandcat)
                
    def _GetSceneLayers(self,sceneD,lp,bandcat):
        '''
        '''
        paramL = ['source', 'product', 'folder', 'band', 'prefix', 'suffix',
                  'dataunit', 'celltype','cellnull', 'scalefac', 'offsetadd', 'measure', 'fileext', 'pattern', 'hdfgrid','filecat']

  
        queryD = {'source':sceneD['source'],'product':sceneD['product'],
                  'collcat':sceneD['collcat'],'collnr':sceneD['collnr'],
                  'proccat':sceneD['proccat'],'required':'Y'}
        if bandcat != None:
            queryD['bandcat'] = bandcat

        layercomps = self.session._SelectSceneTemplate(queryD,paramL)
        
        if len(layercomps) == 0:
            exitStr = 'No templates found for bandcat: %s' %(bandcat)
            print (exitStr)
            SNULLE
            exit(exitStr)
        layerL = []
        for layercomp in layercomps: 
            layerL.append( dict(zip(paramL,layercomp)) )  
        if lp == 'dbupdate':
            return layerL
        self._SelectSrcLayers(sceneD,layerL)
        
        dstD = {}

        #dstD['dntotoa'] = {'folder':'toarfi','prefixsuffix':'toa'}
        #dstD['dntokelvin'] = {'folder':'surfacetemp','prefixsuffix':'kelvin'}
        dstD['dntotoa'] = {'folder':'toarfi','prefixL':[]}
        dstD['dntokelvin'] = {'folder':'surfacetemp','prefixL':[]}
        dstD['dntosrfi'] = {'folder':'dospng','prefixL':['histogram','dospath']}
        dstD['dostosrfi'] = {'folder':'dospng','prefixL':['histogram','dospath']}
        dstD['srfifromdos'] = {'folder':'srfi'}
        dstD['rgb'] = {'folder':'rgb'}
        dstD['metamask'] = {'folder':'mask'}

        if lp == 'srfifromdos': 
            if self.process.params.suffix == 'auto':
                self.suffix = self._SetSRFImodelSuffix()
            else:
                self.suffix = self.process.params.suffix
            allExist = self._SetDstLayers(sceneD,layerL,dstD[lp]['folder'],'srfi',self.suffix)

            if allExist and not self.process.overwrite:
                #Update all the the bands and scene status
                for band in self.dstLayerD:

                    bandD = {'lsatprodid':self.lsatScene.lsatprodid,
                             'folder':self.dstLayerD[band].comp.folder,
                             'band':self.dstLayerD[band].comp.band,
                             'prefix':self.dstLayerD[band].comp.prefix,
                             'suffix':self.dstLayerD[band].comp.suffix} 

                    self.session._InsertBand(bandD)
                    self.session._InsertLayer(self.dstLayerD[band],self.process.overwrite,self.process.delete)

                self.lsatScene._UpdateSceneStatus(self.session,{'column':'srfi', 'status': 'Y'})
                return
        elif lp == 'rgb': 
            #Select the input bands
            queryD = {'lsatprodid':self.lsatScene.lsatprodid,'folder':self.process.params.theme}
            suffix = False
            if self.process.params.suffix != 'auto':
                queryD['suffix'] = self.process.params.suffix
                suffix = self.process.params.suffix
            paramL = ['folder','band','prefix','suffix']
            bandRecs = self.session._SelectBandsFromProdId(queryD,paramL)
            self.srcLayerD = {}
 
            for br in bandRecs:
                bandComp = dict(zip(paramL,br))
                if not suffix:
                    suffix = bandComp['suffix']
                elif suffix != bandComp['suffix']:
                    continue
                    
                bandComp['source'] = sceneD['source']
                bandComp['product'] = sceneD['product']
                bandComp['collcat'] = sceneD['collcat']
                bandComp['collnr'] = sceneD['collnr']
                bandComp['proccat'] = sceneD['proccat']
  
                comp = Composition(bandComp, self.process.system.srcsystem, self.process.system.srcdivision)
                acqdatestr = mj_dt.DateToStrDate(sceneD['acqdate'])
                datumD = {'acqdatestr': acqdatestr, 'acqdate':sceneD['acqdate']}
                #Set the locus
                wrsD = ConvertLandsatScenesToStr((sceneD['wrspath'],sceneD['wrsrow']))
                #Construct the locus dictionary
                locusD = {'locus':wrsD['prstr'], 'wrspath': sceneD['wrspath'], 'wrsrow':sceneD['wrsrow'], 'path':wrsD['prpath']}
                
                #Create a standard layer (that it is raster does not matter, we only want the path)
                layer = RasterLayer(comp, locusD, datumD, self.process.srcpath)
                if os.path.isfile(layer.FPN):
                    self.srcLayerD[bandComp['band']] = layer
                    #self.srcLayerD[bandComp['band']].layer = layer
                    #self.srcLayerD[bandComp['band']].FPN = layer.FPN
                else:
                    BALLE
            #Create the dst layer
            if self.process.params.theme == 'dn':
                bandComp['folder'] = 'rgb'
                bandComp['band'] = 'rgb-dn'
                bandComp['prefix'] = 'rgb-dn'
            elif self.process.params.theme == 'toa':
                bandComp['folder'] = 'rgb'
                bandComp['band'] = 'rgb-toa'
                bandComp['prefix'] = 'rgb-toa'
            elif self.process.params.theme == 'srfi':
                bandComp['folder'] = 'rgb'
                bandComp['band'] = 'rgb-srfi'
                bandComp['prefix'] = 'rgb-srfi'
            elif self.process.params.theme == 'tercor':
                bandComp['folder'] = 'rgb'
                bandComp['band'] = 'rgb-tc'
                bandComp['prefix'] = 'rgb-tc'
            elif self.process.params.theme == 'deveg':
                bandComp['folder'] = 'rgb'
                bandComp['band'] = 'rgb-deveg'
                bandComp['prefix'] = 'rgb-deveg'
            
            comp = Composition(bandComp, self.process.system.srcsystem, self.process.system.srcdivision)

            self.dstLayer = RasterLayer(comp, locusD, datumD, self.process.srcpath)
           
            self.dstBandD = {'lsatprodid':self.lsatScene.lsatprodid,
                'folder':self.dstLayer.comp.folder,
                'band':self.dstLayer.comp.band,
                'prefix':self.dstLayer.comp.prefix,
                'suffix':self.dstLayer.comp.suffix}
            if self.dstLayer._Exists() and not self.process.overwrite:
                self._RegisterRGB()
                return 

        elif lp == 'metamask':
            self.suffix = self.process.params.suffix
            allExist = self._SetDstLayers(sceneD,layerL,dstD[lp]['folder'],'mask',self.suffix)
            
            self.dstBandD = {'lsatprodid':self.lsatScene.lsatprodid,
                         'folder':self.dstLayerD['bqa'].comp.folder,
                         'band':self.dstLayerD['bqa'].comp.band,
                         'prefix':self.dstLayerD['bqa'].comp.prefix,
                         'suffix':self.dstLayerD['bqa'].comp.suffix} 

            if self.dstLayerD['bqa']._Exists() and not self.process.overwrite:
                #Update all the the bands and scene status

                self.session._InsertBand(self.dstBandD)
                self.session._InsertLayer(self.dstLayerD['bqa'],self.process.overwrite,self.process.delete)

                self.lsatScene._UpdateSceneStatus(self.session,{'column':'qamask', 'status': 'Y'})
                return
        else:
            self._SetPngLayers(sceneD,layerL,dstD[lp]['folder'], dstD[lp]['prefixL'] )
        #All of the processes below requires that the bands are availabel:

        for band in self.srcLayerD:
            if not os.path.isfile(self.srcLayerD[band].FPN):
                #SNULLE
                #lsatScene._UpdateSceneStatus(self.session,{'column':'organized', 'status': 'N'})
                #lsatScene._UpdateSceneStatus(self.session,{'column':'exploded', 'status': 'N'})

                return

        if lp == 'dntotoa':
            self._SetDNtoTOARFI(sceneD)
        elif lp =='dntokelvin':
            self._SetDNtoKelvin(sceneD) 
        elif lp == 'dntosrfi':
            self._SetDNtoSRFI(sceneD)
        elif lp == 'dostosrfi':
            self._SetDOStoSRFI(sceneD)
        elif lp == 'srfifromdos':
            self._SetSRFIfromDOS(sceneD)
        elif lp == 'rgb':
            self._ModifiedBroveyRGB()
        elif lp == 'metamask':
            self._MetaMask()
        else:
            SNULLE

    def _SetSRFImodelSuffix(self):
        if self.process.params.powfacdef == self.process.params.darkedgedef:
            return (self.process.params.powfacdef)
        else:
            suffix = '%s-%s' %(self.process.params.powfacdef,self.process.params.darkedgedef)
            return suffix
        
    def _DbUpdate(self,sceneD):
        '''
        '''
        #Check for the following stuff:
        #    downloaded character(1) DEFAULT 'N',
        #    organized character(1) DEFAULT 'N',
        #    exploded character(1) DEFAULT 'N',
        #    reflcalmeta character(1) DEFAULT 'N',
        #    emiscalmeta character(1) DEFAULT 'N',
        #    atmcal character(1) DEFAULT 'N',
        #    toarfi character(1) DEFAULT 'N',
        #    srfi character(1) DEFAULT 'N',
        #    emissivity character(1) DEFAULT 'N',
        #    redundant character(1) DEFAULT 'N',
        #    deleted character(1) DEFAULT 'N',
        #    metastatus character(1) DEFAULT 'N',
        #    maskstatus character(1) DEFAULT 'N',
        #    tgnote varchar(32) DEFAULT 'ok',
        lsatScene = self._ConstructScene(sceneD,self.process.srcpath)
        lsatScene._GetMetaDataFromProduct(sceneD['lsatprodid'])
        sceneFP = os.path.join('/Volumes', self.process.srcpath.volume, 'landsat-bulk-download')
        sceneFN = '%(id)s.tar.gz' %{'id':sceneD['lsatprodid']}
        sceneFPN = os.path.join(sceneFP,sceneFN)
        
        if os.path.exists(sceneFPN):
            print (sceneFPN)
            print ('scene is downloaded')
            lsatScene._UpdateSceneStatus(self.session,{'column':'downloaded', 'status': 'Y'})
            lsatScene._UpdateSceneStatus(self.session,{'column':'organized', 'status': 'N'})

        elif os.path.exists(lsatScene.FPN):
            print (lsatScene.FPN)
            print ('scene is  organized')
            lsatScene._UpdateSceneStatus(self.session,{'column':'organized', 'status': 'Y'})

        else:
            print ('scene is neither downloaded nor organized')
            lsatScene._UpdateSceneStatus(self.session,{'column':'downloaded', 'status': 'N'})
            lsatScene._UpdateSceneStatus(self.session,{'column':'organized', 'status': 'N'})

        #Check if scene is exploded
        bandcat = None
        layerL = self._GetSceneLayers(sceneD,'dbupdate',bandcat)

        self._SelectDstLayers(sceneD,layerL)
        exploded = True
        for band in self.dstLayerD:
            if not os.path.exists(self.dstLayerD[band].FPN):
                exploded = False
            print ('band', band, self.dstLayerD[band].FPN)
        if exploded:
            print ('scene is exploded')
            lsatScene._UpdateSceneStatus(self.session,{'column':'exploded', 'status': 'Y'})
        else:
            print ('scene is not exploded')
            lsatScene._UpdateSceneStatus(self.session,{'column':'exploded', 'status': 'N'})

        #Check is reflection calibration meta is extracted
        bandcat =  'reflectance'
        reflLayerL = self._GetSceneLayers(sceneD,'dbupdate',bandcat)

        self._SelectDstLayers(sceneD,reflLayerL)
        bandExistL = list()
        bandL = list(self.dstLayerD.keys())

        for band in self.dstLayerD:
            if os.path.isfile(self.dstLayerD[band].FPN):
                bandExistL.append(band)

        if all(band in bandExistL for band in bandL):
            print ('scene reflectance meta calibration is extractad')
            lsatScene._UpdateSceneStatus(self.session,{'column':'reflcalmeta', 'status': 'Y'})
        else:
            print ('scene reflectance meta calibration is not extractad')
            lsatScene._UpdateSceneStatus(self.session,{'column':'reflcalmeta', 'status': 'N'})

        #Get the bands recorded in the db table dnreflcal
        queryD = {'lsatprodid':sceneD['lsatprodid']}
        paramL = ['band']
        calL = self.session._SelectBandDos(queryD,paramL)
        calL = [item[0] for item in calL]

        if all(band in calL for band in bandL):
            print ('scene reflectance atmospheric calibration is done')
            lsatScene._UpdateSceneStatus(self.session,{'column':'atmcal', 'status': 'Y'})
        else:
            print ('scene reflectance atmospheric calibration is not done')
            lsatScene._UpdateSceneStatus(self.session,{'column':'atmcal', 'status': 'N'})
        
        #Check is emissivity calibration meta is extracted
        bandcat =  'emissivity'
        emisLayerL = self._GetSceneLayers(sceneD,'dbupdate',bandcat)

        self._SelectDstLayers(sceneD,emisLayerL)
        bandL = list(self.dstLayerD.keys())
        
        calL = self.session._SelectEmissivityCalibration(queryD,paramL)
        calL = [item[0] for item in calL]
        if all(band in calL for band in bandL):
            print ('scene emissivity meta calibration is extractad')
            lsatScene._UpdateSceneStatus(self.session,{'column':'emiscalmeta', 'status': 'Y'})
        else:
            print ('scene emissivity meta calibration is not extractad')
            lsatScene._UpdateSceneStatus(self.session,{'column':'emiscalmeta', 'status': 'N'})
            
        ##################################
        ### Check for the actual bands ###
        
        ###TOA Bands

        self._SetDstLayers(sceneD,reflLayerL,'toarfi','toa',self.process.params.toasuffix)
        dstReflBandL = list(self.dstLayerD.keys()) 
        #print (dstReflBandL)
        dstReflExistL = []
        for band in self.dstLayerD:
            if os.path.isfile(self.dstLayerD[band].FPN):
                dstReflExistL.append(band)

        if all(band in dstReflExistL for band in dstReflBandL):
            print ('scene Top Of Atmosphere reflectance band files exists')
            lsatScene._UpdateSceneStatus(self.session,{'column':'toarfi', 'status': 'Y'})
        else:
            print ('scene Top Of Atmosphere reflectance band files des not exists')
            lsatScene._UpdateSceneStatus(self.session,{'column':'toarfi', 'status': 'N'})
            
        ###SRFI Bands
        self._SetDstLayers(sceneD,reflLayerL,'srfi','srfi',self.process.params.srfisuffix)
        dstReflBandL = list(self.dstLayerD.keys()) 
        #print (dstReflBandL)
        dstReflExistL = []
        for band in self.dstLayerD:
            print (self.dstLayerD[band].FPN)
            if os.path.isfile(self.dstLayerD[band].FPN):
                dstReflExistL.append(band)

        if all(band in dstReflExistL for band in dstReflBandL):
            print ('scene SRFI band files exists')
            lsatScene._UpdateSceneStatus(self.session,{'column':'toarfi', 'status': 'Y'})
        else:
            print ('scene SRFI band files does not exists')
            lsatScene._UpdateSceneStatus(self.session,{'column':'toarfi', 'status': 'N'})
            
        ###Emissivity Bands
        self._SetDstLayers(sceneD,emisLayerL,'surfacetemp','kelvin',self.process.params.emissuffix)
        dstEmisBandL = list(self.dstLayerD.keys()) 
        #print (dstReflBandL)
        dstEmisExistL = []
        for band in self.dstLayerD:

            if os.path.isfile(self.dstLayerD[band].FPN):
                dstEmisExistL.append(band)

        if all(band in dstEmisExistL for band in dstEmisBandL):
            print ('scene Emissivity band files exists')
            lsatScene._UpdateSceneStatus(self.session,{'column':'emissivity', 'status': 'Y'})
        else:
            print ('scene Emissivity band files does not exists')
            lsatScene._UpdateSceneStatus(self.session,{'column':'emissivity', 'status': 'N'})
                        
    def _MetaDataScene(self,sceneD):
        '''
        '''
        paramL = ['folder', 'prefix', 'suffix']
        queryD = {'lsatprodid':sceneD['lsatprodid'], 'band':'mtl'}
        mtlLayer = self.session._SelectSceneBand(queryD,paramL)        
        compD = sceneD
        compD['band'] = 'mtl'
        compD['folder'] = mtlLayer[0]
        compD['prefix'] = mtlLayer[1]
        compD['suffix'] = mtlLayer[2]
        comp = Composition(compD, self.process.system.srcsystem, self.process.system.srcdivision)
        acqdatestr = mj_dt.DateToStrDate(sceneD['acqdate'])
        datumD = {'acqdatestr': acqdatestr, 'acqdate':sceneD['acqdate']}
        #Set the locus
        wrsD = ConvertLandsatScenesToStr((sceneD['wrspath'],sceneD['wrsrow']))
        #Construct the locus dictionary
        locusD = {'locus':wrsD['prstr'], 'wrspath': sceneD['wrspath'], 'wrsrow':sceneD['wrsrow'], 'path':wrsD['prpath']}
        filepath = lambda: None
        filepath.volume = self.process.srcpath.volume; filepath.hdrfiletype = 'txt'
        #Create a standard layer (that it is raster does not matter, we only want the path)
        layer = RasterLayer(comp, locusD, datumD, filepath)
        if os.path.isfile(layer.FPN):
            paramL = ['bandnr','band']
            queryD = {'source':sceneD['source'],'product':sceneD['product'],
                  'collcat':sceneD['collcat'],'collnr':sceneD['collnr'],
                  'proccat':sceneD['proccat'],'bandcat':'reflectance'}
            
            sensorBands = self.session._SelectSceneTemplate(queryD,paramL)
            sensorBandD = {}
            for item in sensorBands:
                sensorBandD[item[0]] = item[1]
            #Change the query to emissivity bands
            queryD['bandcat'] = 'emissivity'
            print ('before')
            thermalBands = self.session._SelectSceneTemplate(queryD,paramL)
            print ('thermalbands',thermalBands)
            thermalBandD = {}
            for item in thermalBands:
                thermalBandD[item[0]] = item[1]
            print ('readingMRL...')
            self._ReadMetaMTL(layer.FPN,sceneD,sensorBandD,thermalBandD)
            print ('..Done')
        else:
            warnstr = '    MTLlayer missing: %s' %(layer.FPN)
            print (warnstr)
            warnstr ='    Landsat product not properly exploded: %s' %(sceneD['lsatprodid'])
            print (warnstr)
            BALLE
            
    def _ReadMetaMTL(self,FPN,sceneD,sensorBandD, thermalBandD):
        """A function to extract the relelvant metadata from the
        USGS control file. Returns dicionaries with LMAX, LMIN,
        QCAL_LMIN and QCAL_LMAX for each of the bands of interest."""
        print ('reading',FPN)
        ######  
  
        groupD = {}
        queryD = {'lsatprodid':sceneD['lsatprodid']} 
        print ('reading IMAGE_ATTRIBUTES')
        groupD['IMAGE_ATTRIBUTES'] = {'CLOUD_COVER':'cloudcov',
                                      'CLOUD_COVER_LAND':'cloudcovland',
                                      'SUN_AZIMUTH':'sunazim', 
                                      'SUN_ELEVATION':'sunelev',
                                      'EARTH_SUN_DISTANCE':'esdist',
                                      'GEOMETRIC_RMSE_MODEL':'geormse'}
        for g in groupD:
            with open(FPN, 'r') as metafile:
                for line in metafile:
                    if g in line:
                        for line in metafile:
                            if g in line:
                                break # stop this inner for loop; outer loop picks up on the next line
                            s = line.split("=")

                            if s[0].strip() in groupD[g]:
                                queryD[ groupD[g][ s[0].strip() ]] = float(s[-1])

        self.session._InsertImageAttributes(queryD)
        
        ###### 
    
        reflmax = {} # Dicts to store constants
        reflmin = {}
        radmax = {} 
        radmin = {}
        qclmax = {}
        qclmin = {}
        radgain = {} 
        radbias = {}
        reflgain = {} 
        reflbias = {}

        groupD = {}
        print ('reading REFLECTANCE CALIBRATION')
        groupD['MIN_MAX_REFLECTANCE'] = {'REFLECTANCE_MAXIMUM_BAND':reflmax, 
                                           'REFLECTANCE_MINIMUM_BAND':reflmin}
        groupD['MIN_MAX_RADIANCE'] = {'RADIANCE_MAXIMUM_BAND':radmax, 
                                           'RADIANCE_MINIMUM_BAND':radmin}
        
        groupD['MIN_MAX_PIXEL_VALUE'] = {'QUANTIZE_CAL_MAX_BAND':qclmax, 
                                           'QUANTIZE_CAL_MIN_BAND':qclmin}
        
        groupD['RADIOMETRIC_RESCALING'] = {'RADIANCE_MULT_BAND':radgain, 
                                           'RADIANCE_ADD_BAND':radbias,
                                           'REFLECTANCE_MULT_BAND':reflgain, 
                                           'REFLECTANCE_ADD_BAND':reflbias}
        for g in groupD:
            with open(FPN, 'r') as metafile:
                for line in metafile:
                    if g in line:
                        for line in metafile:
                            if g in line:
                                break # stop this inner for loop; outer loop picks up on the next line
                            s = line.split("=")
                            paramL = s[0].split('_')
                            paramN = '_'.join(paramL[0:len(paramL)-1]).strip()
                            #paramN.strip()
                            #Get the position of the keyword BAND
                            bandPos = paramL.index('BAND')
                            
                            #bandnr = int(paramL[len(paramL)-1])
                            bandnr = int(paramL[bandPos+1])
                            if bandnr in sensorBandD:
                                #print ('paramN',paramN, bandnr, paramL,s[0])
                                #print ()
                                groupD[g][paramN][sensorBandD[bandnr]] = float(s[-1])
        #convert the dicts to per band - that is how the db expects the data
        calD = {'reflmax':reflmax,'reflmin':reflmin,'radmax':radmax,
                'radmin':radmin,'qclmax':qclmax, 'qclmin':qclmin,
                'radgain':radgain, 'radbias':radbias, 'reflgain': reflgain,'reflbias':reflbias}

        bandCalD = {}
        for band in sensorBandD:
            bandCalD[sensorBandD[band]] = {'lsatprodid':sceneD['lsatprodid'], 'band':sensorBandD[band]}
            for c in calD:
                bandCalD[sensorBandD[band]][c] = calD[c][sensorBandD[band]]

        for b in bandCalD:
            self.session._InsertReflectanceCalibration(bandCalD[b])
        self.lsatScene._UpdateSceneStatus(self.session,{'column':'reflcalmeta', 'status': 'Y'}) 
        
        
        ##### Emissivity data
        radmax = {} 
        radmin = {}
        qclmax = {}
        qclmin = {}
        radgain = {} 
        radbias = {}
        k1 = {}
        k2 = {}
        
        groupD = {}
        print ('reading EMISSIVITY CALIBRATION')
        groupD['MIN_MAX_RADIANCE'] = {'RADIANCE_MAXIMUM_BAND':radmax, 
                                           'RADIANCE_MINIMUM_BAND':radmin}
        
        groupD['MIN_MAX_PIXEL_VALUE'] = {'QUANTIZE_CAL_MAX_BAND':qclmax, 
                                           'QUANTIZE_CAL_MIN_BAND':qclmin}
        
        groupD['RADIOMETRIC_RESCALING'] = {'RADIANCE_MULT_BAND':radgain, 
                                           'RADIANCE_ADD_BAND':radbias}
        groupD['THERMAL_CONSTANTS'] = {'K1_CONSTANT_BAND':k1, 
                                           'K2_CONSTANT_BAND':k2}
        for g in groupD:
            with open(FPN, 'r') as metafile:
                for line in metafile:
                    if g in line:
                        for line in metafile:
                            if g in line:
                                break # stop this inner for loop; outer loop picks up on the next line
                            s = line.split("=")
                            s0 =s[0].replace('_VCID_','')
                            paramL = s0.split('_')
                            paramN = '_'.join(paramL[0:len(paramL)-1]).strip()
                            #paramN.strip()

                            #Get the position of the keyword BAND
                            bandPos = paramL.index('BAND')
                            
                            #bandnr = int(paramL[len(paramL)-1])
                            bandnr = int(paramL[bandPos+1])
                            if bandnr in thermalBandD:
                                #print ('paramN',paramN, bandnr, paramL,s[0])
                                #print ()
                                groupD[g][paramN][thermalBandD[bandnr]] = float(s[-1])
        #convert the dicts to per band - that is how the db expects the data
        calD = {'radmax':radmax,
                'radmin':radmin,'qclmax':qclmax, 'qclmin':qclmin,
                'radgain':radgain, 'radbias':radbias, 'k1':k1, 'k2':k2}

        bandCalD = {}
        for band in thermalBandD:
            bandCalD[thermalBandD[band]] = {'lsatprodid':sceneD['lsatprodid'], 'band':thermalBandD[band]}
            for c in calD:
                bandCalD[thermalBandD[band]][c] = calD[c][thermalBandD[band]]
        

        for b in bandCalD:
            self.session._InsertEmissivityCalibration(bandCalD[b])
        
        self.lsatScene._UpdateSceneStatus(self.session,{'column':'emiscalmeta', 'status': 'Y'})
   
    def _SetPngLayers(self,sceneD,compDL,dstfolder,prefixL):
        '''The dstLayers are png files from matplotlib that shows the DOS transformation
        '''
        self.dstLayerD = {}
        #bandparamL = ['lsatprodid','folder','band','prefix','suffix']
        compIni = compDL[0]

        for prefix in prefixL:
            #Here the comp is changed to a more generic form to allow interprocessing of all Landsat
            compD = {'source':compIni['source'],'product':compIni['product']}
            compD['folder'] = dstfolder
            
            compD['band'] = prefix
            compD['prefix'] = prefix
            #The suffix is replaced by the methods used for determined DOS
            #whihc is not known beforehand and thus added after the analysis
            compD['suffix'] = '0'
            #compD['suffix'] = self.process.params.suffix
            
            comp = Composition(compD, self.process.system.dstsystem, 'scenes')
            
            acqdatestr = mj_dt.DateToStrDate(sceneD['acqdate'])
            datumD = {'acqdatestr': acqdatestr, 'acqdate':sceneD['acqdate']}
            #Set the locus
            wrsD = ConvertLandsatScenesToStr((sceneD['wrspath'],sceneD['wrsrow']))
            #Construct the locus dictionary
            locusD = {'locus':wrsD['prstr'], 'wrspr':wrsD['prstr'], 'wrspath': sceneD['wrspath'], 'wrsrow':sceneD['wrsrow'], 'path':wrsD['prpath']}
            filepath = lambda: None
            filepath.volume = self.process.dstpath.volume; filepath.hdrfiletype = 'png'
            #Create a standard layer (that it is raster does not matter, we only want the path)
            layer = RasterLayer(comp, locusD, datumD, filepath)
            self.dstLayerD[compD['band']] = layer
            
    def _SetDstLayers(self,sceneD,compDL,dstfolder,prefixaddon,suffix):
        '''The dstLayers are png files from matplotlib that shows the DOS transformation
        '''

        self.dstLayerD = {}
        allExists = True
        for compInD in compDL:
            #Here the comp is changed to a more generic form to allow interprocessing of all Landsat
            compD = {'source':compInD['source'],'product':compInD['product']}
            compD['folder'] = dstfolder
            
            compD['band'] = compInD['band']
            compD['prefix'] = '%(b)s-%(pa)s' %{'b':compInD['band'], 'pa':prefixaddon}
            compD['suffix'] = suffix
            if 'dbupdate' not in self.process.proc.processid:
                compD['cellnull'] = self.process.params.cellnull
                compD['celltype'] = self.process.params.celltype
                compD['dataunit'] = self.process.params.dataunit
                compD['scalefac'] = self.process.params.scalefac
                compD['offsetadd'] = self.process.params.offsetadd
                compD['measure'] = self.process.params.measure

            comp = Composition(compD, self.process.system.dstsystem, 'scenes')
            
            acqdatestr = mj_dt.DateToStrDate(sceneD['acqdate'])
            datumD = {'acqdatestr': acqdatestr, 'acqdate':sceneD['acqdate']}
            #Set the locus
            wrsD = ConvertLandsatScenesToStr((sceneD['wrspath'],sceneD['wrsrow']))
            #Construct the locus dictionary
            locusD = {'locus':wrsD['prstr'], 'wrspr':wrsD['prstr'], 'wrspath': sceneD['wrspath'], 'wrsrow':sceneD['wrsrow'], 'path':wrsD['prpath']}
            filepath = lambda: None
            filepath.volume = self.process.dstpath.volume; filepath.hdrfiletype = 'tif'
            #Create a standard layer (that it is raster does not matter, we only want the path)
            layer = RasterLayer(comp, locusD, datumD, filepath)
            self.dstLayerD[compD['band']] = layer
            print (compD['band'], layer.FPN)
            if not os.path.isfile(layer.FPN):
                allExists = False

        return allExists

    def _SelectSrcLayers(self,sceneD,compDL):
        '''This is (almost) a copy from _SearcHExtractLayers
        '''
        self.srcLayerD = {}

        for compD in compDL:
            comp = Composition(compD, self.process.system.srcsystem, self.process.system.srcdivision)
            
            acqdatestr = mj_dt.DateToStrDate(sceneD['acqdate'])
            datumD = {'acqdatestr': acqdatestr, 'acqdate':sceneD['acqdate']}
            #Set the locus
            wrsD = ConvertLandsatScenesToStr((sceneD['wrspath'],sceneD['wrsrow']))
            #Construct the locus dictionary
            locusD = {'locus':wrsD['prstr'], 'wrspath': sceneD['wrspath'], 'wrsrow':sceneD['wrsrow'], 'path':wrsD['prpath']}
            filepath = lambda: None
            filepath.volume = self.process.srcpath.volume; filepath.hdrfiletype = self.process.srcpath.hdrfiletype
            #Create a standard layer (that it is raster does not matter, we only want the path)
            layer = RasterLayer(comp, locusD, datumD, filepath)
            self.srcLayerD[ compD['band'] ] = layer
            
    def _SelectDstLayers(self,sceneD,compDL):
        '''This is (almost) a copy from _SearcHExtractLayers
        '''
        self.dstLayerD = {}

        for compD in compDL:
            comp = Composition(compD, self.process.system.dstsystem, self.process.system.dstdivision)
            
            acqdatestr = mj_dt.DateToStrDate(sceneD['acqdate'])
            datumD = {'acqdatestr': acqdatestr, 'acqdate':sceneD['acqdate']}
            #Set the locus
            wrsD = ConvertLandsatScenesToStr((sceneD['wrspath'],sceneD['wrsrow']))
            #Construct the locus dictionary
            locusD = {'locus':wrsD['prstr'], 'wrspath': sceneD['wrspath'], 'wrsrow':sceneD['wrsrow'], 'path':wrsD['prpath']}
            filepath = lambda: None
            filepath.volume = self.process.dstpath.volume; filepath.hdrfiletype = compD['fileext']
            #Create a standard layer (that it is raster does not matter, we only want the path)
            layer = RasterLayer(comp, locusD, datumD, filepath)
            self.dstLayerD[ compD['band'] ] = layer

    def _SetDNtoTOARFI(self,sceneD):
        #Select the calibration parameters for this scene
        calD = {}
        paramL = ['band','esun']            
        queryD = {'source':sceneD['source']}
        esuns = self.session._SelectBandWaveLengths(queryD,paramL)
        esunD = {}
        for esun in esuns:
            esunD[esun[0]] = dict(zip(paramL,esun))
   
        paramL = ['band','reflmax','reflmin','radmax',
                'radmin','qclmax', 'qclmin',
                'radgain', 'radbias', 'reflgain','reflbias']
        queryD = {'lsatprodid':sceneD['lsatprodid']}
        calrecs = self.session._SelectReflectanceCalibration(queryD,paramL)
        if len(calrecs) == 0:
            RERUN

        for calrec in calrecs:
            bandD = dict(zip(paramL,calrec))
            calD[bandD['band']] = bandD
            calD[bandD['band']]['esun'] = esunD[bandD['band']]['esun']

            
        #Get the required metadata
        paramL = ['sunazim','sunelev','esdist']            
        queryD = {'lsatprodid':sceneD['lsatprodid']}
        recs = self.session._SelectImageAttributes(queryD,paramL)
        imgattr = dict(zip(paramL,recs))

        #Loop over all the bands
        for band in self.dstLayerD:
            if self.dstLayerD[band]._Exists():
                continue
            if not os.path.isfile(self.srcLayerD[band].FPN):
                warnstr = 'Input band missing %s' %(self.srcLayerD[band].FPN)
                print (warnstr)
                continue
            print ('    Converting DN to TOA for %s' %(self.srcLayerD[band].FPN))
            print (self.dstLayerD[band].FPN)
            print (calD[band])
            
            self._DNtoTOARFI(sceneD,self.srcLayerD[band],self.dstLayerD[band],calD[band],imgattr)
       
    def _DNtoTOARFI(self,sceneD, srcLayer, dstLayer, calD, imgAttr):
        '''
        '''
        #Read the src layer
        srcLayer.ReadRasterLayer()   
        
        dstLayer.layer = lambda:None
        #Set the np array as the band
        dstLayer.layer.NPBAND = np.zeros(srcLayer.layer.NPBAND.shape)
        
        #esdist = EartSunDist(DOY)
        factor =  np.pi * imgAttr['esdist']**2 / (calD['esun'] * math.sin(math.radians(imgAttr['sunelev'] )))

        if srcLayer.layer.cellnull == None:
            srcLayer.layer.cellnull = srcLayer.comp.cellnull
            
        if self.process.params.method == 'paris':
            pass
            exit('The Jack Paris methods for converson of DN to TOA not implemented')
        elif self.process.params.method == 'reflgainbias':
            REFLECTANCE = (calD['reflgain']*srcLayer.layer.NPBAND + calD['reflbias']) / \
                        math.sin(math.radians(imgAttr['sunelev']))
        elif self.process.params.method == 'radqcal':
            RAD = (( calD['radmax']-calD['radmin']) / \
                    (calD['qclmax']-calD['qclmin'] )) * \
                    (srcLayer.layer.NPBAND - calD['qclmin']) + \
                    calD['radmin']
            
            REFLECTANCE = RAD * factor
        else:
            exit('methods for TOARFILandsat not recognized')
        
        dstLayer.layer.NPBAND = REFLECTANCE * self.process.params.factor
            
        SingleBandMasking(srcLayer, dstLayer)
            
        '''
        if sceneD['source'] in ['LT05','LE07']:
             
            RAD = (( calD['radmax']-calD['radmin']) / \
                    (calD['qclmax']-calD['qclmin'] )) * \
                    (srcLayer.layer.NPBAND - calD['qclmin']) + \
                    calD['radmin']
            
            REFLECTANCE = RAD * factor
  
            REFLECTANCE = (calD['reflgain']*srcLayer.layer.NPBAND + calD['reflbias']) / \
                        math.sin(math.radians(imgAttr['sunelev'])) 
   
            dstLayer.layer.NPBAND = REFLECTANCE * 10000
            
            SingleBandMasking(srcLayer, dstLayer)
            
        elif sceneD['source'] == 'LC08':
            
            #RAD = calD['radgain']*srcLayer.layer.NPBAND + calD['radbias']
            
            REFLECTANCE = (calD['reflgain']*srcLayer.layer.NPBAND + calD['reflbias']) / \
                        math.sin(math.radians(imgAttr['sunelev'])) 
                        
            dstLayer.layer.NPBAND = REFLECTANCE * 10000
            
            SingleBandMasking(srcLayer, dstLayer)
                           
        else:
            SNULLE
        '''
        dstLayer.CopyGeoformatFromSrcLayer(srcLayer.layer)
        #write the results
        dstLayer.CreateDSWriteRasterArray()

        #Register the layer
        self.session._InsertLayer(dstLayer, self.process.overwrite, self.process.delete)
        dstLayer.layer.NPBAND = None
        dstLayer = None
        srcLayer.layer.NPBAND = None
        srcLayer = None  
        #grabage collect
        gc.collect()
 
    def _SetDNtoKelvin(self,sceneD):  
        #Select the calibration parameters for this scene
        calD = {}

        paramL = ['band','radmax',
                'radmin','qclmax', 'qclmin',
                'radgain', 'radbias', 'k1','k2']
        queryD = {'lsatprodid':sceneD['lsatprodid']}
        calrecs = self.session._SelectEmissivityCalibration(queryD,paramL)
        if len(calrecs) == 0:
            RERUN
            
        for calrec in calrecs:
            bandD = dict(zip(paramL,calrec))
            calD[bandD['band']] = bandD

        #Loop over all the bands
        for band in self.dstLayerD:
            if self.dstLayerD[band]._Exists():
                continue
            if not os.path.isfile(self.srcLayerD[band].FPN):
                warnstr = 'Input band missing %s' %(self.srcLayerD[band].FPN)
                print (warnstr)
                continue
            print ('    Converting DN to Kelvin for %s' %(self.srcLayerD[band].FPN))
            #print (self.dstLayerD[band].FPN)
            #print (calD[band])
            
            self._DNtoKelvin(sceneD,self.srcLayerD[band],self.dstLayerD[band],calD[band])

    def _DNtoKelvin(self, sceneD, srcLayer, dstLayer, calD):
        '''
        '''
        #Read the src layer
        srcLayer.ReadRasterLayer()   
        
        dstLayer.layer = lambda:None
        #Set the np array as the band
        dstLayer.layer.NPBAND = np.zeros(srcLayer.layer.NPBAND.shape)
        
        if srcLayer.layer.cellnull == None:
            srcLayer.layer.cellnull = srcLayer.comp.cellnull
            
        if self.process.params.method == 'lut':
            pass
            exit('The LUT methods for converson of DN to TOA not implemented')
        elif self.process.params.method == 'emisgainbias':
            RAD = calD['radgain']*srcLayer.layer.NPBAND + calD['radbias'] 
            
            KELVIN = calD['k2'] / (np.log( ( calD['k1']/ RAD)+1))

        elif self.process.params.method == 'radqcal':
            RAD = (( calD['radmax']-calD['radmin']) / \
                    (calD['qclmax']-calD['qclmin'] )) * \
                    (srcLayer.layer.NPBAND - calD['qclmin']) + \
                    calD['radmin']
                      
            KELVIN = calD['k2'] / (np.log( ( calD['k1']/ RAD)+1))
        else:
            exit('methods for TOARFILandsat not recognized')
        
        dstLayer.layer.NPBAND = KELVIN * self.process.params.factor
            
        SingleBandMasking(srcLayer, dstLayer)
            
        dstLayer.CopyGeoformatFromSrcLayer(srcLayer.layer)
        #write the results
        dstLayer.CreateDSWriteRasterArray()

        #Register the layer
        self.session._InsertLayer(dstLayer, self.process.overwrite, self.process.delete)
        dstLayer.layer.NPBAND = None
        dstLayer = None
        srcLayer.layer.NPBAND = None
        srcLayer = None  
        #grabage collect
        gc.collect()
         
    def _SetDNtoSRFI(self,sceneD): 
        '''
        ''' 
        res = self._SelectCalPal(sceneD)
        if not res:
            return
        calD ,wlD, imgattr = res
        #Loop over all the bands
        #for band in self.dstLayerD:
        #    print (self.dstLayerD[band].FPN)
            
        #Set the atmospheric correction
        atmcorr = AtmCorr(self.process, self.session, self.verbose, self.lsatScene)
        
        atmcorr._DNtoSRFI(sceneD, self.srcLayerD, self.dstLayerD, calD ,wlD, imgattr)
        
    def _SetDOStoSRFI(self,sceneD): 
        '''
        ''' 
        calD ,wlD, imgattr = self._SelectCalPal(sceneD)            
        #Set the atmospheric correction
        atmcorr = AtmCorr(self.process, self.session, self.verbose)
        atmcorr._DOStoSRFI(sceneD, self.srcLayerD, calD ,wlD, imgattr)
             
    def _SetSRFIfromDOS(self,sceneD): 
        '''
        ''' 

        calD ,wlD, imgattr = self._SelectCalPal(sceneD)
         
        #Set the atmospheric correction
        atmcorr = AtmCorr(self.process, self.session, self.verbose, self.lsatScene)
        atmcorr._SRFIfromDOS(sceneD, self.srcLayerD, self.dstLayerD,calD ,wlD, imgattr, self.suffix)
        self.lsatScene._UpdateSceneStatus(self.session,{'column':'srfi', 'status': 'Y'})
        
    def _SelectCalPal(self,sceneD):
        calD = {}
        wlD = {}
        paramL = ['band','esun','wlcenter']            
        queryD = {'source':sceneD['source']}
        esuns = self.session._SelectBandWaveLengths(queryD,paramL)
        esunD = {}
        for esun in esuns:
            esunD[esun[0]] = dict(zip(['band','esun'],esun[0:2]))
            wlD[esun[0]] = dict(zip(['band','wl'],[esun[0],esun[2]]))
   
        paramL = ['band','reflmax','reflmin','radmax',
                'radmin','qclmax', 'qclmin',
                'radgain', 'radbias', 'reflgain','reflbias']
        queryD = {'lsatprodid':sceneD['lsatprodid']}
        calrecs = self.session._SelectReflectanceCalibration(queryD,paramL)
        if len(calrecs) == 0:
            self.lsatScene._UpdateSceneStatus(self.session,{'column':'reflcalmeta', 'status': 'N'})
            self._MetaDataScene(sceneD)
            calrecs = self.session._SelectReflectanceCalibration(queryD,paramL)
            if len(calrecs) == 0:
                return False

        print (list(calD.keys()))
        print (list(esunD.keys()))
        for calrec in calrecs:
            print ('calcrec',calrec)
            bandD = dict(zip(paramL,calrec))
            print (bandD)
            calD[bandD['band']] = bandD
            calD[bandD['band']]['esun'] = esunD[bandD['band']]['esun']

            
        #Get the required metadata
        paramL = ['sunazim','sunelev','esdist']            
        queryD = {'lsatprodid':sceneD['lsatprodid']}
        recs = self.session._SelectImageAttributes(queryD,paramL)
        imgattr = dict(zip(paramL,recs))
        return (calD ,wlD, imgattr)
        
    def _DNtoSRFIOLD(self,sceneD, srcLayer, dstLayer, calD):
        '''
        '''
        #Read the src layer
        srcLayer.ReadRasterLayer()   
        
        dstLayer.layer = lambda:None
        #Set the np array as the band
        dstLayer.layer.NPBAND = np.zeros(srcLayer.layer.NPBAND.shape)
        
        if srcLayer.layer.cellnull == None:
            srcLayer.layer.cellnull = srcLayer.comp.cellnull
            
        if self.process.params.method == 'lut':
            pass
            exit('The LUT methods for converson of DN to TOA not implemented')
        elif self.process.params.method == 'emisgainbias':
            RAD = calD['radgain']*srcLayer.layer.NPBAND + calD['radbias'] 
            
            KELVIN = calD['k2'] / (np.log( ( calD['k1']/ RAD)+1))

        elif self.process.params.method == 'radqcal':
            RAD = (( calD['radmax']-calD['radmin']) / \
                    (calD['qclmax']-calD['qclmin'] )) * \
                    (srcLayer.layer.NPBAND - calD['qclmin']) + \
                    calD['radmin']
                      
            KELVIN = calD['k2'] / (np.log( ( calD['k1']/ RAD)+1))
        else:
            exit('methods for TOARFILandsat not recognized')
        
        dstLayer.layer.NPBAND = KELVIN * self.process.params.factor
            
        SingleBandMasking(srcLayer, dstLayer)
            
        dstLayer.CopyGeoformatFromSrcLayer(srcLayer.layer)
        #write the results
        dstLayer.CreateDSWriteRasterArray()

        #Register the layer
        self.session._InsertLayer(dstLayer, self.process.overwrite, self.process.delete)
        dstLayer.layer.NPBAND = None
        dstLayer = None
        srcLayer.layer.NPBAND = None
        srcLayer = None  
        #grabage collect
        gc.collect()
        
    def _ModifiedBroveyRGB(self):
        import gdal
        '''
        powD = {'r':0.9, 'g':0.89, 'b':0.83}
        limitD = {'r':10000, 'g':10800, 'b':6300}
        '''
        powD = {'r':self.process.params.redpower, 
                'g':self.process.params.greenpower, 
                'b':self.process.params.bluepower}
        limitD = {'r':self.process.params.redlimit, 
                  'g':self.process.params.greenlimit, 
                  'b':self.process.params.bluelimit}
        
        stretchD = {'r':32767/limitD['r'], 'g':32767/limitD['g'], 'b':32767/limitD['b']}

        #bandD = {'bl':bl,'gl':gl, 'rl':rl, 'nir':na, 'mir':mb}
        bandD = {'bl':self.srcLayerD['bl'],
                 'gl':self.srcLayerD['gl'], 
                 'rl':self.srcLayerD['rl'], 
                 'nir':self.srcLayerD['na'], 
                 'mir':self.srcLayerD['mb']}
        arrD = {}
        #Read the metadata (projection and geotrans) orm the rl layer
        self.srcLayerD['rl'].ReadRasterLayer()
        #Open the input files
        for b in bandD:
            #print ('opening',bandD[b].FPN)
            bandfile = gdal.Open(bandD[b].FPN)
            arrD[b] = bandfile.GetRasterBand(1).ReadAsArray()
            bandfile = None
        
        
         
        #Set the Brovey denominator
        denom = (arrD['gl'] + arrD['rl']  + arrD['nir']);
        denom[denom < 0] = 1
        #Calculate red
        r = 10*((arrD['rl']+arrD['mir'])/2) / denom * arrD['rl']
        #Calculate green
        g = 10*((arrD['gl']+arrD['nir'])/2) / denom * arrD['rl']
        #Calcualte blue
        b = 10*((arrD['bl']+arrD['gl'])/2) / denom * arrD['rl']
        #Avoid negative errors in the power function
        r[(arrD['rl']<1) | (arrD['mir']<1)] = 1
        g[(arrD['gl']<1) | (arrD['nir']<1) ] = 1
        b[(arrD['bl']<1) | (arrD['gl']<1) ] = 1
        
        rgbD = {'r':r,'g':g,'b':b}
    
        rpow = np.power(r, powD['r'])

        gpow = np.power(g, powD['g'])

        bpow = np.power(b, powD['b'])

        rgbPowD = {'r':rpow*stretchD['r']/39,
                   'g':gpow*stretchD['g']/39,
                   'b':bpow*stretchD['b']/39}
    
        rgbAL = []
    
        for b in rgbD:
            
            bandData = rgbPowD[b]+1
            bandData[rgbPowD[b] < 0] = 1
            bandData[rgbPowD[b] > 255] = 255
            bandData[arrD['rl'] == -32768] = 0

            bandData = bandData.astype('uint8')

            rgbAL.append(bandData)

        
        rgb_bands = np.asarray([i for i in rgbAL])
        rgb_bands = rgb_bands.transpose([1, 2, 0])
        cols,rows,layers = rgb_bands.shape
        #print (rgb_bands.shape)
        dst_ds = gdal.GetDriverByName('GTiff').Create(self.dstLayer.FPN, rows, cols, layers, gdal.GDT_Byte)
        #dst_ds = gdal.GetDriverByName('JPEG').Create(self.dstLayer.FPN, rows, cols, layers, gdal.GDT_Byte)

        dst_ds.GetRasterBand(1).WriteArray(rgb_bands[:,:,0])   # write r-band to the    raster
        dst_ds.GetRasterBand(2).WriteArray(rgb_bands[:,:,1])   # write g-band to the raster
        dst_ds.GetRasterBand(3).WriteArray(rgb_bands[:,:,2])   # write b-band to the raster
        dst_ds.SetGeoTransform( self.srcLayerD['rl'].layer.geotrans )
        dst_ds.SetProjection( self.srcLayerD['rl'].layer.projection  )
        dst_ds.FlushCache()                     # write to disk
        dst_ds = None
        self._RegisterRGB()
        #Convert to jpg on the fly
        jpgFPN = self.dstLayer.FPN.replace('.tif','.jpg')
        GdalOfToJpeg = GDALstuff(self.dstLayer.FPN, jpgFPN, self.process.params)
        GdalOfToJpeg.TransformOF('JPEG')
    
    def _RegisterRGB(self):
        #Register as band (this is not a layer )
        print ('bandD',self.dstBandD)
    
        self.session._InsertBand(self.dstBandD)
        if self.process.params.theme == 'dn':
            self.lsatScene._UpdateSceneStatus(self.session,{'column':'rgbdn', 'status': 'Y'})
            
        elif self.process.params.theme == 'toa':
            self.lsatScene._UpdateSceneStatus(self.session,{'column':'rgbtoa', 'status': 'Y'})

        elif self.process.params.theme == 'srfi':
            self.lsatScene._UpdateSceneStatus(self.session,{'column':'rgbsrfi', 'status': 'Y'})

        elif self.process.params.theme == 'tercor':
            self.lsatScene._UpdateSceneStatus(self.session,{'column':'rgbtc', 'status': 'Y'})

        elif self.process.params.theme == 'deveg':
            self.lsatScene._UpdateSceneStatus(self.session,{'column':'rgbdeveg', 'status': 'Y'})

    def _MetaMask(self):
        '''
        '''
        print (self.srcLayerD['bqa'].FPN)
        print (self.dstLayerD['bqa'].FPN)
        
 

        conf = LandsatConfidence.medium
    
        # Get mask indicating cloud pixels with high confidence
        cumulative = True
        collection = self.sceneD['collnr']
        masker = LandsatMasker(self.srcLayerD['bqa'].FPN, collection=1)
    
        cloudmask = masker.get_cloud_mask(conf,cumulative)

        shadowmask = masker.get_cloud_shadow_mask(conf,cumulative)

        snowmask = masker.get_snow_mask(conf,cumulative)
        
        fillmask = masker.get_fill_mask()
        
        shadowmask *= 2
        snowmask *= 3
        fillmask *= 4
        
        mask = cloudmask+shadowmask+snowmask+fillmask
        #Reset the cloudmask to let it rule over cloud shadows and snow
        mask[cloudmask == 1] = 1
        
        
        # save the result
            
        masker.save_tif(mask, self.dstLayerD['bqa'].FPN)
        
        #Register the results
        self.session._InsertBand(self.dstBandD)
        self.session._InsertLayer(self.dstLayerD['bqa'],self.process.overwrite,self.process.delete)

        self.lsatScene._UpdateSceneStatus(self.session,{'column':'qamask', 'status': 'Y'})


def DownloadLandsatESPAbulkOrder(htmlfn,volume):
    '''
    '''
    if htmlfn == 'all':
        FP = os.path.join('/Volumes',volume,'Landsat-ESPA-order') 
        FL = [ f for f in os.listdir(FP) if os.path.isfile(os.path.join(FP,f)) ]
    else:
        FPN = os.path.join('/Volumes',volume,'Landsat-ESPA-order',htmlfn)
        if not os.path.exists(FPN):
            exitstr = 'ESPA bulk download file does not exist: %s'  %(FPN)
            exit(exitstr)
        else:
            FL = [FPN]
        
    dstFP = os.path.join('/Volumes',volume,'Landsat-ESPA-download')
    if not os.path.exists(dstFP ):
        os.makedirs(dstFP)
    for FPN in FL:
        if FPN[0] == '.':
            continue
        ext = os.path.splitext(FPN)[1]
        if ext.lower() == '.html':
            #open and read the html file
            #FPN = os.path.join(eartexplFP,F)
            for line in open(FPN):
                if '>Download</a>' in line:
                    line = line.rstrip('\n')
                    url = line.replace('<a href="','') 
                    url = url.replace('">Download</a>','')
                    dstFN = os.path.split(url)[1]
                    print ('dstFN', dstFN)

                    dstFPN = os.path.join(dstFP,dstFN)
                    if not os.path.exists(dstFPN):
                        CurlCmd = 'curl -o %(local)s %(url)s' %{'local':dstFPN,'url':url}
                        print (CurlCmd)
                        os.system( CurlCmd )
                    else:
                        print ('already done',dstFN)
        
if __name__ == "__main__":
    '''
    pass
    htmlfn = '0101905297395.html'
    volume = 'wetland'
    DownloadLandsatESPAbulkOrder(htmlfn,volume)
    '''
    pass
    
    