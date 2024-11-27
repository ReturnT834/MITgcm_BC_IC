import argparse



def argument():
    parser = argparse.ArgumentParser(description = '''
    Generates 8 files, one for each meteo forcing condition. Hourly data.
    File list is: 
    BC_atemp,
    BC_aqh,
    BC_uwind, 
    BC_vwind, 
    BC_apress, 
    BC_swflux, 
    BC_lwflux, 
    BC_precip.
    
    ''')
    
    parser.add_argument(   '--outmaskfile', '-m',
                                type = str,
                                required = True,
                                help = '/some/path/outmask.nc')

    parser.add_argument(   '--outputdir',"-o",
                                type = str,
                                required = True,
                                help = '/some/path/')

    parser.add_argument(    '--timelist',"-t", 
                                type = str,
                                required = True,
                                help = ''' Path of the file containing times meteo data.''' )

    parser.add_argument(   '--inputdir', '-i',
                                type = str,
                                default = None,
                                required = True,
                                help ='''/some/path/  Directory containg files to interpolate. 
                                '''
                                )
      
    parser.add_argument(   '--nativemask',
                                type = str,
                                default = None,
                                required = True,
                                help = '''NetCDF File name of the mask on meteo data are defined.
                                ''')    
    return parser.parse_args()

args = argument()

from general import *
from bitsea.commons import genUserDateList as DL
import netCDF4
from scipy.interpolate import griddata
try:
    from mpi4py import MPI
    comm  = MPI.COMM_WORLD
    rank  = comm.Get_rank()
    nranks =comm.size
    isParallel = True
except:
    rank   = 0
    nranks = 1
    isParallel = False

def getFilename(INPUTDIR, var, date17):
    '''Returns: 
    filename if filename exist, None if not.
    Searches for names with hours <=24.
    If file not found it searches for hours > 24, in forecast
    Second loop searches increments hours by 1, in order to find file at hour=0
    when previous bulletin file is missing.
     _ALP-3_ECMWF-ERAINT_evaluation_r1i1p1_ICTP-RegCM4-7_fpsconv-x2yn2-v1_1hr_200001010030-200012312330.nc
    '''
    basename='_ALP-3_ECMWF-ERAINT_evaluation_r1i1p1_ICTP-RegCM4-7_fpsconv-x2yn2-v1_1hr_'
    D= DL.readTimeString(date17)
    year=D.year
    if var == 'pr':
        date_str=str(year)+'0101000000-'+str(year+1)+'0101000000'
    else:
        date_str=str(year)+'01010030-'+str(year)+'12312330'
    filename=INPUTDIR+var+'/'+var+basename+date_str+'.nc'
    return filename

#   date0=datetime(int(year),1,1,0,30,0)
#   delta=D-date0
#   frame_index=int(delta.total_seconds()/3600.)
#   return filename,frame_index

def loadMaskRegCM(infile):
    with netCDF4.Dataset(infile, 'r') as dset:
        dims = dset.dimensions
        jpi = dims['x'].size
        jpj = dims['y'].size
        lon=np.array(dset.variables['lon'][:,:])
        lat=np.array(dset.variables['lat'][:,:])
        tmask=np.array(dset.variables['sftlf'][:,:])
        # we are interested in sea points
        tmask[tmask<1.0]=1
        tmask[tmask>1.0]=0
    return jpi,jpj,lon,lat,tmask.astype(bool)

def writeCheckFile():
    checkfile = OUTPUTDIR +  "CHECK/" + var + "." + time + ".nc"
    Mcheck = Minterp.copy()
    missing_value=1.e+20
    Mcheck[~tmask2]=missing_value
            
    NCout =NC.netcdf_file(checkfile,"w")    
    NCout.createDimension("Lon"  ,Mask2.Lon.size)
    NCout.createDimension("Lat"  ,Mask2.Lat.size)
    
    ncvar = NCout.createVariable(var,'f',('Lat','Lon'))
    setattr(ncvar, 'missing_value', missing_value) 
    ncvar[:] = Mcheck
    NCout.close()

#Mask1 = mask(args.nativemask)
jpi,jpj,LON1,LAT1,tmaskRegCM=loadMaskRegCM(args.nativemask)
Mask2 = mask(args.outmaskfile)

nWP = tmaskRegCM.sum()
Points1 = np.zeros((nWP,2), np.float32)
Points1[:,0] = LON1[tmaskRegCM]
Points1[:,1] = LAT1[tmaskRegCM]

LON2,LAT2 = np.meshgrid(Mask2.Lon,Mask2.Lat)
LON2=LON2.astype(np.float32)
LAT2=LAT2.astype(np.float32)
tmask2 = Mask2.tmask[0,:,:]

TIMELIST=file2stringlist(args.timelist)
INPUTDIR=addsep(args.inputdir)
OUTPUTDIR=addsep(args.outputdir)

os.system("mkdir -p " + OUTPUTDIR)
os.system("mkdir -p " + OUTPUTDIR + "CHECK") 

#VARS=['Lat','Lon','atemp','aqh','uwind', 'vwind', 'apress', 'lwdown', 'swdown', 'precip' ]
VARS=['atemp','aqh','uwind', 'vwind', 'swdown', 'lwdown', 'precip' ] # missing apress

VARS_REGCM={'Lat'   : 'Lat',
            'Lon'   : 'Lon',
            'atemp' : 'tas',  # K
            'aqh'   : 'huss', # Kg/Kg
            'uwind' : 'uas',  # m/s
            'vwind' :'vas',   # m/s
            'swdown': 'rsds', # W/m2
            'lwdown': 'rlus', # W/m2
            'precip': 'pr'}   # m/s
#DType=[]
#for var in VARS: DType.append((var,np.float32))

nFrames = len(TIMELIST)

#M   = np.zeros((nFrames,jpj,jpi),dtype=DType)
M0   = np.zeros((nFrames,jpj,jpi),dtype=float)

for var in VARS[rank::nranks]:
    infile = getFilename(INPUTDIR,VARS_REGCM[var],TIMELIST[0])
    print(infile)
#   open file netcdf
    try:
        with netCDF4.Dataset(infile, 'r') as dset:
            M0 = np.array(dset.variables[VARS_REGCM[var]][0:nFrames,:,:])
    except:
        print('****** Wrong meteo format')
    
    outBinaryFile=OUTPUTDIR + "BC_" + var
    print("writing " + outBinaryFile)
    F = open(outBinaryFile,'wb')
    for it, time  in enumerate(TIMELIST):
        FrameMatrix = M0[it,:,:].copy()
        
        
        #import pylab as pl        
        #pl.figure(1); pl.imshow(FrameMatrix); pl.colorbar(); pl.gca().invert_yaxis() ; pl.show(block=False)
        #pl.figure(2); pl.imshow(Mask1.tmask); pl.colorbar(); pl.gca().invert_yaxis() ; pl.show(block=False)
        #FM = FrameMatrix.copy()
        #FM[~Mask1.tmask] = np.NaN
        #pl.figure(3); pl.imshow(FM); pl.colorbar(); pl.gca().invert_yaxis() ; pl.show(block=False)
        #import sys; sys.exit()
        Nearest = griddata(Points1, FrameMatrix[tmaskRegCM], (LON2,LAT2),method='nearest')
        Linear  = griddata(Points1, FrameMatrix[tmaskRegCM], (LON2,LAT2),method='linear').astype(np.float32)
        Linear_Failed = np.isnan(Linear)
        Linear[Linear_Failed]=Nearest[Linear_Failed]
        Minterp = Linear
        writeCheckFile()
          
        F.write(Minterp)
    F.write(Minterp) # Ask if needed?
    F.write(Minterp)
    F.write(Minterp)
    F.close()



