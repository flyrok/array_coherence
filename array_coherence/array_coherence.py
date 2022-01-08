'''
Class to compute array coherence

'''

# standard system stuff
from pathlib import Path
import logging
import sys
#import types
#from math import *
import configparser

# Numpy
import numpy as np

# Obspy
from obspy import Stream,read,UTCDateTime,Trace
from obspy.core.util import AttribDict
from obspy.signal.cross_correlation  import correlate as xcorr 
from obspy.signal.cross_correlation  import xcorr_max
from obspy.geodetics.base import gps2dist_azimuth as gdist

# Scipy Coherence
from scipy.signal import coherence
#from mtspec import mt_coherence as mtcoh

from array_coherence.plot_coherence import plot_twosta,plot_wigs,make_fig
from array_coherence.setup_log import setup_log

log=logging.getLogger(__name__)

class array_coherence(object):
    def __init__(self,sac_files,ini_file=None,outfile='tmp.png',debug=None):
        '''
        '''
        
        # load the ini configuration files
        self.ini=self.read_ini(ini_file)
        # set sac files# 
        self.sac_files=sac_files

    def run_coherence_scipy(self):
        _name='run_coherence_scipy'
        '''
        Main driver to run coherence using SciPy's coherence function
        And plot all the results.
        '''
       
        log.info(f'{_name}:: Reading SAC files       -----------------')
        st=self.read_sacfiles(self.sac_files)

        log.info(f'{_name}:: Doing signal processing -----------------')
        st=self.do_sigproc(st)

        log.info(f'{_name}: Trimming adn time aligning ---------------')
        st,beam=self.reffield_trim(st)


        log.info(f'{_name}: Plotting record section ')
        zero_off=float(self.ini.get('TIME','plot_startsec'))
        plot_wigs(st,zero_off,beam=beam,outfig=self.ini.get('PLOT','wffig'))
        
        log.info(f'{_name}: Computing Coherence  ---------------------')
        Cxy_results=self.coherence_scipy(st)
         
        log.info(f'{_name}: Computing octave bands -------------------')
        fc,fl,fu=self.octave_bands()

        log.info(f'{_name}: Plotting figures -------------------------')
        make_fig(Cxy_results,fc,self.ini.get('PLOT','cdist_base'),fls=fl,fus=fu,domean=self.ini.get('PLOT','domean'))
        
        return 1

    def ang180(self,angle):
        '''
        Confine angles to 0-180
        input:
        angle (float) in degrees
        return:
        angle (float) in degrees
        '''
        angle =  int(angle) % 180; 
        #// force it to be the positive remainder, so that 0 <= angle < 360  
        angle = (angle + 180) % 180; 
        return angle

    def coherence_scipy(self,st_):
        _name='coherence_scipy'
        '''
        Compute the Coherence-squared between all traces in stream *st_*

        '''
        Cxy_results=[] # output results
        st=st_.copy() # copy the Stream

        # Check configuration file vars
        nperseg=self.ini.get('COHERENCE','nperseg')
        if nperseg == 'None':  
            nperseg=None
        else:
            nperseg=int(nperseg)

        noverlap=self.ini.get('COHERENCE','noverlap')
        if noverlap == 'None':  
            noverlap=None
        else: 
            noverlap=int(noverlap) 

        nfft=self.ini.get('COHERENCE','nfft')
        if nfft == 'None':  
            nfft=None
        else: 
            nfft=int(nfft)

        # Loop over Stream
        for i,tr_i in enumerate(st):
            # Trace ID
            _idi=tr_i.id
            # Sample rate
            dt_i=tr_i.stats.delta

            # Latitude/Longitude of tr_i
            lati=tr_i.stats.coordinates.latitude
            loni=tr_i.stats.coordinates.longitude

            # Cut times of trace
            startcut0=tr_i.stats.coordinates.startcut 
            endcut0=tr_i.stats.coordinates.endcut 

            # Check that correlation shift doesn't start before trace time, greatly
            if (startcut0 - tr_i.stats.starttime) < -1*dt_i:
                log.warn(f'{_name}: Startcut( {_idi}) less than trace starttime\n{startcut0} < {tr_i.stats.starttime}\ntry extending startsec time')
                startcut0=tr_i.stats.starttime

            # Trim out the trace segment for Coherence measure
            tr_i.trim(starttime=startcut0,endtime=endcut0)
            tr_i.detrend(type='linear')
            tr_i.taper(0.05)

            if not nfft: nfft=2**nextpow2(tr_i.stats.npts)
            data_i=self.pad(tr_i.data,nfft)
            st.remove(tr_i)

            for tr_j in st:
                _idj=tr_j.id

                latj=tr_j.stats.coordinates.latitude
                lonj=tr_j.stats.coordinates.longitude
                gd=gdist(lati,loni,latj,lonj)
                m=gd[0]/1000 # distance (km)
                a0=gd[1] # azimuth sta1->sta2
                a1=gd[2] # azimuth sta2->sta1
                a2=self.ang180(a0) # interstation azimuth (0-180)


                startcut1=tr_j.stats.coordinates.startcut 
                endcut1=tr_j.stats.coordinates.endcut 
                tr_j.trim(starttime=startcut1,endtime=endcut1)
                tr_j.detrend(type='linear')
                tr_j.taper(0.05)
                data_j=self.pad(tr_j.data,nfft)
                fs=tr_j.stats.sampling_rate
                
                f, Cxy = coherence(data_i, data_j, fs=fs ,nperseg=nperseg, noverlap=noverlap, nfft=nfft)
                plot_twosta(f,Cxy,_idi,_idj,m)
                Cxy_results.append([m,a2,_idi,_idj,f,Cxy])
                log.info(f'{_idi:>13s}->{_idj:>13s} D: {m:06.3f}km Az: {a2:03d} MaxC: {Cxy.max():04.2f} Nfrq: {len(f)}')
        return Cxy_results 

    def bform_tshifts(self,st_):
        _name='bform_tshifts'
        '''
        Get the time shifts to time align traces
        for stacking, plotting, coherence 
        '''
        time_shifts={}
        # sort wfs by distance
        try:
            st_.traces.sort(key=lambda x: x.stats.sac.dist, reverse=False)
        except Exception as e:
            log.error(f'{_name} Sorting traces by distance failed. Is sac.dist set?')
            pass

        tshift=0.0
        for n,tr in enumerate(st_):
            if n == 0:
                dt0=tr.stats.delta
                len0=tr.stats.npts
                mdata=tr.data
                # set shift to zero if doing xcorr
            else:
                if not dt0 == tr.stats.delta:
                    log.error(f'{_name} Sample mismatch!!!')
                if not len0 == tr.stats.npts:
                    log.error(f'{_name} NPTS mismatch {sta0}:{len0} {tr.stats.station}:{tr.stats.npts}') 

                # do xcorr or phase align
                if self.ini.getboolean('TIME','xcorr'):
                    nshift=self.xcorr_fft(mdata,tr.data)
                    tshift=nshift*dt0
                else:
                    try:
                        tshift= 0.0 
                    except Exception as e:
                        log.error(f'{_name}: Crash and burn\n\t{e}')
            log.info(f'{_name}: {n:02d} {tr.id} = {tr.stats.sac.dist:0.3f} km tshift = {tshift:0.4f} sec')
            time_shifts[tr.id]=tshift
        return time_shifts 

    def pad(self,data,n):
        _name='pad'
        '''
        Zero-pad array to length *n*
        '''
        if len(data) < n:
            result= np.zeros(n)
            result[:data.shape[0]]=data
            return result
        else:
             return data

    def xcorr_fft(self,x, y):
        _name='xcorr_fft'
        '''
        Compute the freq domain cross correlation 
        between to array.
        If the arrays are not the same length, then 
        the short one is zero-padded to the length of the
        longer one.
        Returns the shift in samples
        '''

        # pad to same len.
        lengths=[len(x),len(y)]
        N=np.argmax(lengths)
        x=self.pad(x,lengths[N])
        y=self.pad(y,lengths[N])

        # do the xcorr ----
        f1 = np.fft.fft(x)
        # flip the signal of y
        f2 = np.fft.fft(np.flipud(y))
        assert len(f1) == len(f2) # this will throw an exception if not true
        xc = np.fft.fftshift(np.real(np.fft.ifft(f1 * f2)))
        zero_index = int(len(x) / 2) - 1
        shift =  np.argmax(xc) - zero_index
        return shift
            

    def reffield_trim(self,st):
        _name='reffield_trim'
        '''
        SAC time is confusing
        return the epoch start/end time of the desired window.
        This function windows the data to:
            tstart=reffield + startec
            tend=tstart+duration

        input parameters read from the ini file
        '''
        parms={'reffield':'b','startsec':0.0,'duration':1.0,'plot_startsec':0.0,'plot_dursec':0.0}

        # loop through parms, and lookup each parms.key in the ini file
        # grab vale if it exists, else use defaults
        for i in parms: 
            try:
                parms[i]=self.ini.get('TIME',i)
            except Exception as e:
                log.error(f'{_name}: problem with ini={i}, {e}')
                pass

        # deep copy
        st_=st.copy()

        if not parms['reffield']:
            log.error(f'{_name}: No reffield, exiting')
            return 0 
        else:
            reffield=parms['reffield'].lower()

        # get time_shifts -----
        for tr_ in st_: 
            # For each trace, compute
            # 1. SAC's zero_time in absolute time
            # 2. the INI's reference field time, which is relative 
            #    to zero_time
            zero_time, zero_offset = self.get_zerotime(tr_,reffield)
            tcoh_start= zero_time + zero_offset + float(parms['startsec'])
            tcoh_end=tcoh_start+float(parms['duration'])
            tr_.trim(starttime=tcoh_start,endtime=tcoh_end)
            tr_.detrend(type='linear')
            
            # attach some important values to sac attribute dict.

        time_shifts=self.bform_tshifts(st_)

        n=0
        for key,val in time_shifts.items():
            tr=st.select(id=key)[0]
            zero_time, zero_offset = self.get_zerotime(tr,reffield)

            # coherence measurement bounds
            tbase=zero_time + zero_offset + val
            tcoh_start= tbase + float(parms['startsec'])
            tcoh_end= tcoh_start+float(parms['duration'])

            # cut/plot bounds
            tstart= tbase + float(parms['plot_startsec'])
            tend=tstart+float(parms['plot_dursec'])
            tr.trim(starttime=tstart,endtime=tend)

            if n == 0:
                beam=tr.copy()
                bdata=tr.data
                beam.stats.station='beam'
                n+=1
            else:
                beam.data+=tr.data
                n+=1
            tr.stats.coordinates = AttribDict({
                'latitude': tr.stats.sac.stla,
                'elevation': tr.stats.sac.stel,
                'longitude':tr.stats.sac.stlo,
                'startcut':tcoh_start,
                'endcut':tcoh_end})

        beam.data/=len(st)

        return st,beam

    def read_sacfiles(self,sacfiles):
        _name=f"read_sacfiles"
        '''
        Read in a list of sac files
        input: list of sac files
        out: ObsPy Stream object
        '''
        st=Stream()
        try:
            for i in sacfiles:
                st+=read(i,type='SAC')
        except Exception as e:
            log.error(f'{_name}: Problem reading {i} ... \n{e}')
            log.error(f'{_name}: Continuing for better/worse')
            pass
            #sys.exit(0)

        log.debug(f'{_name}: Read {st} ')
        return st


    def get_zerotime(self,tr,reffield):
        _name='get_zerotime'
        try:
            if reffield in 'a' or reffield in 'o' or reffield in 't0' or reffield in 't1':
                zero_time=tr.stats.starttime - tr.stats.sac['b'] # sac zero time
                zero_offset=tr.stats.sac[reffield]
            if reffield in 'b':
                zero_time=tr.stats.starttime
                zero_offset=0.0
            log.debug(f'{_name}: {tr.id} ref={reffield} zero_offset={zero_offset}')
            log.debug(f'{_name}: {tr.id} starttime{tr.stats.starttime} zero_time={zero_time} reftime={zero_time+zero_offset}')
        except Exception as e:
            log.error(f"Exiting, problem with reffield ({reffield}) for \n\t{tr}\n\t{e}")
            sys.exit(0)
        return zero_time,zero_offset


    def do_sigproc(self,st):
        _name='do_sigproc'
        '''
        Perform some signal processing
        Detrend, taper, filter,normalize
        No check to see if filter parms are ok
        '''

        # Important ini file parameters, with default setting
        parms={'detrend':'linear','taper':0.0,'hp':0.0,'lp':0.0,'npoles':0.0,'passes':1}
        for i in parms:
            try:
                parms[i]=self.ini.get('SIG_PROC',i)
            except Exception as e:
                log.error(f'{_name}: {i}: Error {e}')
                pass

        # merge traces, in case there are small gaps
        st.merge()

        # detrend
        st.detrend(type=parms['detrend'])

        # taper
        if parms['taper']:
            st.taper(float(parms['taper']))

        # filter parameters
        hp=float(parms['hp'])
        lp=float(parms['lp'])
        npoles=int(parms['npoles'])
        passes=int(parms['passes']) # 0=minphase,1=zerophase

        # apply filter
        if hp and lp:
            log.debug(f"{_name} Bandpass:{hp} to {lp} Hz")
            st.filter('bandpass',freqmin=hp,freqmax=lp,corners=npoles,zerophase=passes)
        elif hp:
            log.debug(f"{_name} Highpass: {hp} Hz")
            st.filter('highpass',freq=hp,corners=npoles,zerophase=passes)
        elif lp:
            log.debug(f"{_name} Lowpass: {lp} Hz")
            st.filter('lowpass',freq=lp,corners=npoles,zerophase=passes)
        else:
            log.debug(f'{_name} Not filtering')

        # normalize data to the power in the traces
        for tr in st:
            tr.data=tr.data/np.sqrt(np.sum(tr.data** 2)) 

        return st

    def octave_bands(self):
        _name=f'octave_bands'
        '''
        return center freqs, fmin,fmax
        e.g.
        fmin=1
        fmax=12
        bw=1/3
        dt=1/samps_per_sec
        '''
       
        fmin=self.ini.getfloat('PLOT','fmin')
        fmax=self.ini.getfloat('PLOT','fmax')
        bw=self.ini.getfloat('PLOT','bw')

        octs=np.log2(fmax/fmin); # number of octaves
        bmax=np.ceil(octs/bw); # maximum number of bands of bw

        fc=fmin*2**(np.arange(0,bmax)*bw) # center freqs
        fl=fc*2**(-bw/2); # freq_low
        fu=fc*2**(+bw/2); # freq_high
        log.info(f'{_name}: Octave nbands={len(fc)}\n{" ".join(str(x) for x in np.round(fc,4))} Hz')

        return fc,fl,fu


    def read_ini(self,ini_file):
        _name=f'{__name__}.read_pf'
        '''
        read the parm/ini file. 
        This function doesn't test if parms are missing or accurate
        str: pf_file
        configparser configuration files
        returns
        object: configparser or -1
        '''
        if not ini_file:
            return None
        if not Path(ini_file).exists():
            log.error(f"({_name}) {ini_file} doesn't exist")
            return None
        try:
            ini = configparser.ConfigParser()
            ini.read(ini_file)
            return ini 
        except Exception as e:
            log.error(f"{_name} {e}")
            sys.exit(-1)


    def find_nearest(self,array, value):
        _name='find_nearest'
        '''
        find the index of the element of *array* that
        is closest to *value*. This returns a single index
        '''
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    def get_idx(self,arr,b,e):
        _name='get_idx'
        '''
        Retern a list of indexes to array *arr* that meet the
        criteria that the *arr* values are between *b* and *e*
        '''

        if b >= arr[0] and e <= arr[-1]:
            log.debug(f'{_name}: indx from {b} to {e}')
            _ind=np.nonzero((arr >= b) & (arr <= e))[0]
        elif b >= arr[0] and e >= arr[-1]:
            log.debug(f'{_name}: indx from {b} to end')
            _ind=np.nonzero(arr >= b)[0];
        elif b <= arr[0] and e <= arr[-1]:
            log.debug(f'{_name}: indx from {e} to begin')
            _ind=np.nonzero(arr <= e)[0];
        else:
            log.debug(f'{_name}: indx from all')
            _ind=np.nonzero(arr > -1)[0]
        return _ind


    def nextpow2(self,number): # 
        '''
            returns (n) the next power of two, ie. 2**n

        '''
        return np.ceil(np.log2(number))
    
