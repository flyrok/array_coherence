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
from obspy import Stream,read,UTCDateTime
from obspy.core.util import AttribDict
from obspy.signal.cross_correlation  import correlate as xcorr 
from obspy.signal.cross_correlation  import xcorr_max
from obspy.geodetics.base import gps2dist_azimuth as gdist

# Scipy Coherence
from scipy.signal import coherence
from mtspec import mt_coherence as mtcoh

from array_coherence.plot_coherence import plot_twosta,plot_wigs


class array_coherence(object):
    def __init__(self,sac_files,ini_file=None,outfile='tmp.png',debug=0):
        '''
        '''
        self.log=self.setup_log(debug)

        self.ini=self.read_ini(ini_file)
        self.sac_files=sac_files

       

    def run_coherence_scipy(self):
        _name='run_coherence_scipy'
       
        self.log.info(f'{_name}:: Reading SAC files -----------------')
        st=self.read_sacfiles(self.sac_files)

        self.log.info(f'{_name}:: Doing signal processing -----------------')
        st=self.do_sigproc(st)

        self.log.info(f'{_name}: Getting start/end timesd---------------')
        st=self.reffield_trim(st)

        self.log.info(f'{_name}: Computing time shifts to align traces')
        st,beam=self.time_align(st)

        plot_wigs(st,beam)
    
        self.coherence_scipy(st)

       # self.log.info(f'{_name}: Computing stockwell ---------------------')
       # ans=self.coher(st)
        
#      #  self.log.info(f'{_name}: Plotting figure -------------------------')
       # self.make_fig(ans,freqs,endsec,outfile)
        #self.make_fig(st[0],stockw,outfile,reffield,startsec,hp,lp,vmax,dbrng)
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
        _name=f'{__name__}.coherence_scipy'
        ans=[]
        st=st_.copy()
        for i,tr_i in enumerate(st):
            # Trace ID
            _idi=tr_i.id

            # Latitude/Longitude of tr_i
            lati=tr_i.stats.coordinates.latitude
            loni=tr_i.stats.coordinates.longitude

            # Cut times of trace
            startcut0=tr_i.stats.coordinates.startcut - tr_i.stats.coordinates.tshift
            endcut0=tr_i.stats.coordinates.endcut - tr_i.stats.coordinates.tshift

            # Check that correlation shift doesn't start before trace time
            if startcut0 < tr_i.stats.starttime:
                self.log.error(f'{_name} Badnessssssss with startcut {_idi} {startcut0} {tr_i.stats.starttime}')
                startcut0=tr_i.stats.starttime

            # Trim out the trace segment for Coherence measure
            tr_i.trim(starttime=startcut0,endtime=endcut0)
            data_i=tr_i.data

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


                startcut1=tr_j.stats.coordinates.startcut - tr_j.stats.coordinates.tshift
                endcut1=tr_j.stats.coordinates.endcut - tr_j.stats.coordinates.tshift
                tr_j.trim(starttime=startcut1,endtime=endcut1)
                data_j=tr_j.data
                
                f, Cxy = coherence(data_i, data_j,fs=tr_j.stats.sampling_rate,nperseg=32,noverlap=8)
                plot_twosta(f,Cxy,_idi,_idj,m)
                ans.append([m,a2,_idi,_idj,f,Cxy])
#                print(m,Cxy[0:5],f[0:5])
                print(f'{_idi} {_idj} {m:0.3f} km {a0:0.1f} {a1:0.1f} {a2:0.1f}')
        return ans

    def beam_form_xcorr(self,st):
        _name='beam_form_xcorr'
        '''
        '''

        # sort wfs by distance
        try:
            st.traces.sort(key=lambda x: x.stats.sac.dist)
        except Exception as e:
            self.log.error(f'{_name} Sorting traces by distance failed. Is sac.dist set?')
            pass

        st_=st.copy()
        for n,tr in enumerate(st_):
            pre_npts=tr.stats.npts
            
            tr.trim(starttime=tr.stats.coordinates.startcut,endtime=tr.stats.coordinates.endcut)
            
            tr.taper(0.05)
            post_npts=tr.stats.npts
            self.log.debug(f'{_name} pre_npts:{pre_npts} post_npts:{post_npts}')
            if n > 0:
                if not dt0 == tr.stats.delta:
                    self.log.error(f'{_name} Sample mismatch!!!')
                if not len0 == tr.stats.npts:
                    self.log.error(f'{_name} NPTS mismatch {sta0}:{len0} {tr.stats.station}:{tr.stats.npts}') 
#                nshift=self.xcorr_fft(mdata,tr.data)
                xc=xcorr(tr0,tr,tr.stats.npts // 2)
                nshift,va=xcorr_max(xc)
                tshift=nshift*dt0
                st[n].stats.starttime+=tshift
                self.log.debug(f'{_name} Time between {sta0}->{tr.stats.station} = {tshift:0.4f}')
            else:
                tr0=tr
                dt0=tr.stats.delta
                sta0=tr.stats.station
                len0=tr.stats.npts
                mdata=tr.data
                tshift=0.0
                nshift=0
        st.stack(time_tol=2,npts_tol=50)
        for tr in st:
            startcut=float(self.ini.get('TIME','startsec'))
            tr.stats.coordinates = AttribDict({
                'startcut':0,
                'endcut':0,
                'nshift':0,
                'tshift':0.0})

        return st.stack()

    def time_align(self,st):
        _name='time_align'
        '''
        '''

        # make a beam
        st_=st.copy()
        beam=self.beam_form_xcorr(st_)
        t=beam[0].times() + float(self.ini.get('TIME','plot_startsec'))
        s=float(self.ini.get('TIME','startsec'))
        e=float(self.ini.get('TIME','startsec'))+float(self.ini.get('TIME','duration'))
        inds=self.get_idx(t,s,e)
        bdata=beam[0].data[inds[0]:inds[-1]]
        # sort wfs by distance
        try:
            st.traces.sort(key=lambda x: x.stats.sac.dist)
        except Exception as e:
            self.log.error(f'{_name} Sorting traces by distance failed. Is sac.dist set?')
            pass

        # loop through traces and xcorr with beam
        st_=st.copy()
        for n,tr in enumerate(st_):

            pre_npts=tr.stats.npts
            tr.trim(starttime=tr.stats.coordinates.startcut,endtime=tr.stats.coordinates.endcut)
            
            tr.taper(0.05)

            post_npts=tr.stats.npts
            self.log.debug(f'{_name} {tr.id} data_len:{pre_npts} xcorr_len:{post_npts}')
            if n > 0: 
                if not dt0 == tr.stats.delta:
                    self.log.error(f'{_name} Sample mismatch {sta0} vs {tr.stats.station}!!!')
                if not len0 == tr.stats.npts:
                    self.log.error(f'{_name} NPTS mismatch {sta0}:{len0} {tr.stats.station}:{tr.stats.npts}') 

                # Do xcorr
                nshift=int(self.xcorr_fft(bdata,tr.data))
                xc=xcorr(tr0,tr,tr.stats.npts // 2)
                nshift,va=xcorr_max(xc)
                tshift=nshift*dt0
                self.log.debug(f'{_name} {sta0}->{tr.stats.station} xcorr_shift={tshift:0.4f}')
            else:
                tr0=tr
                dt0=tr.stats.delta
                sta0=tr.stats.station
                len0=tr.stats.npts
                mdata=tr.data
                tshift=0.0
                nshift=0

            st[n].stats.coordinates['tshift']=tshift
            st[n].stats.coordinates['nshift']=nshift
#            st_+=tmp_tr
        return st,beam


    def xcorr_fft(self,x, y):

        # pad to same len
        n=len(x)
        m=len(y)
        if n > m:
            result= np.zeros(x.shape)
            result[:y.shape[0]]=y
            y=result
        elif m > n:
            result= np.zeros(y.shape)
            result[:x.shape[0]]=x
            x=result

        f1 = np.fft.fft(x)
        # flip the signal of y
        f2 = np.fft.fft(np.flipud(y))
        assert len(f1) == len(f2)
        xc = np.fft.fftshift(np.real(np.fft.ifft(f1 * f2)))
        zero_index = int(len(x) / 2) - 1
        shift = zero_index + np.argmax(xc)
        return shift
            

    def reffield_trim(self,st):
        _name='reffield_trim'
        '''
        SAC time is confusing
        retun the epoch start/end time of the desired window.
        This windows
        tstart=reffield + startec
        tend=tstart+duration
        e.g.
        tstart=tr.starttime
        '''
        parms={'reffield':'o','startsec':0.0,'duration':1.0,'plot_startsec':0.0,'plot_dursec':0.0}
        for i in parms:
            try:
                parms[i]=self.ini.get('TIME',i)
            except Exception as e:
                self.log.error(f'{_name} problem with ini={i}, {e}')
                pass
        # trim the waveforms if required
        new_st=Stream()
        if not parms['reffield']:
            self.log.debug(f'{_name}: No reffield, skipping trim')
            return 1 
        else:
            reffield=parms['reffield'].lower()

        for tr in st: 
            try:
                rcut=tr.stats.sac[reffield]  
                b=tr.stats.sac['b'] # referenced to starttime 
                # If the reffield is not the sac O field, then it is a pick referenced
                # to O time. 
                if reffield in 'a' or reffield in 't0' or reffield in 't1':
                    o=tr.stats.sac['o']
                else:
                    o=0.0
                self.log.debug(f'{_name}: b={b:0.5f} o={o:0.5f} reftime:{rcut:0.5f}')
            except Exception as e:
                self.log.error(f"Exiting, problem with reffield ({reffield}) for \n\t{tr}\n\t{e}")
                sys.exit(0)

            # Trime traces to plot windows
            plot_startsec=float(parms['plot_startsec'])
            plot_dursec=float(parms['plot_dursec'])
            starttime=tr.stats.starttime
            #t0=starttime + ((o - b) + rcut + plot_startsec)
            t0=starttime  + rcut + plot_startsec - b
            t1=t0+plot_dursec
            tr.trim(starttime=t0,endtime=t1)

            # Figure out coherence window
            startsec=float(parms['startsec'])
            duration=float(parms['duration'])
            shift = (o - b) + rcut + startsec
            self.log.debug(f'{_name}: Cutting {tr.stats.starttime} by {shift} seconds')
            tstart=starttime + shift
            if not duration:
                tend=tr.stats.endtime
            else:
                tend=tstart+duration

            # attach some sac attributes
            tr.stats.coordinates = AttribDict({
                'latitude': tr.stats.sac.stla,
                'elevation': tr.stats.sac.stel,
                'longitude':tr.stats.sac.stlo,
                'startcut':tstart,
                'endcut':tend,
                'nshift':0,
                'tshift':0.0})

            new_st+=tr
        new_st.sort(keys=['starttime'],reverse=False)
        self.log.info(f'{_name}: cut st {new_st}')
        return new_st

    def read_sacfiles(self,sacfiles):
        _name=f"read_sacfiles"
        '''
        '''
        st=Stream()
        try:
            for i in sacfiles:
                st+=read(i,type='SAC')
        except Exception as e:
            self.log.error(f'{_name}: Problem reading {i} ... \n{e}')
            sys.exit(0)

        self.log.debug(f'{_name}: Read {st} ')
        return st


    def find_nearest(self,array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    def get_idx(self,arr,b,e):
        _name='get_idx'
        if b >= arr[0] and e <= arr[-1]:
            self.log.debug(f'{_name}: indx from {b} to {e}')
            _ind=np.nonzero((arr >= b) & (arr <= e))[0]
        elif b >= arr[0] and e >= arr[-1]:
            self.log.debug(f'{_name}: indx from {b} to end')
            _ind=np.nonzero(arr >= b)[0];
        elif b <= arr[0] and e <= arr[-1]:
            self.log.debug(f'{_name}: indx from {e} to begin')
            _ind=np.nonzero(arr <= e)[0];
        else:
            self.log.debug(f'{_name}: indx from all')
            _ind=np.nonzero(arr > -1)[0]
        return _ind
        

    def nextpow2(self,number): # 
        return np.ceil(np.log2(number))
    

    def do_sigproc(self,st):
        _name='do_sigproc'
        '''
        Perform some signal process.
        Detrend, taper, filter,normalize
        No check to see if filter parms are ok
        '''

        # Important ini file parameters, with default setting
        parms={'detrend':'linear','taper':0.0,'hp':0.0,'lp':0.0,'npoles':0.0,'passes':1}
        for i in parms:
            try:
                parms[i]=self.ini.get('SIG_PROC',i)
            except Exception as e:
                self.log.error(f'{_name}: {i}: Error {e}')
                pass

        st.merge()
        st.detrend(type=parms['detrend'])

        if parms['taper']:
            st.taper(float(parms['taper']))

        hp=float(parms['hp'])
        lp=float(parms['lp'])
        npoles=int(parms['npoles'])
        passes=int(parms['passes']) # 0=minphase,1=zerophase

        if hp and lp:
            self.log.debug(f"{_name} Bandpass:{hp} to {lp} Hz")
            st.filter('bandpass',freqmin=hp,freqmax=lp,corners=npoles,zerophase=passes)
        elif hp:
            self.log.debug(f"{_name} Highpass: {hp} Hz")
            st.filter('highpass',freq=hp,corners=npoles,zerophase=passes)
        elif lp:
            self.log.debug(f"{_name} Lowpass: {lp} Hz")
            st.filter('highpass',freq=lp,corners=npoles,zerophase=passes)
        else:
            self.log.debug(f'{_name} Not filtering')

        # normalize amplitudes
        for tr in st:
            tr.data=tr.data/np.sqrt(np.sum(tr.data** 2)) 

        return st

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
            self.log.error(f"({_name}) {ini_file} doesn't exist")
            return None
        try:
            ini = configparser.ConfigParser()
            ini.read(ini_file)
            return ini 
        except Exception as e:
            self.log.error(f"{_name} {e}")
            sys.exit(-1)


    def setup_log(self,debug):

        ''' Helper function to set up logging
            at the right debug level
        '''
        # INFO,DEBUG,WARN
        if debug == 0:
            loglevel="WARN"
        elif debug == 1:
            loglevel="INFO"
        else:
            loglevel="DEBUG"

        logging.basicConfig(filename='/dev/null', level=loglevel,
            datefmt="%Y-%j %H:%M:%S", format="%(asctime)s-%(levelname)s %(message)s")
        log=logging.getLogger(__name__)
        ch = logging.StreamHandler()
        ch.setFormatter(logging.Formatter(datefmt="%Y-%j %H:%M:%S",fmt="%(levelname)s %(asctime)s %(message)s"))
        log.addHandler(ch)
        log.setLevel(loglevel)
        print(f' {__name__} loglevel set to {loglevel}:{log.level}')
        logging.getLogger('matplotlib.font_manager').disabled = True
        return log

