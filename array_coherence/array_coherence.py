'''
Plot the wiggle and its stockwell transform
Requires the stockwell.smt package

'''

# standard system stuff
from pathlib import Path
import logging
import sys
import types
from math import *

# matplot lib stuff
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('Agg')
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
    AutoMinorLocator)
from matplotlib.ticker import FuncFormatter

from matplotlib.colors import LinearSegmentedColormap
import matplotlib.dates as mdates


import numpy as np
from obspy import Stream,read,UTCDateTime
from obspy.core.util import AttribDict
from obspy.signal.array_analysis import array_processing
from obspy.geodetics.base import gps2dist_azimuth as gdist
from obspy.imaging.cm import obspy_sequential
from scipy.signal import coherence
color={"color": "0.9"}

class array_coherence(object):
    #def __init__(self,sac_file,reffield=None,startsec=None,endsec=None, hp=None,lp=None,npoles=3,outfile='tmp.png',debug=0):
    def __init__(self,sac_file,reffield=None,startsec=None,endsec=None, freqs=None,outfile='tmp.png',debug=0):
        '''
        '''
        self.debug=debug
        self.log=self.setup_log(self.debug)
       
        # just do it 
        self.run(sac_file,reffield,startsec,endsec,freqs,outfile)

    #def run(self,sacfile,reffield,startsec,endsec,hp,lp,npoles,vmax,dbrng,outfile):
    def run(self,sacfiles,reffield,startsec,endsec,freqs,outfile):
        _name=f'{__name__}.run'
       
        self.log.info(f'{_name}: Reading sacfiles ------------------------')
        st=self.read_sacfiles(sacfiles)

        self.log.info(f'{_name}: Trimming wiggle, if needed---------------')
        st=self.reffield_trim(st,reffield,startsec,endsec)

        self.log.info(f'{_name}:: Doing signal processing -----------------')
        st=self.do_sigproc(st,False)

        self.log.info(f'{_name}: Computing stockwell ---------------------')
        ans=self.coher(st)
        
#        self.log.info(f'{_name}: Plotting figure -------------------------')
        self.make_fig(ans,freqs,endsec,outfile)
        #self.make_fig(st[0],stockw,outfile,reffield,startsec,hp,lp,vmax,dbrng)
    def ang180(self,angle):
        angle =  int(angle) % 180; 

        #// force it to be the positive remainder, so that 0 <= angle < 360  
        angle = (angle + 180) % 180; 
        return angle

    def coher(self,st):
        _name=f'{__name__}.coher'
        stime=st.sort(keys=['starttime'],reverse=True)[0].stats.starttime
        etime=st.sort(keys=['endtime'],reverse=False)[0].stats.endtime
        ans=[]
        for i,tri in enumerate(st):
            lati=tri.stats.coordinates.latitude
            loni=tri.stats.coordinates.longitude
            data_i=tri.data
            _idi=tri.id
            st.remove(tri)
            for j,trj in enumerate(st):
                latj=trj.stats.coordinates.latitude
                lonj=trj.stats.coordinates.longitude
                _idj=trj.id
                gd=gdist(lati,loni,latj,lonj)
                m=gd[0]/1000
                a0=gd[1]
                a1=gd[2]
                a2=self.ang180(a0) # interstation azimuth (0-180)
                if m == 0.0:
                   continue 
                #f, Cxy = coherence(data_i, trj.data,fs=trj.stats.sampling_rate,nperseg=32,nfft=64,noverlap=16)
                f, Cxy = coherence(data_i, trj.data,fs=trj.stats.sampling_rate,nperseg=32,noverlap=8)
                self.plot_twosta(f,Cxy,_idi,_idj,m)
                ans.append([m,a2,_idi,_idj,f,Cxy])
#                print(m,Cxy[0:5],f[0:5])
                print(_idi,_idj,m,a0,a1,a2)
        return ans

    def plot_twosta(self,f,Cxy,id1,id2,m):
        fig=plt.figure(figsize=(6,6),dpi=200)
        gs=fig.add_gridspec(1,1)
        ax=fig.add_subplot(gs[0])
        ax.plot(f,Cxy,linewidth=2,c='black',alpha=.6)
        #ax.scatter(f,Cxy,marker='o',linewidth=.1,c='red',edgecolor='black')

        ax.set_xlim(0,20)
        xmajor=5
        xminor=1
        ax.xaxis.set_major_locator(MultipleLocator(xmajor))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
        ax.xaxis.set_minor_locator(MultipleLocator(xminor))
        ax.xaxis.grid(b=True, which="minor", **color)
        ax.xaxis.grid(b=True, which="major", **color)
        ax.set_xlabel('Frequency(Hz')

        ax.set_ylim(0,1.)
        xmajor=0.25
        xminor=0.05
        ax.yaxis.set_major_locator(MultipleLocator(xmajor))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
        ax.yaxis.set_minor_locator(MultipleLocator(xminor))
        ax.yaxis.grid(b=True, which="minor", **color)
        ax.yaxis.grid(b=True, which="major", **color)
        ax.set_ylabel('Coherence')
        ax.set_title(f'Coherence between {id1} and {id2} ({m:0.2} km)',fontdict={'fontsize':8},loc='left')

        outfile=f'cxy_{id1}_{id2}.png'
        plt.savefig(outfile,bbox_inches='tight')
        plt.close()

    def reffield_trim(self,st,reffield,startsec,endsec):
        _name=f'{__name__}.reffield_trim'
        '''
        SAC time is confusing
        '''
        # trim the waveforms if required
        new_st=Stream()
        if not reffield:
            self.log.debug(f'{_name}: No reffield, skipping trim')
            return st
        if not startsec:
            self.log.debug(f'{_name}: No startsec, using waveform starttime as 0 time')
            startsec=0
        for tr in st: # there should be only 1 trace in original Stream

            # attach some sac attributes
            tr.stats.coordinates = AttribDict({
                'latitude': tr.stats.sac.stla,
                'elevation': tr.stats.sac.stel,
                'longitude':tr.stats.sac.stlo})
            try:
                rcut=tr.stats.sac[reffield.lower()]  # This will become 0 time 
                b=tr.stats.sac['b'] # referenced to starttime 
                # If the reffield is not the sac O field, then it is a pick referenced
                # to O time. I think
                if reffield.lower() in 'a' or ref in 't0' or ref in 't1':
                    o=tr.stats.sac['o']
                else:
                    o=0.0
                self.log.debug(f'{_name}: b={b} o={o} reftime:{rcut}')
            except Exception as e:
                self.log.error(f"Exiting, problem with reffield ({reffield}) for \n\t{tr}\n\t{e}")
                sys.exit(0)

           # tstart=tr.stats.starttime - b + rcut +startsec
            shift = (o - b) + rcut + startsec
            self.log.debug(f'{_name}: Cutting {tr.stats.starttime} by {shift} seconds')
            tstart=tr.stats.starttime + shift
            if not endsec:
                tend=tr.stats.endtime
            else:
                tend=tstart+endsec
            tr.trim(tstart,tend)
            new_st+=tr
        self.log.info(f'{_name}: cut st {new_st}')
        new_st.sort(keys=['starttime'],reverse=False)
        return new_st

    def read_sacfiles(self,sacfiles):
        _name=f"{__name__}.read_sacfiles"
        '''
        '''
        st=Stream()
        try:
            for i in sacfiles:
                st+=read(i,type='SAC')
        except Exception as e:
            self.log.error(f'{_name}: Problem reading {i} ... \n{e}')
            sys.exit(0)

        self.log.debug(f'{_name}: Read {i} ')
        return st



    def plot_coherdist(self,ax,ans,freq,winl):
        _name=f'{__name__}.coherdist'
        color={"color": "0.9"}
        colors = [(1,1,1), (0, 0, 1), (0, 1, 0), (1, 0, 0)]
        colors = [(0, 0, 1), (0, 1, 0), (1, 0, 0)]
        cmap = LinearSegmentedColormap.from_list('my_colors', colors, N=6)


        dists=[]
        Cxys=[]
        azs=[]
        freqs=ans[0][4]
        _idx=self.find_nearest(freqs,freq)
        
        for i in ans:
            dists.append(i[0])
            azs.append(i[1])
            Cxys.append(i[5][_idx])
#ax1.scatter(x, data, c=wts, alpha=0.6, edgecolors='none', cmap=cmap)

        scat=ax.scatter(dists,Cxys,c=azs,alpha=0.6,linewidth=0.35,marker='o',s=50,edgecolor='black',cmap=cmap)

    
        # xaxis stuff
        xmajor,xminor=self.tick_stride(np.min(dists),np.max(dists),base=1,prec=2)
        ax.set_xlim(0,np.max(dists)*1.05)
        ax.xaxis.set_major_locator(MultipleLocator(xmajor))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
        ax.xaxis.set_minor_locator(MultipleLocator(xminor))
        ax.xaxis.grid(b=True, which="minor", **color)
        ax.xaxis.grid(b=True, which="major", **color)
        ax.set_xlabel('Interstation Distance (m)')
#        ax.tick_params(labelbottom=False)    
#    
#        # yaxis stuff
        ax.set_ylim(0,1.1)
        ymajor,yminor=self.tick_stride(0,1,base=.1,prec=2)
        ax.yaxis.set_major_locator(MultipleLocator(ymajor))
        ax.yaxis.set_minor_locator(MultipleLocator(yminor))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax.yaxis.grid(b=True, which="major", **color)
        ax.set_ylabel(f'Coherence')

    
        # plot colorbar
        scat.set_cmap(cmap)
        pos1 = ax.get_position()
        newax=[pos1.x0+pos1.width+.01, pos1.y0, 0.01,pos1.height]
        fig=ax.get_figure()
        ax1=fig.add_axes(newax)
        cbar=fig.colorbar(scat, cax=ax1)
        cbar.set_label(f'Inter-sensor azimuth')
        ax1.invert_yaxis()
#        cbar_maj,cbar_min = self.tick_stride(vmin,vmax,base=1,prec=1)
        ax1.yaxis.set_major_locator(MultipleLocator(30))
#
#
#        # Title
        ax.set_title(f'Signal coherence at center frequency: {freqs[_idx]} Hz',fontdict={'fontsize':8},loc='center')
#        ax.set_title(f'Ref:{true_start}',fontdict={'fontsize':8},loc='right')
        return 

    def plot_q(self,in_data,slow):
        from matplotlib.colorbar import ColorbarBase
        from matplotlib.colors import Normalize
        _name=f'{__name__}.plot_p'
       # make output human readable, adjust backazimuth to values between 0 and 360
        t, rel_power, abs_power, baz, slow = in_data.T
        baz[baz < 0.0] += 360

        colors = [(1,1,1), (0, 0, 1), (0, 1, 0), (1, 0, 0)]
        cmap = LinearSegmentedColormap.from_list('my_colors', colors, N=50)

        # choose number of fractions in plot (desirably 360 degree/N is an integer!)
        N = 36
        N2 = 50
        abins = np.arange(N + 1) * 360. / N
        sbins = np.linspace(0, slow, N2 + 1)
        import types
        print(len(abins),len(sbins),type(abins),type(sbins))
        # sum rel power in bins given by abins and sbins
        hist, baz_edges, sl_edges = \
            np.histogram2d(baz, slow, bins=[len(abins), len(sbins)], weights=rel_power)

        # transform to radian
        baz_edges = np.radians(baz_edges)

        # add polar and colorbar axes
        fig = plt.figure(figsize=(6, 6))
        cax = fig.add_axes([0.85, 0.2, 0.05, 0.5])
        ax = fig.add_axes([0.10, 0.1, 0.70, 0.7], polar=True)
        ax.set_theta_direction(-1)
        ax.set_theta_zero_location("N")

        dh = abs(sl_edges[1] - sl_edges[0])
        dw = abs(baz_edges[1] - baz_edges[0])

        # circle through backazimuth
        for i, row in enumerate(hist):
            #bars = ax.bar(left=(i * dw) * np.ones(N2), height=dh * np.ones(N2), width=dw, bottom=dh * np.arange(N2), color=cmap(row / hist.max()))
            bars = ax.bar((i * dw) * np.ones(N2), height=dh * np.ones(N2), width=dw, bottom=dh * np.arange(N2), color=cmap(row / hist.max()))

#        ax.set_xticks(np.linspace(0, 2 * np.pi, 12, endpoint=False))
        ax.set_xticks(np.linspace(0, 2 * np.pi, 36, endpoint=False))
        ax.set_xticklabels(['0°','','','30°','','','60°','','', '90°','','', '120°','','','150°','','', '180°','','','210°','','','240°','','','270°','','','300°','','','330°''','',])

        # set slowness limits
#        ax.set_ylim(0, slow)
#        [i.set_color('grey') for i in ax.get_yticklabels()]
        ColorbarBase(cax, cmap=cmap,
             norm=Normalize(vmin=hist.min(), vmax=hist.max()))     
        # plot colorbar
        plt.savefig('junk2.png',bbox_inches='tight') 
    def plot_p(self,in_data,slow):
        from matplotlib.colorbar import ColorbarBase
        from matplotlib.colors import Normalize
        _name=f'{__name__}.plot_p'
       # make output human readable, adjust backazimuth to values between 0 and 360
        t, rel_power, abs_power, baz, slow = in_data.T
        baz[baz < 0.0] += 360
        print(len(baz),len(t));
        colors = [(1,1,1), (0, 0, 1), (0, 1, 0), (1, 0, 0)]
        cmap = LinearSegmentedColormap.from_list('my_colors', colors, N=50)

        # choose number of fractions in plot (desirably 360 degree/N is an integer!)
        N = 36
        N2 = 25
        abins = np.arange(N + 1) * 360. / N
        sbins = np.linspace(0, slow, N2 + 1)
        import types
        print(len(abins),len(sbins),type(abins),type(sbins))
        # sum rel power in bins given by abins and sbins
        hist, baz_edges, sl_edges = \
            np.histogram2d(baz, slow, bins=[len(abins), len(sbins)], weights=rel_power)

        # transform to radian
        baz_edges = np.radians(baz_edges)

        # add polar and colorbar axes
        fig = plt.figure(figsize=(6, 6))
        cax = fig.add_axes([0.85, 0.2, 0.05, 0.5])
        ax = fig.add_axes([0.10, 0.1, 0.70, 0.7], polar=True)
        ax.set_theta_direction(-1)
        ax.set_theta_zero_location("N")

        dh = abs(sl_edges[1] - sl_edges[0])
        dw = abs(baz_edges[1] - baz_edges[0])

        # circle through backazimuth
        for i, row in enumerate(hist):
            #bars = ax.bar(left=(i * dw) * np.ones(N2), height=dh * np.ones(N2), width=dw, bottom=dh * np.arange(N2), color=cmap(row / hist.max()))
            bars = ax.bar((i * dw) * np.ones(N2), height=dh * np.ones(N2), width=dw, bottom=dh * np.arange(N2), color=cmap(row / hist.max()))

#        ax.set_xticks(np.linspace(0, 2 * np.pi, 12, endpoint=False))
        ax.set_xticks(np.linspace(0, 2 * np.pi, 36, endpoint=False))
        ax.set_xticklabels(['0°','','','30°','','','60°','','', '90°','','', '120°','','','150°','','', '180°','','','210°','','','240°','','','270°','','','300°','','','330°''','',])

        # set slowness limits
#        ax.set_ylim(0, slow)
#        [i.set_color('grey') for i in ax.get_yticklabels()]
        ColorbarBase(cax, cmap=cmap,
             norm=Normalize(vmin=hist.min(), vmax=hist.max()))     
        # plot colorbar
        plt.savefig('junk2.png',bbox_inches='tight') 
   


    def make_fig(self,ans,freqs,winl,outfile):
        _name=f'{__name__}.make_fig'

        for i in freqs:
            fig = plt.figure(figsize=(8, 6),dpi=200)

            gs=fig.add_gridspec(1,1)
            ax=fig.add_subplot(gs[0])
            self.plot_coherdist(ax,ans,i,winl)
            gs.update(wspace=0.05, hspace=0.20)


#        self.plot_q(ans,slow)
            outfile=f'Cxy_dist_f{i}.png'
            plt.savefig(outfile,bbox_inches='tight')
            self.log.debug(f'{_name}: Done plotting')

    def find_nearest(self,array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return idx

    def get_idx(self,arr,b,e):
        _name=f'{__name__}.get_idx'
        if b and e:
            self.log.debug(f'{_name}: indx from {b} to {e}')
            _ind=np.nonzero((arr > b) & (arr < e))[0]
        elif b:
            self.log.debug(f'{_name}: indx from {b} to end')
            _ind=np.nonzero(arr >= b)[0];
        elif e:
            self.log.debug(f'{_name}: indx from {e} to begin')
            _ind=np.nonzero(arr <= e)[0];
        else:
            self.log.debug(f'{_name}: indx from all')
            _ind=np.nonzero(arr > -1)[0]
        return _ind
        
    def plot_stockwell(self,ax,stockw,tr,tshift,hp,lp,vmax,dbrng):
        _name=f'{__name__}.plot_stockwell'

        # white to blue to green to red
        colors = [(1,1,1), (0, 0, 1), (0, 1, 0), (1, 0, 0)]
        cmap = LinearSegmentedColormap.from_list('my_colors', colors, N=50)
 
        nfreqs,ntimes=np.shape(stockw)
        self.log.debug(f'{_name} nfreqs:{nfreqs} ntimes:{ntimes}')
        t=tr.times()+tshift; # time vector
        df=1/(tr.stats.delta*ntimes)

        # figure out freq vector
        if not hp:
            hp=0.0
        freqs=np.cumsum(np.zeros(nfreqs)+df)+hp # freq vector
        self.log.debug(f'{_name} nfreqs={len(freqs)} freqmin={freqs[0]},freqmax={freqs[-1]}')      

        # figure out bounds
        #halfbin_time = tr.stats.delta/2.
        #halfbin_freq = df/2.

        xmin=t[0]
        xmax=t[-1]
#        ymin=freqs[0] - halfbin_freq
        ymin=freqs[0] 
        ymax=freqs[-1]
        extent = (xmin,xmax,ymin,ymax)
        self.log.debug(f'{_name}: setting extent to {extent}')
        # adjust y-axis limits
#        yminax=1/((x[-1]-x[0])/2) # require X number of cycles vs 1/2 window leng
#        msg=f"Freq min: {freqs[0]:0.2f} Halfbin: {halfbin_freq:0.2f} Ymin:{yminax:0.2f}"
#        print(msg)
    
        ampdata=stockw
        if vmax:
            vmax=vmax
            vmin=vmax-60
        else:
            vmax=10*np.log10(np.max(ampdata))
            vmin=vmax - 60 # 50 dB range
#        vmin=np.log10(np.max(ampdata))
#        #vmin=np.log10(np.min(ampdata))*.8
        self.log.info(f'Stockwell mesh Vmin/Vmax range: {vmin}/{vmax}')
    
        colorplot = ax.imshow(10*np.log10(ampdata), interpolation='none',vmax=vmax,vmin=vmin,
                               extent=extent, cmap=cmap,aspect="auto")
    
        # plot colorbar
        colorplot.set_cmap(cmap)
        pos1 = ax.get_position()
        newax=[pos1.x0+pos1.width+.01, pos1.y0, 0.01,pos1.height]
        fig=ax.get_figure()
        ax1=fig.add_axes(newax)
        cbar=fig.colorbar(colorplot, cax=ax1)
        cbar.set_label(r'$10\cdot\log_{10}\left(Amp^{2}/Hz\right)$')
        ax1.invert_yaxis()
        cbar_maj,cbar_min = self.tick_stride(vmin,vmax,base=1,prec=1)
        self.log.debug(f'{_name}: cbar_maj/minor: {cbar_maj}/{cbar_min}')
        ax1.yaxis.set_major_locator(MultipleLocator(cbar_maj))
        ax1.yaxis.set_major_formatter(FormatStrFormatter('%0d'))
        ax1.yaxis.set_minor_locator(MultipleLocator(cbar_min))

        # X-axis
        xmajor,xminor=self.tick_stride(t[0],t[-1])
        ax.xaxis.set_major_locator(MultipleLocator(xmajor))
        ax.xaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
        ax.xaxis.set_minor_locator(MultipleLocator(xminor))
        ax.set_xlabel('Time (s)')

    # Y-axis
        ymajor,yminor=self.tick_stride(ymin,ymax)
        ax.set_ylim(ymin,ymax)
        ax.yaxis.set_major_locator(MultipleLocator(ymajor))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
        ax.yaxis.set_minor_locator(MultipleLocator(yminor))
        ax.set_ylabel('Frequency (Hz)')
    # going log
        formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
        ax.yaxis.set_major_formatter(formatter)

    def nextpow2(self,number): # 
        return np.ceil(np.log2(number))
    
    #def do_sigproc(self,st,norm,hp,lp,npoles):
    def do_sigproc(self,st,norm):
        _name=f'{__name__}.do_sigproc'
        st.merge()
        st.detrend(type='linear')
        for tr in st:
            if norm:
                d=tr.data/np.sqrt(np.sum(tr.data** 2))
                tr.data=d
#            if hp and lp:
#                self.log.debug(f"Bandpass:{hp} to {lp} Hz")
#                tr.filter('bandpass',freqmin=hp,freqmax=lp,corners=3)
#            elif hp:
#                self.log.debug(f"Highpass: {hp} Hz")
#                tr.filter('highpass',freq=hp,corners=3)
#            elif lp:
#                self.log.debug(f"Lowpass: {lp} Hz")
#                tr.filter('lowpass',freq=lp,corners=3)
            tr.taper(0.05)
        return st


    def tick_stride(self,xmin,xmax,prec=2, base=1):
        '''
        Figure out a reasonable tick stride
        '''
        def myround(x, prec=2, base=1):
              return round(base * round(float(x)/base),prec)

        # Need to figure out a good date algo
        if isinstance(xmin,UTCDateTime) and isinstance(xmax,UTCDateTime):
            mrange=xmax.matplotlib_date - xmin.matplotlib_date # mdate days
            fmajor=myround(mrange/5,prec=prec,base=base)
            fminor=fmajor/5
        else:
            mrange=myround((xmax-xmin),prec=prec,base=base)
            fmajor=myround(mrange/4,prec=prec,base=base)
            fminor=fmajor/4
            return fmajor,fminor

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

