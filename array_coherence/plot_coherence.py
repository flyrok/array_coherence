'''

'''

# standard system stuff
import sys
from math import *

# Numpy
import numpy as np
# Obspy
from obspy import Stream,read,UTCDateTime
#from obspy.core.util import AttribDict
#from obspy.signal.array_analysis import array_processing
from obspy.geodetics.base import gps2dist_azimuth as gdist
#from obspy.imaging.cm import obspy_sequential


# matplot lib stuff
import matplotlib as mpl
import matplotlib.pyplot as plt
mpl.use('Agg')
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
    AutoMinorLocator)
from matplotlib.ticker import FuncFormatter
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.dates as mdates

from matplotlib.colorbar import ColorbarBase
from matplotlib.colors import Normalize

color={"color": "0.9"}
fontdict_axis={'fontsize':12 }
fontdict_title={'fontsize':12 }
axis_tick_size=12


def plot_wigs(st,zero_off,beam=None,zscl=0.0,env=False,outfig='recsection.png'):
    '''
    '''


    fig=plt.figure(figsize=(6,6),dpi=200)
    gs=fig.add_gridspec(1,1)
    ax=fig.add_subplot(gs[0])

    if beam:
        yb=10*beam.data/np.sqrt(np.sum(beam.data** 2)) *.5
        xb=beam.times()+zero_off
    
    st.traces.sort(key=lambda x: x.stats.sac.dist, reverse=False)
    for n,tr in enumerate(st):

        # segment for coherence and xcorr
        tr_=tr.copy()
        tr_=tr_.slice(starttime=tr.stats.coordinates.startcut,endtime=tr.stats.coordinates.endcut)
        

        x_ = (tr_.stats.starttime - tr.stats.starttime) + tr_.times()    + zero_off
        x = tr.times() + zero_off

        inds=get_idx(x,x_[0],x_[-1])
        
        dist=tr.stats.sac.dist

        if zscl == 0.0:
            y=5*tr.data/np.sqrt(np.sum(tr.data** 2)) 
        else:
            y=(tr.data*zscl)

        if env:
            yenv=-1*filt.envelope(y)
            yenv+=dist

        y+=n
        if n == 0:
            xmin=x[0]
            xmax=x[-1]
            ymin=min(y)

        ax.plot(x,y,linewidth=0.25,color='midnightblue')
        ax.plot(x[inds],y[inds],linewidth=0.25,color='red')
#        ax.text(x=xmin, y=n,s=f'{tr_.id}' , ha="right",fontsize=10)
        ax.annotate(f'{tr_.id}', xy=(xmin,n), xytext=(-5, 0), textcoords="offset points", ha='right', va='center')
        ax.annotate(f'{tr_.stats.sac.dist:0.2f} km', xy=(xmax,n), xytext=(+5, 0), textcoords="offset points", ha='left', va='center')


        if env:
            ax.plot(x,yenv,linewidth=0.5,linestyle='--',color='red')

    ymax=n+0.5
    if beam:
        ymax+=1
        ax.plot(xb,yb+n+1,linewidth=0.25,color='black')
        ax.annotate(f'beam', xy=(xmin,n+1), xytext=(-5, 0), textcoords="offset points", ha='right', va='center')


    
    ax.set_xlim(xmin,xmax)
    xmajor=5
    xminor=1
    xmajor,xminor=tick_stride(xmin,xmax,base=1,prec=2)
    
    ax.xaxis.set_major_locator(MultipleLocator(xmajor))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
    ax.xaxis.set_minor_locator(MultipleLocator(xminor))
    ax.xaxis.grid(b=True, which="minor", **color)
    ax.xaxis.grid(b=True, which="major", **color)
    ax.set_xlabel('Time (sec)',fontdict_axis)

    ax.set_ylim(ymin,ymax)
    ax.tick_params(labelleft=False,labelsize=axis_tick_size) 
    ax.yaxis.set_major_locator(MultipleLocator(1))
    ax.yaxis.grid(b=True, which="minor", **color)
    ax.yaxis.grid(b=True, which="major", **color)
    ax.set_title(f'Record Section',fontdict=fontdict_title,loc='left')

    fig.savefig(outfig,format='png',bbox_inches="tight")



def plot_twosta(f,Cxy,id1,id2,m):
    fig=plt.figure(figsize=(6,6),dpi=200)
    gs=fig.add_gridspec(1,1)
    ax=fig.add_subplot(gs[0])
    ax.scatter(f,Cxy,marker='o',s=10,linewidth=.1,c='red',edgecolor='black')
    ax.plot(f,Cxy,linewidth=2,c='midnightblue',alpha=.6)

    ax.set_xlim(f[0],f[-1])
    xmajor=5
    xminor=1
    ax.xaxis.set_major_locator(MultipleLocator(xmajor))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
    ax.xaxis.set_minor_locator(MultipleLocator(xminor))
    ax.xaxis.grid(b=True, which="minor", **color)
    ax.xaxis.grid(b=True, which="major", **color)
    ax.set_xlabel('Frequency(Hz')

    ax.set_ylim(0.0,1.)
#    ax.set_yscale('log')
    xmajor=0.25
    xminor=0.05
    ax.yaxis.set_major_locator(MultipleLocator(xmajor))
#    formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))
#    ax.yaxis.set_major_formatter(formatter)

    ax.yaxis.set_major_formatter(FormatStrFormatter('%0.1g'))
    ax.yaxis.set_minor_locator(MultipleLocator(xminor))
    ax.yaxis.grid(b=True, which="minor", **color)
    ax.yaxis.grid(b=True, which="major", **color)
    ax.set_ylabel('Coherence')
    ax.set_title(f'Coherence between {id1} and {id2} ({m:0.2} km)',fontdict={'fontsize':8},loc='left')

    outfile=f'cxy_{id1}_{id2}.png'
    plt.savefig(outfile,bbox_inches='tight')
    plt.close()

def make_fig(ans,fcs,outfile,fls=None,fus=None,domean=True):
    _name='make_fig'

    for fc,fl,fu in zip(fcs,fls,fus):
        fig = plt.figure(figsize=(8, 6),dpi=200)

        gs=fig.add_gridspec(1,1)
        ax=fig.add_subplot(gs[0])

        fc_actual=plot_coherdist(ax,ans,fc,fl=fl,fu=fu,domean=domean)

        gs.update(wspace=0.05, hspace=0.20)


        outpng=f'{outfile}{fc:07.4f}.png'
        plt.savefig(outpng,bbox_inches='tight')
        plt.close()


def plot_coherdist(ax,ans,freq,fl=None,fu=None,domean=None):
    _name=f'{__name__}.coherdist'
    color={"color": "0.9"}
    colors = [(1,1,1), (0, 0, 1), (0, 1, 0), (1, 0, 0)]
    colors = [(0, 0, 1), (0, 1, 0), (1, 0, 0)]
    cmap = LinearSegmentedColormap.from_list('my_colors', colors, N=6)


    dists=[]
    Cxys=[]
    azs=[]
    freqs=ans[0][4]
    if domean:
        _idx0=find_nearest(freqs,fl)
        _idx1=find_nearest(freqs,fu)
    _idx=find_nearest(freqs,freq)

    
    for i in ans:
        dists.append(i[0])
        azs.append(i[1])
        if domean:
            mn=np.mean(i[5][_idx0:_idx1])
            Cxys.append(mn)
        else:
            Cxys.append(i[5][_idx])
#ax1.scatter(x, data, c=wts, alpha=0.6, edgecolors='none', cmap=cmap)

    scat=ax.scatter(dists,Cxys,c=azs,alpha=0.6,linewidth=0.35,marker='o',s=50,edgecolor='black',cmap=cmap)


    # xaxis stuff
    xmajor,xminor=tick_stride(np.min(dists),np.max(dists),base=1,prec=2)
    ax.set_xlim(0,np.max(dists)*1.05)
    ax.xaxis.set_major_locator(MultipleLocator(xmajor))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%0.1f'))
    ax.xaxis.set_minor_locator(MultipleLocator(xminor))
    ax.xaxis.grid(b=True, which="minor", **color)
    ax.xaxis.grid(b=True, which="major", **color)
    ax.set_xlabel('Interstation Distance (km)')
#        ax.tick_params(labelbottom=False)    
#    
#        # yaxis stuff
    ax.set_ylim(0,1.05)
    ymajor,yminor=tick_stride(0,1,base=.1,prec=2)
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
    cbar.set_label(f'Interstation Azimuth $^\deg$')
    ax1.invert_yaxis()
    ax1.yaxis.set_major_locator(MultipleLocator(30))
#
#
#        # Title
    ax.set_title(f'Signal coherence at center frequency: {freq:0.4f} Hz',fontdict={'fontsize':8},loc='center')
#        ax.set_title(f'Ref:{true_start}',fontdict={'fontsize':8},loc='right')

    return  freqs[_idx]


def tick_stride(xmin,xmax,prec=2, base=1):
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

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx


def get_idx(arr,b,e):
    _name='get_idx'
    if b > arr[0] and e < arr[-1]:
        _ind=np.nonzero((arr >= b) & (arr <= e))[0]
    elif b > arr[0] and e > arr[-1]:
        _ind=np.nonzero(arr >= b)[0];
    elif b < arr[0] and e < arr[-1]:
        _ind=np.nonzero(arr <= e)[0];
    else:
        _ind=np.nonzero(arr > -1)[0]
    return _ind


