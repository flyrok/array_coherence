def template_config():
    template=f'''

# PARM FILE----CUT HERE -----

[TIME]
## SAC reference time field. ##
#   *reffield*
#   All times are relative to the value stored in this sac header. Typical
#   SAC header values to chose are o,a,t0,t1. Case insensitive.  These need
#   to be set in the SAC header prior to running this code.  SAC header
#   values a, t0, are referenced to *o*, which is usually the zero time
#   in the sac file, which is typicaly the event origin time. By
#   convention, the following is usually true
#   a->P_time, t0->S_time, t1->some_other_phase, o=origin_time
#   (str)
reffield=a
## Coherence and xcorr time window ##
#   seconds are referenced to reffield, as:
#   *start_window=reffield_time + startsec*
#   *end_window=start_window+duration*
# (floats)
startsec=-0.5
duration=4.0
## plot window ##
#   Seconds are referenced to reffield
#   These values need to define a larger time window than *startsec* &
#   *duration( paramters above, this is for plotting
#   (floats)
plot_startsec=-5.
plot_dursec=25.0
# Set to *True* to align time series by a x-corr timeshift.
# By default, each time series is x-correlated to the closest station
# Otherwise, set to *False* to align time series to *reffield*
# (boolean)
xcorr=False

[COHERENCE]
## The power spectra are estimated using Welch's method. ##
#   As such, the following paramters relate to Welch's method.
#   The following parameters must be set to either None or an integer
#   (ints)
#   Length of each segment in Welch's (in samples)
nperseg=32
#   Overlap number (in samples) between segments. Defaults to 50% overlap
noverlap=8
#   Number of fft points. Data are zero padded if NFFT longer than nperseg
#   if set to None, nfft=2**n
nfft=128

[SIG_PROC]
## Signal processing ##
#    By default, A linear trend is removed each time series, a taper is applied, 
#    and optionally a lowpass, highpass, or bandpass is applied.
#    (str)
detrend=linear
#    Raised Cosine taper. Percent taper 
#    (float)
taper=0.05
#    high pass, set to 0.0 to turn off
#    (float)
hp=1.0
#    low pass, set to 0.0 to turn off
#    (float)
lp=0.0
#    number of filter poles
npoles=3
#    number of passes
#    int 0=minimum phase, 1=zerophase
passes=0

[PLOT]
#    output record section file name
wffig=recsection.png
#    output basename for coherence vs distance plots
cdist_base=SPE-5_Coh_v_dist_
#    The follow parameters define the scope of the 
#    coherence versus distance plots.
#    fmin & fmax are the range of center frequencies for the plots
#    bw defines the octave freq spacing between fmin & fmax
#    bw values could be 1/12,1/8,1/4,1/3,1/2,1 etc
#    The larger the bw, the larger the averaging
#    note: 
#    (floats
fmin=1.0
fmax=20.0
bw=0.5
#     If this is set to True, then coherence values are 
#     averaged over each octave banded, Otherwise no averaging
domean=True

# ------------CUT HERE ------
    ''' 
    return template

