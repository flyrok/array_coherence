def template_config():
    template=f'''

# PARM FILE----CUT HERE -----
[TIME]
## SAC reference time field. ##
#   *reffield*
#   All times are relative to the value stored in this sac header i.
#   Typical SAC header values to chose are o,a,t0,t1. Case insensitive.
#   These need to be set in the SAC header prior to running this code.
#   SAC header values a, t0, are reference to *o* which is usually the 
#   zero time in the sac file, which is typicaly the event origin time.
#   For more reference,
#   a->P_time, t0->S_time, t1->some_other_phase, o=origin_time,
#   (str)
reffield=a
## Coherence and xcorr time window ##
#   seconds are referenced to reffield, as
#   start_window=reffield_time + startsec, 
#   end_window=start_window+duration
# (floats)
startsec=-1.0
duration=5.0
## plot window ##
#   seconds are reference to reffield
#   These values need to define a larger time window
#   than startsec & duration paramters above, this is for plotting
#   floats
plot_startsec=-10.
plot_dursec=40.0

[COHERENCE]
## The power spectra are estimated using Welch's method. ##
#   As such, the paramters relate to Welch's method.
#   The follower parameters must be set to either None or an integer
#   Length of each segment in Welch's (in samples)
nperseg=64
#   Overlap number (in samples) between segments. Defaults to 50% overlap
noverlap=32
#   Number of fft points. Data are zero padded if longer than nperseg
nfft=64

[SIG_PROC]
detrend=linear
# Raised Cosine taper. Percent taper 
taper=0.05
# high pass
hp=2.5
# low pass
lp=0
# number of filter poles
npoles=3
# number of passes
# int 0=minimum phase, 1=zerophase
passes=0

[PLOT]
# output record section file name
wffig=recsection.png
# output basename for coherence vs distance plots
cdist_base=DAG-2_Coh_v_dist
# The follow parameters define the scope of the 
# coherence versus distance plots.
# fmin & fmax are the min/max center frequencies for the plots
# bw defines the octave freq spacing between fmin & fmax
# bw values could be 1/8, 1/4,1/3,1/4, etc
fmin=1.0
fmax=20.0
bw=0.333
# If this is set to True, then 
domean=True

# ------------CUT HERE ------
    ''' 
    return template

