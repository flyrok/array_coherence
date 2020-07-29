def template_config():
    template=f'''

# PARM FILE----CUT HERE -----
[TIME]
## SAC reference time field. ##
#   Time windows relative to this value
#   Typical header value choices are o,b,a,t0,t1. Case insensitive.
#   str
reffield=a
## Coherence time window ##
#   seconds are referenced to reffield, as
#   t0=reffield_time + startsec, t1=t0+duration
# floats
startsec=-0.25
duration=5
## plot window 
#   seconds are reference to this time
#   floats
plot_startsec=-10.
plot_dursec=40.0

[SIG_PROC]
detrend=linear
taper=0.05
hp=2.0
lp=0.00
npoles=3
# int 0=minimum phase, 1=zerophase
passes=0

[PLOT]
do_wf=True
outfile=DAG-2_wigs.png
# ------------CUT HERE ------
    ''' 
    return template

