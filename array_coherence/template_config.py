def template_config():
    template=f'''

# PARM FILE----CUT HERE -----
[TIME]
# SAC reference time field. Time windows relative to this value
# Typical choices are o,b,a,t0,t1. Case insensitive.
reffeild=o
t0=0.0
t1=2.5

[SIG_PROC]
detrend=linear
taper=0.05
hp=0.0
lp=0.0
npoles=2
passes=1



[PLOT]
do_wf=True
outfile=S12_template.png

# ------------CUT HERE ------
    ''' 
    return template

