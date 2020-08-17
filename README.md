## array_coherence ##
Compute the Coherence Squared between a suite of SAC files.  


### Purpose/Scope ###
This script computes the coherence of an event signal between a set of 
seismic array elements. The data are expected to be in SAC format, and several 
header vales need to be set. The values are, at minimum,
*stla,stlo,evla,evlo*, and some time reference *(e.g. a,t0,o)*.  Although
the header value (*b*) may also be used. 

This version only uses scipy.coherence, which uses Welch's method to estimate
the spectra. The multitaper estimation method is currently turned off. 

The script makes a suite of plots:  
1. A record section of each time series sorted by distance, along with the array beam  
2. The Two-Station coherence as a function of frequency for each station pair  
3. The Coherence as a function of interstation distance at specific  
   frequency values.

## Install ##

Clone source package  
`git clone http://github.com/flyrok/array_coherence`

Install with pip after download  
`pip install .`

Install in editable mode  
`pip install -e .`


## Python Dependencies ##

python>=3.6 (script uses f-strings)  
obspy (https://github.com/obspy/obspy/wiki)  


## Usage/Examples ##
The main driver is a INI configuration file. Use the `-h` option to
print an example. It is recommended to run the script with the `-v` to
understand some of the output. It is very verbose. Turn on debugging 
with `-vv` if things go wrong. 

To see help:  
`array_coherence -h`    

To see version:  
`array_coherence --version`    

To run it:  
`array_coherence -f *.sac -i example.ini

