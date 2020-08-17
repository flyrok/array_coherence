## array_coherence ##
Compute the Coherence Squared between a suite of SAC files.  




### Purpose/Scope ###
This script computes the coherence of an between signal between a set of 
seismic array elements. The data is expected to be in SAC format, and several 
header vales are required to be set. The values are, at minimum,
stla,stlo,evla,evlo

The script makes a series of plots:  
1. record section sorted by distance of all the data  
2. The Two-Station coherence as a function of frequency  
3. The Coherence as a function of interstation distance as specific  
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
The main driver is a INI configuration files. Use the `-h` option to
print an example. It is recommend to run the script with the `-v` to
understand what it is doing. It is very verbose. Turn on debugging `-vv`
if things go wrong to get hints of the possible issue.

To see help:  
`array_coherence --help`    

To see version:  
`array_coherence --version`    

To run it:  
`array_coherence -f *.sac -i example.ini

