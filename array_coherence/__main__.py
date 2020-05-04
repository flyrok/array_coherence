#!/usr/bin/env python3

import argparse
from pathlib import Path
from datetime import datetime
import sys
from array_coherence import array_coherence

here = Path(__file__).resolve().parent
exec(open(here / "version.py").read())

def main():
    '''
    Routine to compute/plot array coherence
    '''

    parser = argparse.ArgumentParser(prog=progname,
            formatter_class=CustomFormatter,
            description= '''
            ''',
            epilog="")

    parser.add_argument("-r","--reffield", type=str,
        required=False, help="If set, and SAC input, then this defines zero time. Startsec/endsec are relative to this.")

    parser.add_argument("-s","--startsec", type=float,
        required=False, help="Number of seconds relative to reffield start to computation. Defaults to 0, aka 1st sample in input.")

    parser.add_argument("-e","--endsec", type=float,
        required=False, help="If not set, then takes the last sample of the input trace. Otherwise its the duration")

    parser.add_argument("-f","--wfile",type=str, nargs='*',
        required=True, help="Input seismic data file. SAC, for now.")

    parser.add_argument("--lp", type=float,default=None,
        required=False, help="Lowpass filter.")

    parser.add_argument("--hp", type=float,default=None,
        required=False, help="Highpass filter.")

    parser.add_argument("--npoles", type=int, default=None,
        required=False, help="Npoles of filter.")

    parser.add_argument("-o","--outpng", type=str, default='tmp.png',
        required=False, help="output png file")

    parser.add_argument("-v", "--verbose", action="count",default=0,
        help="increase debug spewage spewage (e.g. -v, -vv, -vvv)")

    parser.add_argument('--version', action='version',
        version='%(prog)s {version}'.format(version=__version__))

    args = parser.parse_args()
    startsec=args.startsec
    endsec=args.endsec
    reffield=args.reffield
    wfile=args.wfile
    outpng=args.outpng
    lp=args.lp
    hp=args.hp
    npoles=args.npoles
    debug=args.verbose

    # Do it all 
    array_coherence(wfile,reffield=reffield,startsec=startsec,endsec=endsec,hp=hp,lp=lp,npoles=npoles,outfile=outpng,debug=debug)


class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
    '''
    re-class ArgDefaults formatter to also print things pretty. Helps printing out the config file example
    '''
    def _get_help_string(self, action):
        help = action.help
        if '%(default)' not in action.help:
            if action.default is not argparse.SUPPRESS:
                defaulting_nargs = [argparse.OPTIONAL, argparse.ZERO_OR_MORE]
                if action.option_strings or action.nargs in defaulting_nargs:
                    if type(action.default) == type(sys.stdin):
                        print( action.default.name)
                        help += ' (default: ' + str(action.default.name) + ')'
                    else:
                        help += ' (default: %(default)s)'
        return help

if __name__ == '__main__':
    main()

