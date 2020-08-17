#!/usr/bin/env python3

import argparse
from pathlib import Path
from datetime import datetime
import sys
from array_coherence import array_coherence
from array_coherence import template_config
from array_coherence import progname,__version__

def main():
    '''
    Routine to compute/plot array coherence.
    '''

    parser = argparse.ArgumentParser(prog=progname,
            formatter_class=CustomFormatter,
            description= '''
            ''',
            epilog=template_config())

    parser.add_argument("-f","--wfiles",type=str, nargs='*',
        required=True, help="Input seismic data files. SAC, for now.")

    parser.add_argument("-i","--ini",type=str, default=None,
        required=True, help="Input initialization file. See -h for example.")


    parser.add_argument("-v", "--verbose", action="count",default=0,
        help="increase debug spewage spewage (e.g. -v, -vv, -vvv)")

    parser.add_argument('--version', action='version',
        version='%(prog)s {version}'.format(version=__version__))

    args = parser.parse_args()
    sac_files=args.wfiles
    ini_file=args.ini
    debug=args.verbose

    # Do it all 
    arr_obj=array_coherence(sac_files=sac_files,ini_file=ini_file,debug=debug)
    arr_obj.run_coherence_scipy()


# This class overwrite makes help messages prettier
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

