import logging
import sys
def setup_log(debug):

    ''' Helper function to set up stdout logging
        at the right debug level
    '''
    # INFO,DEBUG,WARN
    if debug == 0:
        loglevel="WARN"
    elif debug == 1:
        loglevel="INFO"
    elif debug > 1:
        loglevel="DEBUG"

    log=logging.getLogger()
    stream=sys.stdout,
    logging.basicConfig(stream=sys.stdout,datefmt="%Y-%jT%H:%M:%S",
        format="%(asctime)s-%(levelname)s %(message)s")
    log.setLevel(loglevel)
    if debug:
        log.info(f'{__name__} loglevel set to {loglevel}:{log.level}')
    logging.getLogger('matplotlib.font_manager').disabled = True
    logging.getLogger('matplotlib.colors').disabled = True
    logging.getLogger('matplotlib.ticker').disabled = True
    return log

