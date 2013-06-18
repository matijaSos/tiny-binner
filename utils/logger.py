import logging
import __builtin__

class Logger(object):

    def __init__(self, log_filename):

        __builtin__.log = logging.getLogger()
        ch = logging.StreamHandler(file(log_filename, 'w'))
        ch.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
        ch.setFormatter(formatter)
        log.addHandler(ch) 
        log.setLevel(logging.DEBUG)
