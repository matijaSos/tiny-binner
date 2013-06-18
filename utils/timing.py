# -*- coding: utf-8 -*-
'''
Created on Apr 27, 2013

@author: marin
'''

from datetime import datetime
import dateutil.relativedelta as relativedelta

def start():
    return datetime.now()

def end(start):
    return relativedelta.relativedelta(datetime.now(), start)

def humanize(download_delta):
    return '%02d:%02d:%02d' % (download_delta.hours, download_delta.minutes, download_delta.seconds)