'''
Created on December 21, 2015
Description: This module will be launching functions as threads
@author: Ronak H Shah
'''
'''
::Inputs::
threadID: Id for the thread
name: name for the thread
counter: a counter to keep track of the thread
'''
'''
Taken from:
http://www.tutorialspoint.com/python/python_multithreading.htm

'''

import logging
import threading
import time
class myThread (threading.Thread):
    def __init__(self, threadID, name, counter):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.counter = counter
    def run(self):
        logging.info( "Starting %s" , self.name )
        # Get lock to synchronize threads
        threadLock.acquire()
        print_time(self.name, self.counter, 3)
        # Free lock to release next thread
        threadLock.release()