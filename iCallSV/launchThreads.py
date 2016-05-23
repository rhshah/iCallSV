"""
Created on December 21, 2015
Description: This module will be launching functions as threads
@author: Ronak H Shah
"""
"""
::Inputs::
threadID: Id for the thread
name: name for the thread
counter: a counter to keep track of the thread
"""
"""
Taken from:
http://www.tutorialspoint.com/python/python_multithreading.htm

"""

import logging
import threading
import time

logger = logging.getLogger(__name__)

class myThread (threading.Thread):

    def __init__(self, threadID, name, counter):
        threading.Thread.__init__(self)
        self.threadID = threadID
        self.name = name
        self.counter = counter

    def run(self):
        
        logger.info("Starting %s", self.name)
        # Get lock to synchronize threads
        threadLock.acquire()
        print_time(self.name, self.counter, 3)
        # Free lock to release next thread
        threadLock.release()


def print_time(threadName, delay, counter):
    while counter:
        time.sleep(delay)
        print "%s: %s" % (threadName, time.ctime(time.time()))
        counter -= 1

threadLock = threading.Lock()
threads = []

# Create new threads
thread1 = myThread(1, "Thread-1", 1)
thread2 = myThread(2, "Thread-2", 2)

# Start new Threads
thread1.start()
thread2.start()

# Add threads to thread list
threads.append(thread1)
threads.append(thread2)

# Wait for all threads to complete
for t in threads:
    t.join()
print "Exiting Main Thread"
