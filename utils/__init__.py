import time

def enum(**enums):
    return type('Enum', (), enums)

def timeit(function, *args):
    start = time.time()
    function(*args)
    stop = time.time()
    print stop-start
