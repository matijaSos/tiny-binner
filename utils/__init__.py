import time

def enum(**enums):
    return type('Enum', (), enums)

def timeit(function, *args):
    start = time.time()
    ret_value = function(*args)
    stop = time.time()
    print stop-start
    return ret_value
