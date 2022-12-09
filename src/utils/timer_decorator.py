import time 

def timer(func):
    "returns output and execution time of 'func'"
    def wrap_func(*args, **kwargs):
        t1 = time.time()
        result = func(*args, **kwargs)
        t2 = time.time()
        print(f'Function {func.__name__!r} executed in {(t2-t1):.4f}s')
        return result, round(t2-t1,2)
    return wrap_func