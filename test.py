import multiprocessing as mp

def square(x):
    return x*x

if __name__ == '__main__':
    p = mp.Pool(5)
    vals = range(100000)
    ret = p.map(square, vals)
    print (ret)
    