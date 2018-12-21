import numpy as np
import h5py
import sys
import getopt

def compress_dataset(fname):
    with h5py.File(fname) as f:
        for key in list(f.keys()):
            N =  np.argmax(np.all(f[key][:] == 0, axis=0))
            dat = f[key][:, :N]
            del f[key]
            f[key] =dat


def main(argv):
    try:
        opts, args = getopt.getopt(argv,"i:",
                                    ["infile="])
    except getopt.GetoptError:
        print(getopt.GetoptError())
        print('error')
        return
    
    infile = None
    for opt, arg in opts:
        if opt in ("-i", "--infile"):
            infile = arg

    compress_dataset(infile)

if __name__ == "__main__":
    main(sys.argv[1:])
