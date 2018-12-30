import numpy as np
import h5py
import sys
import getopt

def compress_dataset(fname):
    with h5py.File(fname, 'r') as f:
        stat = f["status"][:]
        esc_frac = np.sum(stat==1, axis=-1)/float(stat.shape[-1])
        
        non_valid = np.sum(np.sum(f["position_start"][:],axis=0),axis=-1)==0
        esc_frac[non_valid] = np.nan

#        esc_loc = np.mean(np.sqrt(np.sum(f["position_end"]**2,axis=0)),axis=-1)
#        esc_loc[non_valid] = np.nan

        esc_vel = np.mean(np.sqrt(np.sum(f["velocity_end"][:]**2,axis=0)),axis=-1)
        esc_vel[non_valid] = np.nan

    lin = np.linspace(-4,4,esc_frac.shape[0])*3390.0
    xmesh, ymesh, zmesh = np.meshgrid(lin,lin,lin, indexing='ij')

    with h5py.File(fname[:-3]+"_compressed.h5", 'w') as f:
        f.create_dataset("escape_fraction", data=esc_frac)
        f.create_dataset("escape_velocity", data=esc_vel)
        f.create_dataset("xmesh", data=xmesh) 
        f.create_dataset("ymesh", data=ymesh) 
        f.create_dataset("zmesh", data=zmesh) 

        f.attrs.create('radius', 3390.0)


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
