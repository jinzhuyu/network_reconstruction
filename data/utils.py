import pickle
import numpy as np

def load_synthetic_data(dyn, net):
    '''
    Load synthetic networks
    '''
    fnm = 'net-{}.pkl'.format(net)
    with open(fnm,'rb') as file:
        dat=pickle.load(file)
    G, ks, A_real, A_hh = dat['G'], dat['ks'], dat['A_real'], dat['A_hh']

    fnm = 'xs-{0}-{1}.pkl'.format(dyn,net)
    with open(fnm,'rb') as file:
        xs = pickle.load(file)

    fs,gs=interactions(dyn, xs, ks)
    return G, xs, ks, A_real, A_hh, fs, gs


def interactions(dyn,xs,ks=None):
    '''
    Obtain the interaction terms (self or mutual interactions)

    Parameters:
    ----------
    dyn:        dynamics
    xs:         steady states
    ks:         degrees

    Returns:f
    --------
    self and mutual interaction terms
    '''
    N = len(xs)
    if dyn == 'ecology':
        B, C, D, E, H, K = 0.1, 1, 5, 0.9, 0.1, 5
        fs = B + xs * (1 - xs / K) * (xs / C - 1)
        gs = np.array([[xs[i] * xs[j] / (D + E * xs[i] + H * xs[j]) for j in range(N)] for i in range(N)])
    elif dyn == 'regulatory':
        B, f, h = 1, 1, 2
        fs = -B * xs ** f
        gs = np.repeat([[xs[j] ** h / (xs[j] ** h + 1) for j in range(N)]], [N], axis=0)
    elif dyn == 'epidemic':
        R, B = 1, 2
        fs = -B * xs
        gs = np.array([[R * (1 - xs[i]) * xs[j] for j in range(N)] for i in range(N)])
    elif dyn == 'pagerank':
        alpha = 0.85
        fs = (1-alpha)/N - xs
        gs = alpha*np.array([[xs[j]/ks[j] for j in range(N)] for i in range(N)])

    for i in range(N):
        gs[i, i] = 0
    return fs, gs


def main():

    import os
    os.chdir('C:/code/network_reconstruction/data')
    
    # os.chdir('/mnt/c/code/network_reconstruction/data')
    
    
    dyn='ecology'
    real = False
    
    net = 'BA-n{0}-m{1}-{2}'.format(100,4,99)
    
    G, xs_real, ks, A_real, A_hh, fs, gs = load_synthetic_data(dyn, net)
    
    # save data
    import numpy as np
    
    # np.savetxt('ks.txt', ks, delimiter = ',')
    # np.savetxt('fs.txt', ks, delimiter = ',')
    # np.savetxt('gs.txt', ks, delimiter = ',')
    
    np.save('gs', gs)
    np.save('ks', ks)
    np.save('fs', fs)
 
    
#####  
if __name__ == "__main__":  
    main()
        



