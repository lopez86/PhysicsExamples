import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.basemap import Basemap
# Coordinates in standard spherical coords
# Z-axis is north pole

# Hanford, WA 
hanford = np.array([0.7578,2.0837])
# Livingston, LA 
livingston = np.array([1.0386,1.5839])
# Near Pisa, Tuscany
virgo = np.array([0.8093,0.1833])

earthRadius = 6371e3 # m
speedOfLight = 3e8 # m/s

def Likelihood(time,uncertainty):
    # We're going to take the difference

    def to_cartesian(th,ph):
        return np.array([
                         np.sin(th)*np.cos(ph),
                         np.sin(th)*np.sin(ph),
                         np.cos(th)
                       ])

    uncertainty = np.sqrt(2) * uncertainty
    hf = to_cartesian(hanford[0],hanford[1])
    ls = to_cartesian(livingston[0],livingston[1])

    npts = 200
    phi = np.linspace(-np.pi,np.pi,npts)
    theta = np.linspace(0,np.pi,npts)
    grid_theta,grid_phi = np.meshgrid(theta,phi)

    likelihood = np.zeros([npts,npts])

    for i in range(npts):
        for j in range(npts):
            u = to_cartesian(grid_theta[i,j],grid_phi[i,j])
            hfCos = u.dot(hf)
            lCos = u.dot(ls)
            dx = earthRadius * (hfCos-lCos)
            dt = dx / speedOfLight            
            #likelihood[i,j] = 1./ np.sqrt(2*np.pi*uncertainty**2)  \
            likelihood[i,j] = np.exp(-0.5*( (dt-time)/uncertainty) **2)
            #print(i,j,dt,t,likelihood[i,j])

    return grid_theta,grid_phi,likelihood

def main(time,uncertainty):
    th,ph,ll = Likelihood(time,uncertainty)
    
    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.pcolor(th * 180./np.pi,ph * 180./np.pi,ll,cmap='gist_heat',
               vmin=0,vmax=1)
    #plt.imshow(ll, interpolation='nearest', cmap='gist_heat'
    #           extent=[0, 180, -180, 180])

    ax.set_xlabel('Theta [deg]')
    ax.set_ylabel('Phi [deg]')

    fig = plt.figure(2)
    ax = fig.add_subplot(111)
    bmap = Basemap(projection='hammer',lon_0=0)
    bmap.drawmeridians(np.arange(0,360,60),color='0.5')
    bmap.drawparallels(np.arange(-90,90,30),color='0.5')
    lats = 90 - 180./np.pi * th
    lons = ph * 180./np.pi
    x,y = bmap(lons,lats)
    
    #print(x)
    #print(y)
    bmap.pcolor(x,y,ll,cmap='gist_heat',vmin=0,vmax=1)
    cb = bmap.colorbar(location='right',label='Likelihood')
    plt.title(r'Hammer Projection, $\Delta t = $%.2e $\pm$ %0.2e'%(time,uncertainty))

    cb.set_ticks([0,0.2,0.4,0.6,0.8,1.0])

    plt.show()

if __name__=='__main__':
    argc = len(sys.argv)
    t = float(sys.argv[1])
    unc = 1e-6
    try:
        unc = float(sys.argv[2])
    except:
        unc = 1e-6
    main(t,unc)
