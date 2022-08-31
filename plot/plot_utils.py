import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import os
os.environ['PROJ_LIB'] = '/home/Youtong.Zheng/.conda/envs/demo/share/proj'
from mpl_toolkits import basemap
from mpl_toolkits.basemap import Basemap,maskoceans 


def plt2D(fig, ax, x, y, z, z1 = [0], cmap=cm.viridis, levels=np.linspace(0, 40, 10 + 1), levels1 = [35., 40., 45.],
          title = '', ylabel = '', zlabel = 'Pressure (hPa)', xlabel = False, 
          cbar = True, extend = 'both', fontsize = 12, 
          yrange = [1000., 700.], xrange = [-6, 6],
         plotregion = False):
    
    oax = ax.contourf(x,y,z, cmap = cmap, levels = levels, extend = extend)
    
    if np.size(z1) > 1:
        CS = ax.contour(x, y, z1, levels = levels1, colors = 'black', linewidths= 0.4)
        plt.clabel(CS, fontsize=7, inline=1,fmt = '%2.0f')
        
    ax.set_xlim(xrange)
    ax.set_ylim(yrange)
    
    ax.set_ylabel(ylabel, fontsize=fontsize) 
    if xlabel:
        ax.set_xlabel('Distance from ice edge \n ($\degree$ latitude)', fontsize=fontsize) 
    
    ax.set_title(title, fontsize=fontsize)
    
    ax.axvline(x=0., color='black', lw = 0.5, linestyle='--')
    
    if plotregion:
        ax.axvline(x=-2., color='grey', lw = 0.5, linestyle='--')
        ax.axvline(x=2., color='grey', lw = 0.5, linestyle='--')
            
    if cbar:
        cbar = fig.colorbar(oax,  ax=ax, orientation='vertical', shrink=0.8, aspect=15, 
                ticks = [levels.min(), levels.max()])
    
        cbar.set_label(zlabel, fontsize=0.8*fontsize)
    
    return oax

def pltpolar(x, y, var, ice, ax, mycmap = cm.Spectral_r, mylevels = np.linspace(0, 100., num=11),
            mylabel = '', myextend = 'both', myfontsize = 12, hemisphere = 'nh'):
    
    if hemisphere[0:2] == 'nh':
        m = Basemap(projection='npstere',boundinglat=70,lon_0=0,resolution='l', ax=ax)
    else:
        m = Basemap(projection='spstere',boundinglat=-55,lon_0=180,resolution='l', ax=ax)        
    m.drawcoastlines(color='grey')
    
    x, y = np.meshgrid(x, y)
    cs = m.contourf(x,y, var, cmap=mycmap, levels=mylevels, latlon = True, extend = myextend)
    m.contour(x,y, ice,levels=[50], latlon = True, linewidths= 2, colors = 'black')
    
    ax.set_title(mylabel)
    
    if hemisphere == 'nh':
        x,y = m(np.linspace(-15,-15),np.linspace(70,90))
        m.plot(x, y, color = 'red', linewidth=1) 
        x,y = m(np.linspace(55,55),np.linspace(70,90))
        m.plot(x, y, color = 'red', linewidth=1) 
        x,y = m(np.linspace(-15,55),np.linspace(70,70))
        m.plot(x, y, color = 'red', linewidth=1) 
    elif hemisphere == 'nh-chukchi':
        x,y = m(np.linspace(170,170),np.linspace(70,90))
        m.plot(x, y, color = 'red', linewidth=1) 
        x,y = m(np.linspace(-165,-165),np.linspace(70,90))
        m.plot(x, y, color = 'red', linewidth=1) 
    elif hemisphere == 'sh':
        x,y = m(np.linspace(-50,-50),np.linspace(-90,-55))
        m.plot(x, y, color = 'red', linewidth=1) 
        x,y = m(np.linspace(0,0),np.linspace(-90,-55))
        m.plot(x, y, color = 'red', linewidth=1) 
        x,y = m(np.linspace(-50,0),np.linspace(-55,-55))
        m.plot(x, y, color = 'red', linewidth=1) 
    else:
        x,y = m(np.linspace(-180,-180),np.linspace(-90,-55))
        m.plot(x, y, color = 'red', linewidth=1) 
        x,y = m(np.linspace(-130,-130),np.linspace(-90,-55))
        m.plot(x, y, color = 'red', linewidth=1) 
        x,y = m(np.linspace(-180,-130),np.linspace(-55,-55))
        m.plot(x, y, color = 'red', linewidth=1) 
        
    return cs

def plt2l(ax, x, y, y1, y2, title = '', label = ['', '',''], fontsize = 12, 
          yrange = [40., 85.], xrange = [-6, 6], xlabel = True, ylabel = 'Cloud cover (%)',legend = False):
    
    ax.plot(x, y, label = label[0])
    ax.plot(x, y1, label = label[1])
    ax.plot(x, y2, label = label[2], color = 'black')
    
    ax.set_xlim(xrange)
    ax.set_ylim(yrange)
    
    ax.set_ylabel(ylabel, fontsize=fontsize) 
    if xlabel:
        ax.set_xlabel('Distance from ice edge \n ($\degree$ latitude)', fontsize=fontsize) 
    
    ax.set_title(title, fontsize=fontsize)
    ax.axvline(x=0., color='black', lw = 0.5, linestyle='--')
    
    if legend:
        ax.legend(loc = 'best',fontsize=0.8*fontsize)