'''
COPYRIGHT 2010 Alexander Crosby

This file is part of SCI-WMS.

    SCI-WMS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    SCI-WMS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with SCI-WMS.  If not, see <http://www.gnu.org/licenses/>.
'''

import netCDF4, sys, os
from sklearn.neighbors import BallTree, KDTree
from datetime import datetime
import numpy as np
try:
    import cPickle as pickle
except:
    import Pickle as pickle

def build_from_nc(filename):
    timer = datetime.now()
    nc = netCDF4.Dataset(filename)
    filename = filename[:-3]
    if nc.grid == 'cgrid':
        lat = nc.variables['lat'][:].flatten()
        lon = nc.variables['lon'][:].flatten()
        #print nc.variables['lon'].shape[0]*nc.variables['lon'].shape[1]
        #nc.close()
        #print np.asarray(zip(lon, lat)).shape
        try:
            latb = np.ma.getmaskarray(lat)
            lonb = np.ma.getmaskarray(lon)
            lat = np.ma.getdata(lat)
            lon = np.ma.getdata(lon)
            lat[latb] = -9999
            lon[lonb] = -9999
        except:
            pass
        #print nc.variables['lon'].shape[0]*nc.variables['lon'].shape[1]
        tree = BallTree(np.asarray(zip(lon, lat)), leaf_size=1000000, metric='euclidean')
        with open(filename+"_nodes.tree", 'w') as f:
            pickle.dump(tree, f)
    else:
        lat = nc.variables['lat'][:]
        lon = nc.variables['lon'][:]
        latc = nc.variables['latc'][:]
        lonc = nc.variables['lonc'][:]
        nv = nc.variables['nv'][:] # (3, long)
        nc.close()
        print np.asarray(zip(lon, lat)).shape
        node_tree = BallTree(np.asarray(zip(lon, lat)), leaf_size=1000000, metric='euclidean')
        with open(filename+"_nodes.tree", 'w') as f:
            pickle.dump(node_tree, f)
        del lon, lat, node_tree
        cell_tree = BallTree(np.asarray(zip(lonc, latc)), leaf_size=1000000, metric='euclidean')
        with open(filename+"_cells.tree", 'w') as f:
            pickle.dump(cell_tree, f)
        del lonc, latc, cell_tree

if __name__ == "__main__":
    filename = sys.argv[1]
    build_from_nc(filename)
