import matplotlib.pyplot as plt
import matplotlib.colors as cols
import matplotlib.ticker as mticker
import matplotlib.patches as mpatches
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cartopy
from cartopy.util import add_cyclic_point
import netCDF4
from pyproj import Transformer, transform, CRS
import matplotlib.tri as tri
from matplotlib.colors import Normalize, TwoSlopeNorm
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import region_config


def dist(i1, i2, xCell, yCell):  # helper distance fn
    dist = ((xCell[i1]-xCell[i2])**2 + (yCell[i1]-yCell[i2])**2)**0.5
    return dist

pii = 3.14159

fname = '/lcrc/group/acme/data/inputdata/ocn/mpas-o/oEC60to30v3wLI/oEC60to30v3wLI60lev.171031.nc'
fmesh = netCDF4.Dataset(fname)

lonCell = fmesh.variables['lonCell'][:]/pii*180.0
latCell = fmesh.variables['latCell'][:]/pii*180.0
bottomDepth = fmesh.variables['bottomDepth'][:]
idx = np.nonzero( (latCell<-65.0) )[0]; size=1; sz=(14,10); #entire SO
data = bottomDepth[idx]


polar_stereo = '+proj=stere +lat_ts=-71.0 +lat_0=-90 +lon_0=0.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=WGS84'  # Note: BEDMAP2 elevations use EIGEN-GL04C geoid


latmin = 50
projection = cartopy.crs.SouthPolarStereo(true_scale_latitude=-71.0)
#projection = cartopy.crs.CRS(polar_stereo)
#data_crs = cartopy.crs.SouthPolarStereo()
data_crs = cartopy.crs.SouthPolarStereo(true_scale_latitude=-71.0)
#data_crs = cartopy.crs.PlateCarree()
fig = plt.figure(figsize=(5, 5), dpi=72)
ax = plt.subplot(111, projection=projection)
ax.coastlines()
ax.set_extent([-1.8e6, -0.4e6, 0.0e6, 2.0e6], crs=projection)


ice_50m = cartopy.feature.NaturalEarthFeature(
              'physical', 'antarctic_ice_shelves_polys', '50m', edgecolor='k',
              facecolor='none', linewidth=0.5)
ax.add_feature(ice_50m, zorder=4)

## BEDMAP2 projection
## make a CRS (coordinate reference system) for projections from Proj string:
crs_in = CRS.from_proj4('+proj=longlat +ellps=WGS84')
crs_out = CRS.from_proj4(polar_stereo)
## define a transformer
ll_to_proj = Transformer.from_crs(crs_in,crs_out)
x,y = ll_to_proj.transform(lonCell[idx], latCell[idx], radians=False)


triang = tri.Triangulation(x, y)
# Maximum distance in m of edges between points.
# Make twice dcEdge to be safe
triMask = np.zeros(len(triang.triangles))
maxDist = 30000.0 * 2.0
for t in range(len(triang.triangles)):
    thisTri = triang.triangles[t, :]
    if dist(thisTri[0], thisTri[1], x, y) > maxDist:
        triMask[t] = True
    if dist(thisTri[1], thisTri[2], x, y) > maxDist:
        triMask[t] = True
    if dist(thisTri[0], thisTri[2], x, y) > maxDist:
        triMask[t] = True
triang.set_mask(triMask)


colormap='YlGnBu'
#colormap='ocean_r'
norm = Normalize(vmin=0.0, vmax=5000.0)
bathymap = ax.tripcolor(triang, data, cmap=colormap, shading='flat', transform=data_crs, norm=norm)
                                      #levels=[0.9999], colors='grey',
                                      #linestyles='solid')
fig.colorbar(bathymap, extend='max')

#transect_coord_lon = [327,335]
#transect_coord_lat = [-76,-72]
transect_name = 'trough_crossshelf_32'
transect_coord_lon, transect_coord_lat = region_config.region_coordbounds[region_config.region_name.index(transect_name)]
transect_x, transect_y = ll_to_proj.transform(transect_coord_lon, transect_coord_lat)
ax.plot(transect_x, transect_y, 'r', linewidth=3, transform=data_crs)
ax.plot(transect_x, transect_y, 'ro', transform=data_crs)

degree_sign = u'\N{DEGREE SIGN}'
#for i in range(2):
#    ax.annotate(f'{transect_coord_lat[i]}{degree_sign} {transect_coord_lon[i]}{degree_sign}', [transect_x[i], transect_y[i]], transform=data_crs)




#gl = ax.gridlines(crs=cartopy.crs.PlateCarree(), draw_labels=True,
#            linewidth=1, color='gray', alpha=0.5, linestyle='--')
#gl.xlabels_top = False
#gl.ylabels_left = False
#gl.xlines = False
##gl.xlocator = mticker.FixedLocator([-180, -45, 0, 45, 180])
#gl.xformatter = LONGITUDE_FORMATTER
#gl.yformatter = LATITUDE_FORMATTER








#cntr1 = ax.scatter(lonCell[idx], latCell[idx], s=8, c=bottomDepth[idx], transform=data_crs)
#plt.axis('equal')
plt.savefig('transect_location_map.pdf')
plt.show()

## northern limit for subplots
#ms=1 
##xcell = fmesh.variables['xCell'][:]
##ycell = fmesh.variables['yCell'][:]
##latcell = fmesh.variables['latCell'][:]
#idx_scope = pick_from_region(region=scope_name, run=run, plot_map=False)
#zmax_scope = np.multiply(-1.,fmesh.variables['bottomDepth'][idx_scope])
#icemask_scope  = fmesh.variables['landIceMask'][0,idx_scope]
##zice_scope     = fmesh.variables['landIceDraft'][idx_scope]
#ycell_scope = ycell[idx_scope]
#xcell_scope = xcell[idx_scope]
#
#fig = plt.figure()
#ax = fig.add_subplot(111)
#cntr1 = plt.scatter(ycell_scope,xcell_scope,
#                    s=loc_ptsize[region_name.index(scope_name)], c=zmax_scope)
##cntr1 = plt.tricontourf(ycell[idx_scope],xcell[idx_scope],zmax, 20, cmap="viridis")
##plt.plot(ycell[idx_scope],xcell[idx_scope], 'o', color = 'white', 
##         markersize = ms, fillstyle = 'none', alpha = 0.5)
##cntr = plt.tricontour(ycell_scope,xcell_scope,
##                      zice, [-10], colors = 'k',linewidth=lw1)
#plt.plot(ycell_scope[icemask_scope==1], xcell_scope[icemask_scope==1],
#         'o', color = 'white', 
#         markersize = 5*ms, fillstyle = 'none')
#plt.plot(ycell[cellidx], xcell[cellidx],
#         'o', color = 'black', 
#         markersize = ms, fillstyle = 'none')
##plt.plot([y0,y1],[x0,x1], 'k-')
#
#cNorm  = Normalize(vmin=-1*pi, vmax=1*pi)
#scalarMap = cmx.ScalarMappable(norm=cNorm, cmap='jet')
##for i in range(len(idx)):
##    for j in range(6):
##        colorVal = scalarMap.to_rgba(fmesh.variables['angleEdge'][edgesOnCell[i,j]])
##        sc = plt.scatter(fmesh.variables['yEdge'][edgesOnCell[i,j]],
##                    fmesh.variables['xEdge'][edgesOnCell[i,j]],s=ms/2,c=colorVal)
#for k in edgeidx:
#    colorVal = scalarMap.to_rgba(fmesh.variables['angleEdge'][edgesOnCell[i,j]])
#    sc = plt.plot([fmesh.variables['yVertex'][fmesh.variables['verticesOnEdge'][k,0]-1],
#                   fmesh.variables['yVertex'][fmesh.variables['verticesOnEdge'][k,1]-1]],
#                  [fmesh.variables['xVertex'][fmesh.variables['verticesOnEdge'][k,0]-1],
#                   fmesh.variables['xVertex'][fmesh.variables['verticesOnEdge'][k,1]-1]],
#                  'b-',linewidth=lw1)#marker='None',linestyle='-','k')
#        #print(fmesh.variables['yVertex'][verticesOnEdge[i,j,1]])
#        #print(fmesh.variables['xVertex'][verticesOnEdge[i,j,1]])
#        #xv = [fmesh.variables['yVertex'][verticesOnEdge[i,j,0]],
#        #      fmesh.variables['yVertex'][verticesOnEdge[i,j,1]]]
#        #yv = [fmesh.variables['xVertex'][verticesOnEdge[i,j,0]],
#        #      fmesh.variables['xVertex'][verticesOnEdge[i,j,1]]]
#        #plt.plot(xv,yv,'b-',linewidth=0.5)
#plt.plot(ypt[idx], xpt[idx], 'k.',markersize=ms)
##    for i,j in enumerate(idx):
##        plt.plot([yv1[i],yv2[i]],[xv1[i],xv2[i]], 'k-')
##ax.set_xlabel('y (m)')
##ax.set_ylabel('x (m)')
#plt.axis('equal')
#plt.ylim([region_xybounds[region_name.index(scope_name)][0,0],
#          region_xybounds[region_name.index(scope_name)][0,1]])
#plt.xlim([region_xybounds[region_name.index(scope_name)][1,0],
#          region_xybounds[region_name.index(scope_name)][1,1]])
#cbar = plt.colorbar(cntr1)
##cbar = plt.colorbar(sc)
#cbar.set_label(r'Depth (m)')    
#print('save ','bathy_' + transect_name)
#plt.savefig(savepath + 'bathy_' + transect_name)
#plt.close()
    

