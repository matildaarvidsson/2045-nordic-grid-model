import geopandas as gpd
import matplotlib.pyplot as plt
from fiona.drvsupport import supported_drivers
from pprint import pprint
import contextily as cx


supported_drivers['KML'] = 'rw'

file = "output\2022-11-24_084013\test.kml"

df = gpd.read_file(file, driver='KML')
df_wm = df.to_crs(epsg=3857)

ax = df_wm.plot(figsize=(30, 30), alpha=0.5, edgecolor='k')
cx.add_basemap(ax, zoom=6)
plt.show()