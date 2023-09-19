import spherical_poly
import var
import util

# 1st exam
# r = var.R_earth
# p = spherical_poly.s_point([0,90], r=r)
#
# print(p.__dict__)

# 2nd exam
# # print([-123.4112,48.4358].__type__)
# esquimalt = spherical_poly.s_point([-123.4112, 48.4358])
# paris = spherical_poly.s_point([2.3522, 48.8566])
# d = util.great_circle_distance(paris, esquimalt)
# # 8019586.401921335 m
# norm_sp = spherical_poly.s_point(util.great_circle_normal(paris.v,esquimalt.v))
# # lon = -60.31224176777391, lat = -21.860629758043608
# m = util.get_midpoint(paris, esquimalt)
# # lon = -60.9962496812035, lat = 68.13795921981125
# bearing_vec = \
# spherical_poly.s_point(util.bearing_great_circle(esquimalt, paris, bearing = False))
# # lon = 14.403637693988331, lat = 33.306676446329135
# bearing_deg = util.bearing_great_circle(esquimalt, paris,bearing = True, geobearing = True, radians = False)



# intersection

# lon = [0, 45, 0, -45]
# lat = [20, 40, 70, 40]
# dlon = 60
# s_polys = []
# for d in range(int(360/dlon)):
#     points = []
#     for j in range(len(lon)):
#         points.append(spherical_poly.s_point([lon[j] + d * dlon, lat[j]]))
#     s_polys.append(spherical_poly.s_poly(points))
# pole_lon = [d * dlon for d in range(int(360/dlon))]
# pole_lat = [60 for d in range(int(360/dlon))]
# points = []
# for j in range(len(pole_lon)):
#     points.append(spherical_poly.s_point([pole_lon[j], pole_lat[j]]))
# s_polys.append(spherical_poly.s_poly(points))
# # find intersections
# for i in range(len(s_polys)):
#     for j in range(i+1, len(s_polys)):
#         C = spherical_poly.poly_intersection(s_polys[i], s_polys[j])


# 4.1 globe and projections
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
def create_projection(lon = -90, lat = 45, proj = ccrs.Orthographic, figsize = None):
    if figsize is not None:
        fig = plt.figure(figsize = figsize)
    else:
        fig = plt.figure()
        projection_lon = lon
        projection_lat = lat
        projection = proj(projection_lon, projection_lat)
        transform = ccrs.Geodetic()
        ax = fig.add_subplot(1, 1, 1, projection = projection)
        ax.set_global()
        ax.gridlines()
        return(fig, ax, transform)

def add_geography(ax, coast = True, land = True, ocean = False):
    if ocean:
        ax.add_feature(cfeature.OCEAN, zorder=0)
    if land:
        ax.add_feature(cfeature.LAND, zorder=0, edgecolor = [0.6, 0.6, 0.6])
    if coast:
        ax.coastlines()

create_projection()
plt.show()