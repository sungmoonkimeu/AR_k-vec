import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import spherical_poly as sp
import var
# polar to cartesian and vice versa
def get_lon_lat(v, radians = True):
    # get the lon and lat of a vertex in x,y,z
    r_xy = np.sqrt(v[0]**2 + v[1]**2)
    lat = np.arctan2(v[2], r_xy)
    lon = np.arctan2(v[1], v[0])
    if not radians:
        lon = rad2deg(lon)
        lat = rad2deg(lat)
    return(lon, lat)
def get_xyz(lon, lat, radians = True):
    # get the cartesian coordinates for the lon, lat entries on a unit sphere
    if not radians:
        lon = deg2rad(lon)
        lat = deg2rad(lat)
    return (np.array([np.cos(lat) * np.cos(lon),
                      np.cos(lat) * np.sin(lon),
                      np.sin(lat)]))

# angular conversions
def deg2rad(theta):
    if isinstance(theta, list):
        return ([t * var.pi_over_180 for t in theta])
    return theta * var.pi_over_180

def rad2deg(theta):
    if isinstance(theta, list):
        return [t / var.pi_over_180 for t in theta]
    return theta / var.pi_over_180

# algebraic functions
def _cross(a, b):
    """
    Faster than `numpy.cross` as it only does the 3D cross product
    """
    c = np.zeros(var.n_dim)
    c[0] = a[1] * b[2] - a[2] * b[1]
    c[1] = a[2] * b[0] - a[0] * b[2]
    c[2] = a[0] * b[1] - a[1] * b[0]
    return c

def IsCCW(u, v, w):
    '''
    Checks if the smaller interior angle for the great circles of the
    planes u-v and v-w is CCW inputs u, v, w are vectors from the
    origin to the surface of a unit sphere details of approach taken
    from https://mathoverflow.net/questions/193569/
    determining-orientation-of-spherical-polygons:
        Another, even simpler solution of the problem is the
        following: under the assumption that the sphere is centered at
        the origin, all points on the sphere have equal distance from
        the origin and can also be interpreted as vectors of equal
        norm.
        Lets further assume we are given three linearly independent
        vectors u,v,w \in R3 ^ ||u||^2 = ||v||^2 = ||w||^2 of equal
        length, then the sequence (u,v,w) is left turn on the sphere
        (and thus indicates counter clockwise traversal) iff ((v-u) x
        (w-v)) \cdot v > 0
    '''
    return (np.dot(_cross(v - u, w - v), v) > 0)

def v_spoint_to_np(v):
    # check to see if object is s_point if so return just its vertex position
    # otherwise return the original object
    if (isinstance(v, sp.s_point)):
        return (v.v)
    else:
        return (v)

def great_circle_distance(v1, v2, r=None, bearing=False, \
                          geobearing=False, radians=True):
    # calculate the distance between v1 and v2
    v1 = v_spoint_to_np(v1)
    v2 = v_spoint_to_np(v2)
    r1 = np.linalg.norm(v1)
    if r is None: r = r1
    if bearing is False:
        r2 = np.linalg.norm(v2)
        if (abs(r1 - r2) / r2 > 0.01):
            raise ('Vertices in great_circle_distance are on spheres of ' + \
                   'different radii.')
        angle = np.arccos(np.dot(v1, v2) / (r ** 2))
        return (angle * r)
    else:
        if geobearing:  # convert to polar
            v2 = geobearing_transform(radians=radians)
        return (v2 * r)

def geobearing_transform(alpha, radians=True):
    # Alternates between geographic bearings (East = 90, North = 0)
    # and a polar, righthanded system (East = 0, North = 90) where x
    # is East and y is North Note that alpha_geo = 90 - alpha_polar
    # and alpha_polar = 90 - alpha_geo so regardless of the transform
    # the same function call can be used
    if radians:
        east = np.pi / 2
    else:
        east = 90
    return (east - alpha)

def segment_intersection(Sa, Sb, n=None):
    '''Find the intersection of two segments
    each segment Sa and Sb

    Sa and Sb are 2-by-3 vectors, Sa[0,:] and Sa[1,:] giving the
    vertices of the line
    '''
    if n is None:
        n = np.zeros((2, var.n_dim))
        n[0, :] = _cross(Sa[0, :], Sa[1, :])
        n[1, :] = _cross(Sb[0, :], Sb[1, :])
    t = _cross(n[0, :], n[1, :])
    t = t / np.linalg.norm(t)
    sa1 = np.dot(_cross(Sa[0, :], n[0, :]), t)
    sa2 = np.dot(_cross(Sa[1, :], n[0, :]), t)
    sb1 = np.dot(_cross(Sb[0, :], n[1, :]), t)
    sb2 = np.dot(_cross(Sb[1, :], n[1, :]), t)
    intersect = np.sign(sa1) == np.sign(-sa2) and \
                np.sign(sa1) == np.sign(sb1) and \
                np.sign(sa1) == -np.sign(sb2)
    if intersect:
        t = np.sign(sa2) * t
        return (t)
    return False

def get_midpoint(v1, v2, return_spoint=True):
    # Returns an s_point midway along the greatcircle segment
    # connecting v1 and v2
    v1 = v_spoint_to_np(v1)
    v2 = v_spoint_to_np(v2)
    # GC distance from v1 to v2
    d = 0.5 * great_circle_distance(v1, v2)
    alpha = d / np.linalg.norm(v1)
    # find the vector perpendicular to v1
    v1_perp = bearing_great_circle(v1, v2)
    midpoint = v1 * np.cos(alpha) + v1_perp * np.sin(alpha)
    if return_spoint:
        midpoint = sp.s_point(midpoint)
    return (midpoint)

def bearing_great_circle(v1, v2, bearing=False, \
                         geobearing=False, radians=True):
    # finds the bearing from v1 to v2 along a great circle
    #
    # v2 may be a s_point or a numpy array
    #
    # if bearing = False, the function returns the vector normal to v1
    # that, with v1, defines the components for a great_circle
    #
    # If bearing is true, it returns the bearing from v1 that will
    # reach v2 on a great circle path.
    #
    # Depending on if geobearing is True or False, the bearing can be
    # in a coordinate system where east is 0 and north is 90 or pi/2,
    # or in a geographic bearing where north is 0 and east is 90 or
    # pi/2.
    #
    # If radians is True, beta will be on [0, 2 * pi), otherwise it is
    # on [0, 360).
    v1 = v_spoint_to_np(v1)
    v2 = v_spoint_to_np(v2)
    x_vec = _cross(_cross(v1, v2), v1)
    x_vec = x_vec / np.linalg.norm(x_vec)
    if bearing is False:
        return (x_vec * np.linalg.norm(v1))
    else:
        (x_local, y_local, z_local) = local_plane(v1)
        x1 = np.dot(x_vec, x_local)
        y1 = np.dot(x_vec, y_local)
        beta = np.arctan2(y1, x1)
        if geobearing:
            beta = geobearing_transform(beta, radians=True)
        if not radians:
            beta = rad2deg(bearing)
        return (beta)

def interior_angle_points(v1, v2, v3):
    # finds the interior angle of the vertices v1---v2---v3
    # (the smaller of the two)
    n_12 = great_circle_normal(v1, v2)
    n_23 = great_circle_normal(v2, v3)
    return (interior_angle_planes(n_12, n_23))

def interior_angle_planes(n1, n2):
    # finds the internal angle of intersection (the smaller of two)
    # between planes n1, and n2
    d = np.dot(n1, n2)
    # Numerically, sometimes abs(d) > 1 in which case just treat it as
    # 1
    if np.abs(d) > 1:
        alpha = np.arccos(-np.sign(d))
    else:
        alpha = np.arccos(-d)
    return (alpha)

def great_circle_normal(v1, v2, normalize=True):
    # returns the normal to the plane of the great circle connecting
    # v1--v2
    v1 = v_spoint_to_np(v1)
    v2 = v_spoint_to_np(v2)
    n = _cross(v1, v2)
    if normalize:
        n = n / np.linalg.norm(n)
    return (n)

def order_vertices(vertices, start=0):
    # orders a set of vertices
    vertice_order = [start, ]
    v_centre = np.mean(vertices, axis=0)
    v_centre = v_centre / np.linalg.norm(v_centre)
    A = []
    for j in range(1, len(vertices)):
        i = (j + start) % len(vertices)
        # order vertices in counter clockwise order, starting with "start"
        A.append(interior_angle_points(vertices[i], v_centre, vertices[start]))
        # if not CCW then got the wrong angle from cosine and needs correcting
        if not IsCCW(vertices[i], v_centre, vertices[start]):
            A[-1] = 2 * np.pi - A[-1]

    vertice_order = np.argsort(A)
    vertices_sorted = np.empty_like(vertices)
    vertices_sorted[0, :] = vertices[start, :]
    for i in range(len(vertice_order)):
        vertices_sorted[i + 1, :] = vertices[start + 1 + vertice_order[i], :]
    return vertices_sorted

def local_plane(point):
    '''returns a local "flat earth" coordinate system where x is East, y
    in North, and z is Up
    '''
    # temp gives the normal to the plane including the North Pole and
    # the location
    v = v_spoint_to_np(point)
    North = sp.s_point([0, 0, 1])
    local_z = v / np.linalg.norm(v)
    temp = _cross(local_z, North.v)
    # check that the point is not the North Pole
    if np.linalg.norm(temp) == 0:
        local_x = np.array([1, 0, 0])
        local_y = np.array([0, 1, 0])
    else:
        temp /= np.linalg.norm(temp)
        # crossing temp with the radar location gives the unit vector
        # pointing north from the radar location
        local_y = _cross(temp, local_z)
        # then we can calculate x and z
        local_x = _cross(local_y, local_z)
    return (local_x, local_y, local_z)

# plotting
def create_projection(lon=-90, lat=90, grid=True, \
                      proj=ccrs.Orthographic, figsize=None, \
                      set_global=True):
    if figsize is not None:
        fig = plt.figure(figsize=figsize)
    else:
        fig = plt.figure()
    projection_lon = lon
    projection_lat = lat
    projection = proj(projection_lon, projection_lat)
    transform = ccrs.Geodetic()
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    if set_global:
        ax.set_global()
    if grid:
        ax.gridlines()
    return (fig, ax, transform)

def add_geography(ax, coast=True, land=True, ocean=False, alpha=1):
    if ocean:
        ax.add_feature(cfeature.OCEAN, zorder=0, alpha=alpha)
    if land:
        ax.add_feature(cfeature.LAND, zorder=0, edgecolor=[0.6, 0.6, 0.6], alpha=alpha)
    if coast:
        ax.coastlines(alpha=alpha)

def get_arc_values(up, vp, N):
    # before use, need to figure out how to implement this best for
    # degrees/xyz and filling versus plotting...
    u = up / np.linalg.norm(up)
    v = vp / np.linalg.norm(vp)
    angle = np.arccos(np.dot(v, u))
    nhat = great_circle_normal(u, v)
    w = _cross(_cross(u, v), u)
    w = w / np.linalg.norm(w)
    points = np.zeros((N, 2))
    for n in range(N):
        alpha = n * angle / (N - 1)
        p = u * np.cos(alpha) + w * np.sin(alpha)
        points[n, :] = get_lon_lat(p, radians=False)
    return (points)

def plot_arc(up, vp, N, *args, **kwargs):
    # plot N points along a great circle segment starting at sp1 and
    # ending at sp2
    u = up / np.linalg.norm(up)
    v = vp / np.linalg.norm(vp)
    angle = np.arccos(np.dot(v, u))
    nhat = great_circle_normal(u, v)
    w = _cross(_cross(u, v), u)
    w = w / np.linalg.norm(w)
    lon = np.zeros(N)
    lat = np.zeros(N)
    for n in range(N):
        alpha = n * angle / (N - 1)
        p = u * np.cos(alpha) + w * np.sin(alpha)
        (lon[n], lat[n]) = get_lon_lat(p, radians=False)
    plt.plot(lon, lat, *args, **kwargs)
