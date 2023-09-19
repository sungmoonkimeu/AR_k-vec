import numpy as np
import warnings
import copy
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.spatial import SphericalVoronoi as SV
import json
import util
import var
#import spherical-geometry as sphere
class s_point():
    # Point on a sphere, defined as longitude and latitude
    def __init__(self, coordinates, lonlatorder = True, r = var.R_earth):
        if type(coordinates) is dict:
            # assume JSON format
            self.lon = coordinates["longitude"]
            self.lat = coordinates["latitude"]
            self.r = coordinates["r"]
            self.phi = coordinates["phi"]
            self.theta = coordinates["theta"]
            self.v = np.array(coordinates["v"])
        else:
            # print("print")
            if len(coordinates) == 2:
                if lonlatorder:
                    # assume lon lat
                    self.lon = coordinates[0]
                    self.lat = coordinates[1]
                else:
                    self.lat = coordinates[0]
                    self.lon = coordinates[1]
                # convert to theta, phi and define location in
                # cartesian coordinates
                self.phi = util.deg2rad(self.lon)
                self.theta = util.deg2rad(self.lat)
                # get unit vector of vertex v = (x, y, z)
                self.v = util.get_xyz(self.phi, self.theta) * r
                self.r = r
            elif len(coordinates) == 3:
                # assume x,y,z
                self.v = np.array(coordinates)
                self.r = np.linalg.norm(self.v)
                (self.phi, self.theta) = util.get_lon_lat(self.v / r, \
                                                          radians = True)
                (self.lon, self.lat) = util.get_lon_lat(self.v / r, \
                                                        radians = False)
    def to_dict(self):
        # creates a dictionary with lon, lat and r to allow saving in JSON
        return({"type": "s_point",
                "longitude": self.lon,
                "latitude": self.lat,
                "r": self.r,
                "phi": self.phi,
                "theta": self.theta,
                "v": self.v.tolist()})
    def save(self, filename):
        with open(filename, 'w') as fid:
            json.dump({"1": self.to_dict()}, fid)
    def __eq__(self, other):
        # Check if another object is the same
        if not isinstance(other, s_point):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return(self.lon == other.lon and \
               self.lat == other.lat and \
               self.r == other.r and \
               self.phi == other.phi and \
               self.theta == other.theta and \
               np.all(self.v == other.v))

class s_poly():
    # Class to manage manipulation and calculations on convex
    # counterclockwise spherical polygons on a sphere
    def __init__(self, points):
        if type(points) is dict:
            self.n_edge = points["n_edge"]
            self.v = np.array(points["v"])
            self.lon = np.array(points["lon"])
            self.lat = np.array(points["lat"])
            self.r = points["r"]
            self.n = np.array(points["n"])
            self.alpha = np.array(points["alpha"])
            self.centre = self.center = s_point(points["centre"])
        else: # points is a list of s_point objects
            self.n_edge = len(points)
            # generate v, a matrix of vertices
            # v[i,:] gives the (x,y,z) coordinates of vertex i
            self.v = np.zeros((len(points), var.n_dim))
            self.lon = np.zeros(len(points))
            self.lat = np.zeros(len(points))
            self.r = np.linalg.norm(points[0].v)
            for p in range(self.n_edge):
                self.v[p,:] = points[p].v
                self.lon[p] = points[p].lon
                self.lat[p] = points[p].lat
            # ensure vertices are listed in CCW order and that
            # spherical poly is convex if it is not convex, the
            # function returns an error ending the program if it is
            # convex but CW, the function returns "False" if it convex
            # and CCW, the function returns "True"
            if not self.IsConvexCCW():
                warnings.warn('Polygon is clockwise, so reversing order' + \
                              'of points for calculations.')
                self.v = np.flipud(self.v)
            self.centre = self.center = s_point(self.moment(normed = True), \
                                                            r = self.r)
            # normals to planes containing great circles
            self.n = np.zeros((len(points), var.n_dim))
            for i in range(self.n_edge):
                j = (i+1) % self.n_edge
            self.n[i,:] = util.great_circle_normal(self.v[i,:], \
                                                   self.v[j,:])
            # internal angles for each vertex on the great circle
            self.alpha = np.zeros(len(points))
            for j in range(self.n_edge):
                i = (j - 1) % self.n_edge
            # minus sign is to deal with possibility of angle
            # being 180 degrees
            self.alpha[i] = util.interior_angle_planes(self.n[i,:], \
                                                       self.n[j,:])
    def to_dict(self):
        return({"type": "s_poly",
                "n_edge": self.n_edge,
                "v": self.v.tolist(),
                "lon": self.lon.tolist(),
                "lat": self.lat.tolist(),
                "r": self.r,
                "n": self.n.tolist(),
                "alpha": self.alpha.tolist(),
                "centre": self.centre.to_dict()
                })

    def save(self, filename):
        with open(filename, 'w') as fid:
            json.dump({"1": self.to_dict()}, fid)

    def __eq__(self, other):
        # Check if another object is the same
        if not isinstance(other, s_poly):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return(np.all(self.lon == other.lon) and \
               np.all(self.lat == other.lat) and \
               self.r == other.r and \
               np.all(self.v == other.v) and \
               self.n_edge == other.n_edge and \
               np.all(self.n == other.n) and \
               np.all(self.alpha == other.alpha) and \
               self.centre == other.centre)

    def moment(self, normed = True):
        # Calculate the centre of the spherical polygon by calculating
        # the moment of the region inside the spherical polygon, following
        # https://stackoverflow.com/questions/19897187/
        # locating-the-centroid-center-of-mass-of-spherical-polygons
        # which references, Victor Alexandrov (2004) Problem 10957,
        # American Mathematical Monthly, 2004, 111, no. 2, p. 172.
        # Based on an elegant solution using Stokes Theorem...(!)
        #
        # From SO: The moment of the region to the left of a closed
        # path on the unit sphere is half the integral of the leftward
        # unit vector as you walk around the path....In particular,
        # for a spherical polygon, the moment is the half the sum of
        # (a x b) / ||a x b|| * (angle between a and b) for each pair
        # of consecutive vertices a,b.
        #
        # Thinking about the line integral in Stokes Theorem, a x b / ||a
        # x b|| is the leftward pointing unit vector and the angle
        # between a and b is the arc length of the side.
        #
        # not normalized gives the moment
        # normalized gives the unit vector in the direction of the moment
        moment = np.zeros(3)
        for i in range(self.n_edge):
            j = (i + 1) % self.n_edge
            N = util._cross(self.v[i,:], self.v[j,:])
            N = N / np.linalg.norm(N)
            alpha = np.arccos(np.dot(self.v[i,:], self.v[j,:]) /
                              (np.linalg.norm(self.v[i,:]) * \
                               np.linalg.norm(self.v[j,:])))
            moment += N * alpha / 2
        if normed:
            moment /= np.linalg.norm(moment) / self.r

        return(moment)

    def IsConvexCCW(self):
        '''
        check that all sides of the polygon have the same orientation and
        that the s_poly is convex
        '''
        alpha_ccw = []
        for j in range(self.n_edge):
            i = (j - 1) % self.n_edge
            k = (j + 1) % self.n_edge
            alpha_ccw.append(util.IsCCW(self.v[i,:], self.v[j,:], self.v[k,:]))
        # check that they are all the same orientation and thus convex
        if alpha_ccw.count(alpha_ccw[0]) == self.n_edge:
            if alpha_ccw[0]:
                return True
            else:
                return False
        else:
            raise('Non-convex series of points sent ' + 'as verties to s_poly class')

    def SphericalExcess(self):
        # returns the spherical excess of the polygon
        return np.sum(self.alpha) - (self.n_edge - 2) * np.pi
    def Area(self, norm = True):
        # returns the surface area of the polygon if normalized (norm
        # = True) it is returned as the fraction of the sphere if not
        # normalized (norm = False) it is returned as the area
        A = self.SphericalExcess()
        if norm:
            A /= 4 * np.pi
        else:
            A *= self.r**2
        return(A)
    def get_lon_lat_point_list(self, N = 2):
        lon = np.array([])
        lat = np.array([])
        for i in range(self.n_edge):
            j = (i + 1) % self.n_edge
            p = util.get_arc_values(self.v[i,:], self.v[j, :], N)
            lon = np.append(lon, p[:,0])
            lat = np.append(lat, p[:,1])
        return(lon, lat)

    def plot(self, N, *args, **kwargs):
        lon, lat = self.get_lon_lat_point_list(N = N)
        plt.plot(lon, lat, *args, **kwargs)
        ## plot N points along each great circle section of polygon
        #for i in range(self.n_edge):
        # j = (i + 1) % self.n_edge
        # util.plot_arc(self.v[i,:], self.v[j,:], N, *args, **kwargs)
    def fill(self, N, *args, **kwargs):
        lon, lat = self.get_lon_lat_point_list(N = N)
        plt.fill(lon, lat, *args, **kwargs)
class small_circle():
    def __init__(self, parameters, N = 30, ax_limits = None):
        # define a small circle of radius (on sphere surface) rho and
        # centred at p.
        # Assume radius of sphere is R = |p|,
        # so arc length is rho / R.
        # If ax_limits is None, the entire small circle is described by a
        # spherical poly composed of N points. If not, the spherical
        # poly object is the small circle segment plus its pole at p
        # are returned.
        if type(parameters) is dict:
            self.N = parameters["N"]
            self.p = np.array(parameters["p"])
            self.rho = parameters["rho"]
            self.R = np.array(parameters["R"])
            self.alpha = parameters["alpha"]
            self.r = parameters["r"]
            self.d = parameters["d"]
            self.c = np.array(parameters["c"])
            self.n = np.array(parameters["n"])
            self.u = np.array(parameters["u"])
            self.v = np.array(parameters["v"])
            self.plane_const = parameters["plane_const"]
            self.poly = s_poly(parameters["poly"])
        else:
            p = parameters[0]
            rho = parameters[1]
            self.N = N
            self.p = p
            self.rho = rho
            self.R = np.linalg.norm(p)
            self.alpha = self.rho / self.R # arc length
            # shortest distance from axis defined by p to sphere edge in plane
            self.r = self.R * np.sin(self.alpha)
            # containing small circle
            self.d = np.sqrt(self.R**2 - self.r**2)
            # vector to "centre" of plane of small circle
            self.c = self.p/self.R * self.d
            # Equation for a plane with normal n, p = x \cdot n =
            # const. We know that c is parallel to the normal of the
            # plane so c \cdot n = |c|. Thus, the equation of the
            # plane is x \cod n = |c|
            self.plane_const = np.linalg.norm(self.c)
            self.n = self.c / self.plane_const
            # Find a vector perpendicular to n.
            # To start, take x,y,z unit vectors and use whichever has
            # the smallest as the dot product

            # First cross product
            u = np.array([1, 0, 0])
            v = np.array([0, 1, 0])
            w = np.array([0, 0, 1])
            if abs(abs(np.dot(self.n, w)) - 1) < 10 ** -4:
                # small circle is at one of the poles set the 0th and
                # 90th meridians as the two orthogonal axes
                u = util._cross(v, self.n)
                self.u = u * self.r / np.linalg.norm(u)
                v = util._cross(self.n, self.u)
                self.v = v * self.r / np.linalg.norm(v)
            else:
                # set east and north on the sphere surface as the two
                # orthogonal axes
                v = util._cross(self.n, w)
                v = v / np.linalg.norm(v)
                self.v = util._cross(v, self.n) * self.r
                self.u = util._cross(self.v, self.n)
            # create s_poly with N points
            vertices = self.points(self.N, ax_limits=ax_limits)
            if ax_limits is not None:
                if (ax_limits[1] != np.pi):
                    vertices.append(s_point(self.p))
            # else:
            # vertices = self.points(N+1, ax_limits = ax_limits)
            # vertices.pop(-1)
            self.poly = s_poly(vertices)

    def to_dict(self):
        return ({"type": "small_circle",
                 "N": self.N,
                 "p": self.p.tolist(),
                 "rho": self.rho,
                 "R": self.R,
                 "alpha": self.alpha,
                 "r": self.r,
                 "d": self.d,
                 "c": self.c.tolist(),
                 "n": self.n.tolist(),
                 "plane_const": self.plane_const,
                 "u": self.u.tolist(),
                 "v": self.v.tolist(),
                 "poly": self.poly.to_dict(),
                 })

    def save(self, filename):
        with open(filename, 'w') as fid:
            json.dump({"1": self.to_dict()}, fid)

    def points(self, N, ax_limits=None):
        if ax_limits is None:
            Omega = np.linspace(0, 2 * np.pi, N, endpoint=False)

        else:
            Omega = np.linspace(ax_limits[0], ax_limits[0] + ax_limits[1], N)
        x = []
        for omega in Omega:
            x.append(s_point(self.c + self.u * np.cos(omega) + \
                             self.v * np.sin(omega)))
        return (x)

    def __eq__(self, other):
        # Check if another object is the same
        if not isinstance(other, small_circle):
            # don't attempt to compare against unrelated types
            return NotImplemented
        return(self.N == other.N and \
               np.all(self.p == other.p) and \
               self.rho == other.rho and \
               self.R == other.R and \
               self.alpha == other.alpha and \
               self.r == other.r and \
               self.d == other.d and \
               np.all(self.c == other.c) and \
               np.all(self.n == other.n) and \
               self.plane_const == other.plane_const and \
               self.poly == other.poly and \
               np.all(self.u == other.u) and \
               np.all(self.v == other.v))

def load(filename, type = dict):
    with open(filename, 'r') as fid:
        input_data = json.load(fid)
    values = {}
    for key in input_data:
        if input_data[key]["type"] == "s_point":
            values[key] = s_point(input_data[key])
        elif input_data[key]["type"] == "s_poly":
            values[key] = s_poly(input_data[key])
        elif input_data[key]["type"] == "small_circle":
            values[key] = small_circle(input_data[key])
    if len(input_data.keys()) == 1:
        return(values[key])
    if type == dict:
        return(values)
    if type == list:
        return[values[key] for key in values.keys()]
def save(values, filename):
    if isinstance(values, dict):
        values = {i: values[i].to_dict() for i in values.keys()}
    if isinstance(values, list):
        values = {i: v.to_dict() for (i,v) in enumerate(values)}
    if isinstance(values, dict):
        with open(filename, 'w') as fid:
            json.dump(values, fid)
    elif isinstance(values, s_point) or \
        isinstance(values, s_poly) or \
        isinstance(values, small_circle):
        values.save(filename)
    else:
        raise('Unable to save object of type %s' % type(values))

def InPoly(object, poly):
    '''
    Determine if an object is in a spherical polygon
    the object can be an instance s_point or s_poly
    '''
    if(isinstance(object, s_point)):
        return(PointInPoly(object.v, poly))
    elif(isinstance(object, s_poly)):
        vertices = []
        for i in range(object.n_edge):
            if PointInPoly(object.v[i,:], poly):
                j = i % object.n_edge # bookkeeping between vertices and edges
                vertices.append(j)
        if len(vertices) > 0:
            return vertices
        else:
            return False
    else:
        raise TypeError('object type %s not supported' % type(object))

def PointInPoly(point, poly):
    '''
    For a convex ccw spherical polygon, normals to all the planes
    should be pointing into the polygon so check that the dot product
    of these normal vectors with the points are all positive.
    '''
    point = util.v_spoint_to_np(point)
    for i in range(poly.n_edge):
        # "On the line is not in!" -every ultimate player ever
        # ...well, it is here...
        inpoly = (np.dot(point, poly.n[i,:]) >= 0)
        if not inpoly:
            return(inpoly)
    return(inpoly)

def poly_intersection(A,B, return_poly = True):
    '''Find polygon intersection of spherical polygons A and B.
    Inputs:
    A, B - s_poly objects descring a pair of spherical polygons
    return_poly - optional boolian variable that defines the type of
    output (default is True)
    Outputs:
    If the flag return_poly is True, then the function returns
    an s_poly object describing that intersection.

    If the flag return_poly is False, then the function retuns
    a boolean value indicating whether the polygons intersect.
    Algorithm Description:
    First, we identify all potential points that overlap.

    All potential types of intersection shown below are considered

    o---------------o
    | |
    o---------------+------o
    |\ | |
    | \ o |
    | \ / |
    | \ / |
    | \ / |
    | \ / |
    | \ / |
    | \ / |
    o--------o-------------o
    Second, the points are ordered so that all points together define
    a CCW convex spherical polygon. Since the two original polygons
    are CCW convex spherical polygons, so too is their intersection.
    The special case where the one polygon is completely inside
    another returns the encased polygon.
    '''
    AinB = InPoly(A,B)
    BinA = InPoly(B,A)
    vertices = []
    # find vertices of A in B
    if(isinstance(AinB, list)):
        if not return_poly:
            return True
        if len(AinB) == A.n_edge:
            intersection = copy.deepcopy(A)
            return intersection
        else:
            vertices += [A.v[i] for i in AinB]
    # find vertices of B in A
    if(isinstance(BinA, list)):
        if not return_poly:
            return True
        if len(BinA) == B.n_edge:
            intersection = copy.deepcopy(B)
            return intersection
        else:
            vertices += [B.v[i] for i in BinA]
    # find intersections of edges
    Sa = np.zeros((2, var.n_dim))
    Sb = np.zeros((2, var.n_dim))
    n_ab = np.zeros((2, var.n_dim))
    for i in range(A.n_edge):
        Sa[0,:] = A.v[i,:]
        Sa[1,:] = A.v[(i+1)%A.n_edge,:]
        n_ab[0,:] = A.n[i,:]
        for j in range(B.n_edge):
            Sb[0,:] = B.v[j,:]
            Sb[1,:] = B.v[(j+1)%B.n_edge,:]
            n_ab[1,:] = B.n[j,:]
            t = util.segment_intersection(Sa, Sb, n=n_ab)

            if isinstance(t, np.ndarray):
                vertices += [t,]
    if len(vertices) == 0:
        if return_poly:
            return None
        else:
            return False
    else:
        if not return_poly:
            return True

    vertices = np.array(vertices)
    vm = np.mean(vertices, axis=0)
    # for testing only - randomize order of vertices
    #np.random.shuffle(vertices)
    vertices_sorted = util.order_vertices(vertices)
    points = [s_point(v) for v in vertices_sorted]
    k = 0
    try:
        intersection = s_poly(points)
        return intersection
    except:
        return None

def voronoi_cells(s_points, radius = var.R_earth, centre = np.array([0,0,0])):
    # Calculates the voronoi points on a sphere
    # based on docs for the scipy package
    # https://docs.scipy.org/doc/scipy/reference/generated/
    # scipy.spatial.SphericalVoronoi.html
    N = len(s_points)
    # create array of points from the s_points,
    points = np.zeros((N, var.n_dim))
    for i in range(len(s_points)):
        points[i,:] = s_points[i].v
    V = SV(points, radius, centre)
    # for some reason there are a number of points at the same
    # location so rather than creating a bunch of s_points, one must
    # make a bunch of s_points for each region, removing the
    # repetitive points

    # create s_poly list of voronoi cells
    s_polys = []
    for i in range(N):
        v_points = np.unique(V.vertices[V.regions[i]], axis = 0)
        point_list = []
        for j in range(len(v_points)):
            point_list.append(v_points[j])
        point_list = util.order_vertices(np.array(point_list))
        s_point_list = [s_point(p) for p in list(point_list)]
        s_polys.append(s_poly(s_point_list))
    return(s_polys)

def interpolate(from_polys, to_polys, from_values, \
                spacing_limit = None, area_threshold = 0.9995, \
                progress_reporting = True):
    '''
    Interpolate values from one grid of spherical polygons to another,
    from the 'from_polys' list to the 'to_polys' list, using values in
    the 'from_values' list len(from_polys) must equal len(from_data).
    Returns a list of length len(to_polys) with the value for each
    element in the list
    INPUTS:
    from_polys - list of s_poly objects defining the original grid
    to_polys - list of s_poly objects defining the new grid
    from_values - list of values for the s_poly objects in from_polys
    NOTE: len(from_polys) and len(from_values) must be the same.
    OPTIONAL INPUTS:
    spacing_limit - if set, will only look for intersections of
    s_polys with centres within this limit
    area_threshold - will not seek intersections with polygons in the
    "from_poly" list if this fraction of their area
   is aready accounted for in the new grid
    progress_reporting - if True, updates progress of the
    interpolation routine since some
    interpolations can take considerable time
     OUTPUT:
     to_values - a list of values for each s_poly in the list to_polys
     '''
    if len(from_polys) != len(from_values):
        raise ('List of s_poly objects and list of values in ' + \
               'spherical_poly.interpolate must be same length.')

    P = len(from_polys)
    from_areas_frac_interpolated = np.zeros(P)

    N = len(to_polys)
    to_values = np.zeros(N)
    to_areas = [to_polys[n].Area() for n in range(N)]
    for n in range(N):
        to_poly_matrix = to_polys[n].v[:, np.newaxis, :]
        if progress_reporting:  # sometimes interpolation takes a while
        # so this just gives progess
            if n % 5 == 0: print('.', end='')
            if (n + 1) % 120 == 0: print('%i/%i' % (n + 1, N))
    for p in range(P):
        from_poly_matrix = from_polys[p].v[np.newaxis, :, :]
        if (np.linalg.norm(to_poly_matrix - from_poly_matrix, \
                       axis=2) < spacing_limit).any():
            I = poly_intersection(to_polys[n], from_polys[p], \
                              return_poly=True)
            if I is not None:
                Iarea = I.Area()
                to_values[n] += from_values[p] * Iarea / to_areas[n]
                from_areas_frac_interpolated[p] += Iarea / to_areas[n]
    return (to_values)
