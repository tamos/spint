''' class for calculating population in specified Radiation
'''

from libpysal.cg import get_points_dist, PointLocator, Point
from libpysal import io
from functools import reduce

class RadiusCalculator:

    def __init__(self, xs, ys, vals):
        points = [Point((x,y)) for x,y in zip(xs,ys)]
        self._points = PointLocator(points)
        self._lookup = {}
        for x,y,val in zip(xs, ys, vals):
            if x not in self._lookup:
                self._lookup[x] = {}
            if y not in self._lookup[x]:
                self._lookup[x][y] = val

    def _get_val(self, coord_tup):
        x,y = coord_tup
        return self._lookup[x][y]

    def _get_points_in_radius(self, point, radius):
        return self._points.proximity(Point(point), radius)

    def _get_vals_in_radius(self, point_coords, radius):
        ''' Fetch the population within `radius` of the given point.
        '''
        return map(lambda x: self._get_val(x),
                    self._get_points_in_radius(point_coords, radius))

    def total_vals_in_radius(self, point_coords, radius):
        vals = self._get_vals_in_radius(point_coords, radius)
        return sum(vals)

    def total_vals_in_radius_between(self, point_1, point_2):
        p1 = Point(point_1)
        dist = get_points_dist(p1, Point(point_2))
        return self.total_vals_in_radius(p1, dist)
