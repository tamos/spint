''' class for calculating population in specified Radiation
'''

from libpysal.cg import get_points_dist, PointLocator, Point
from libpysal import io
from functools import reduce
import numpy as np

class RadiusCalculator:

    def __init__(self, xs, ys, vals):
        self._points = PointLocator([Point((x,y)) for x,y in zip(xs,ys)])
        print("... Point Locator Created ...")
        self._lookup = {}
        npoints = len(vals)
        progress_interval = npoints / 20
        count_of = 0
        progress_benchmark = 0
        print("... Building Lookup Tree ...")
        for x,y,val in zip(xs, ys, vals):
            if x not in self._lookup:
                self._lookup[x] = {}
            if y not in self._lookup[x]:
                self._lookup[x][y] = val
            count_of += 1
            if count_of > progress_benchmark:
                pct_progress = round((count_of / npoints) * 100, 3)
                print("...{}% of points added...".format(pct_progress))
                progress_benchmark += progress_interval


    def _get_val(self, coord_tup):
        nearest = self._points.nearest(coord_tup)
        n_x, n_y = nearest._Point__loc
        return self._lookup[n_x][n_y]

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
