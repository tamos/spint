''' Radiation model of mobility

References:

    - Simini et al. 2013

    - Kang et al. 2015 https://doi.org/10.1371/journal.pone.0143500


'''

from .radius_calc import RadiusCalculator
from itertools import combinations
import numpy as np

class BaseRadiation:

    def __init__(self, locations, pops_x, pops_y, pops_val, cheap = True, Nc = 0.5):
        self._radius_size = None
        self._N = len(pops_x)
        self._Nc = Nc # total number movers
        self._pop_calc = RadiusCalculator(pops_x, pops_y, pops_val)
        self.results = None
        self.locations = locations
        self.cheap = cheap

    def _tot_pop_radius_btwn(self, i, j):
        return self._pop_calc.total_vals_in_radius_between(i,j)

    def _pop_at_loc(self, loc):
        return self._pop_calc._get_val(loc)

    def num_commuters_starting_at(self, i):
        if self.cheap:
            m_i = self._pop_at_loc(i)
            return m_i * (self._Nc / self._N)
        else: # not yet implemented
            return None

    def _model(self):
        return None # each child overloads

    def _run_model(self):
        return None # each child overloads

    def _make_model_func(self):
        #return np.vectorize(self._model)
        return self._model
        # ref https://stackoverflow.com/questions/7701429/efficient-evaluation-of-a-function-at-every-cell-of-a-numpy-array
        # ref https://stackoverflow.com/questions/4231190/python-numpy-tuples-as-elements-of-an-array

    def calculate(self):
        self.results = self._run_model()

    @property
    def locations_xy_list(self):
        return list(map(lambda x: x._Point__loc, self.locations))

class Radiation(BaseRadiation):
    ''' this follows the basic form:
    mean(T_ij) = T_i * [ (m_i n_j) / ( (m_i + s_ij) (m_i + n_j + s_ij) ) ]
    as implemented in simini et al 2012
    '''

    def _model(self, pair):
        i, j = pair
        m_i = self._pop_at_loc(i)
        n_j = self._pop_at_loc(j)
        s_ij = self._tot_pop_radius_btwn(i,j) - m_i - n_j
        T_i = self.num_commuters_starting_at(i)
        result = T_i * ( (m_i * n_j) / ( (m_i + s_ij) * (m_i + n_j + s_ij) ) )
        return (i,j,result)

    def _run_model(self):
        model_func = self._make_model_func()
        pairs = list(combinations(self.locations_xy_list, 2))
        return [model_func(x) for x in pairs]
        #return model_func(pairs)

class Production(BaseRadiation): # not yet implemented
    def __init__(self):
        pass

class Attraction(BaseRadiation): # not yet implemented
    def __init__(self):
        pass
