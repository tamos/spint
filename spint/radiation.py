''' Radiation model of mobility

References:

    - Simini et al. 2013

    - Kang et al. 2015 https://doi.org/10.1371/journal.pone.0143500


'''

from .radius_calc import RadiusCalculator
from itertools import combinations
import numpy as np

class BaseRadiation:

    def __init__(self, locations, pops_x, pops_y, pops_val, Nc = None):
        self._radius_size = None
        self._N = sum(pops_val)
        self._Nc = Nc # total num movers
        print('... Creating RadiusCalculator Object...')
        self._pop_calc = RadiusCalculator(pops_x, pops_y, pops_val)
        self.results = None
        self.locations = locations

    def _tot_pop_radius_btwn(self, i, j):
        return self._pop_calc.total_vals_in_radius_between(i,j)

    def _pop_at_loc(self, loc):
        return self._pop_calc._get_val(loc)

    def num_commuters_starting_at(self, i):
        m_i = self._pop_at_loc(i)
        return m_i * (self._Nc / self._N)

    def _model(self):
        return None # each child overloads

    def _run_model(self):
        return None # each child overloads

    def _make_model_func(self):
        #return np.vectorize(self._model)
        return self._model
        # ref https://stackoverflow.com/questions/7701429/efficient-evaluation-of-a-function-at-every-cell-of-a-numpy-array
        # ref https://stackoverflow.com/questions/4231190/python-numpy-tuples-as-elements-of-an-array

    def calculate_all(self):
        self.results = self._run_model()

    @property
    def locations_xy_list(self):
        return list(map(lambda x: x._Point__loc, self.locations))

class Radiation(BaseRadiation):
    ''' this follows the basic form:
    mean(T_ij) = T_i * [ (m_i n_j) / ( (m_i + s_ij) (m_i + n_j + s_ij) ) ]
    as implemented in simini et al 2012
    '''

    def _model(self, pair, nc = None):
        i, j = pair
        m_i = self._pop_at_loc(i)
        n_j = self._pop_at_loc(j)
        s_ij = max(0, self._tot_pop_radius_btwn(i,j) - m_i - n_j)
        if nc:
            T_i = nc
        elif not self._Nc:
                raise ValueError('Must assign or supply value for Nc')
        else:
            T_i = self.num_commuters_starting_at(i)

        result = T_i * ( (m_i * n_j) / ( (m_i + s_ij) * (m_i + n_j + s_ij) ) )
        return (i,j,result)

    def _run_model(self):
        model_func = self._make_model_func()
        pairs = list(combinations(self.locations_xy_list, 2))
        return map(model_func, pairs)

class BaseGeneralized(BaseRadiation): # not yet implemented
    def __init__(self):
        pass

class ProductionIntervention(BaseGeneralized): # not yet implemented
    def __init__(self):
        pass

class AttractionIntervention(BaseGeneralized): # not yet implemented
    def __init__(self):
        pass

class AttractionCompetition(BaseGeneralized): # not yet implemented
    def __init__(self):
        pass
