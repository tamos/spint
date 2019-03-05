''' Radiation model of mobility

References:

    - Simini et al. 2013

    - Kang et al. 2015 https://doi.org/10.1371/journal.pone.0143500


'''

from .radius_calc import RadiusCalculator
from itertools import combinations
import numpy as np
from scipy.optimize import minimize
from functools import partial

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
        if nc:
            T_i = nc
        elif not self._Nc:
                raise ValueError('Must assign or supply value for Nc')
        else:
            T_i = self.num_commuters_starting_at(i)

        return self._model_base(i,j,T_i)

    def _model_base(self, i, j, T_i):
        m_i = self._pop_at_loc(i)
        n_j = self._pop_at_loc(j)
        s_ij = max(0, self._tot_pop_radius_btwn(i,j) - m_i - n_j)
        result = T_i * ( (m_i * n_j) / ( (m_i + s_ij) * (m_i + n_j + s_ij) ) )
        return (i,j,result)

    def _run_model(self):
        model_func = self._make_model_func()
        pairs = list(combinations(self.locations_xy_list, 2))
        return map(model_func, pairs)


class CommuterPredictingRadiation(Radiation): # solves for N_c iteratively

    def _fit_model(self, plot = False, **kwargs):
        bnds = [(0,None)] * len(self.locations) # must never be negative
        mean_pop = self._N / len(self.locations)
        init_guess = [mean_pop] * len(self.locations)
        rv = minimize(self._objective_func, x0 = init_guess, bounds = bnds,
                        **kwargs)
        if plot:
            x  = [i._Point__loc[0] for i in self.locations]
            y  = [i._Point__loc[1] for i in self.locations]
            z  = rv.x
            return (x,y,z), rv # gives you all you need to plot them
        else:
            return rv

    def _objective_func(self, commuter_guesses):

        locations_dict = dict(zip(self.locations_xy_list, commuter_guesses))

        candidate_T_is = {}
        nlocs = len(commuter_guesses)
        result_dict = {}

        for each_pair in list(combinations(self.locations_xy_list, 2)):
            self._Nc = locations_dict[each_pair[0]]
            i, j, result = self._model(each_pair)
            if i not in result_dict:
                result_dict[i] = {}
            if j not in result_dict[i]:
                result_dict[i][j] = result
            T_i  = self.num_commuters_starting_at(i)
            if i not in candidate_T_is:
                candidate_T_is[i] = T_i

        total_error = 0
        n_results = 0

        for each_i in candidate_T_is.keys():
            total_i = sum( [k for k in result_dict[each_i].values()] )
            total_error += (total_i - candidate_T_is[each_i])**2
            n_results += 1
        print("Error is", total_error/n_results)
        return total_error/n_results


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
