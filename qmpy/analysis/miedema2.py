from qmpy.analysis.miedema import Miedema, params
from itertools import combinations

class Miedema_Ter:
    def __init__(self, composition):
        assert(isinstance(composition, dict))
        self.composition = composition

    def __call__(self):
        dh_ter = 0
        for i, j in combinations(self.composition.keys(), 2):
            x_i = self.composition[i]
            x_j = self.composition[j]
            h_i_in_j = Miedema.get_ij({i : x_i, j : x_j}) #binary
            h_j_in_i = Miedema.get_ij({j : x_j, i : x_j}) #binary
            dh_ij = x_i * x_j * ((x_j * h_i_in_j) + (x_i * h_j_in_i ))
            dh_ter += dh_ij

        A, B, C = self.composition.keys()
        htrans = [params[A][7], params[B][7], params[C][7]]
        D_htrans = sum(htrans)
        dh_ter += D_htrans
        return round(dh_ter * 0.01036427, 2)