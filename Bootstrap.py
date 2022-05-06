from mpmath import mp, det
from sympy.matrices import zeros


class Bootstrap:

    def __init__(self, E, Recurrence_Func, Initial_conition, Matrix_el_Fnc, Default_N = 30):
        
        self.E = E
        self.N = Default_N
        self.rec = Recurrence_Func
        self.computed = Initial_conition
        self.matrix_el_fnc = Matrix_el_Fnc

        # Build Hankel matrix
        self.H = mp.matrix(zeros(self.N + 1, self.N + 1))
        for i in range(self.N + 1):

            for j in range(i + 1):

                try:
                    self.H[i, j] = self.H[j, i] = self.matrix_el_fnc(i, j)
                
                except KeyError:

                    if i + j not in self.computed:
                        self.computed[i + j] = self.rec(i + j)
                    
                    if i + 1 not in self.computed:
                        self.computed[i + 1] = self.rec(i + 1)
                    
                    if j + 1 not in self.computed:
                        self.computed[j + 1] = self.rec(j + 1)
                    
                    self.H[i, j] = self.H[j, i] = self.matrix_el_fnc(i, j)
    
    
    def is_Eigen(self, n: int, allowed = 1, not_allowed = 0):
        
        i = 1

        while det((self.H)[0: i, 0: i]) >=0 and i <= n:
            
            i += 1
        
        if i - 1 == n:
            return allowed
        
        else:
            return not_allowed
    

    def is_Eigen_Multi_K(self, k_list, allowed = 1, not_allowed = 0):

        i = 1
        k_allowed = []
        k_list = sorted(k_list)

        while i <= k_list[-1]:
            
            tmp_minor = det(self.H[0: i, 0: i])
            
            if tmp_minor >= 0 and i in k_list:

                k_allowed.append(allowed)
                i += 1
                
            elif tmp_minor < 0:
                break
            
            else:
                i += 1
        
        return k_allowed + [not_allowed]*(len(k_list) - len(k_allowed))