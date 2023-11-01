# cython: language_level=3
import numpy as np
cimport numpy as np
import cython
 
# type declarations
DTYPE = np.int_
ctypedef int DTYPE_t

@cython.boundscheck(False)
@cython.wraparound(False)
cdef class Packing:
    cdef np.ndarray neighborhood
    cdef np.ndarray vc_solution
    cdef dict satellites
    cdef set reduction 
    cdef int is_violated
    cdef list adj

    def __cinit__(self, int V, list adj, np.ndarray vc_solution):
        self.neighborhood = np.zeros(V, dtype=DTYPE)
        self.vc_solution = vc_solution
        self.reduction = set()
        self.is_violated = 0
        self.adj = list(map(lambda a: list(a), adj))
        self.satellites = dict()

    cpdef np.ndarray get_vc_solution(self):
        return self.vc_solution

    cpdef int packing_is_violated(self):
        return self.is_violated

    cpdef list get_reduction_vertices(self):
        cdef int v, w
        cdef list not_in_cover = list()
        for v in self.reduction:
            w = self.find_reduction_vertice(v)
            if w != -1: not_in_cover.append(w)
        return not_in_cover

    cpdef int has_reduction(self):
        return len(self.reduction)>0

    cpdef void reset_reduction(self):
        self.reduction=set()

    cdef void add_to_reduction(self, int u):
        self.reduction.add(u)
    
    cdef void remove_from_reduction(self, int u):
        self.reduction.discard(u)

    cdef int find_reduction_vertice(self, int u):
        cdef int w = -1
        cdef np.ndarray not_in_cover
        cdef list candy = list()
        if u in self.satellites:
            for s in self.satellites[u]:
                if self.vc_solution[s]==0: 
                    if w == -1:
                        w = s
                    else: candy.append(s)
            if len(candy)>0: print(candy)
        else:
            not_in_cover = np.where(self.vc_solution[self.adj[u]]==False)[0]
            if not_in_cover.size!=1:
                w = 0
                print("# Error in packing reduction", not_in_cover)
            else: 
                w = self.adj[u][not_in_cover[0]]
        return w

    cpdef void init_satellites(self, int v, list neighbors, list adj_list):
        self.update_constraint(neighbors)
        cdef set n_plus
        for w in neighbors:
            n_plus = adj_list[w].difference(adj_list[v])
            n_plus.discard(v)  
            self.neighborhood[w] = len(n_plus)+1
            self.satellites[w]=n_plus

    cpdef void revert_init_satellites(self,neighbors):
        self.neighborhood[neighbors]=0
        cdef int v
        for v in neighbors:
            del self.satellites[v]
        self.revert_update_constraint(neighbors)

    cpdef void update_constraint_v(self, int v, set neighbors):
        self.neighborhood[v] = len(neighbors)+1 # |N[v]| init 
        self.satellites[v]=neighbors
        for u in self.adj[v]:
            if self.neighborhood[u]>1:
                self.neighborhood[u] -=1
                if self.neighborhood[u] == 1: self.is_violated = 1
                elif self.neighborhood[u] == 2: self.add_to_reduction(u)

    cpdef void revert_update_constraint_v(self, int v, set neighbors):
        cdef int u
        self.neighborhood[v] = 0 # delete constraint 
        del self.satellites[v]
        for u in self.adj[v]:
            if self.neighborhood[u]>0:
                if self.neighborhood[u] == 2: self.remove_from_reduction(u)
                self.neighborhood[u] +=1
        self.is_violated = 0

    cpdef void update_constraint(self, list cover_vertices):
        cdef int v
        for v in cover_vertices:
            for u in self.adj[v]:
                if self.neighborhood[u]>1:
                    self.neighborhood[u] -=1
                    if self.neighborhood[u] == 1: self.is_violated = 1
                    elif self.neighborhood[u] == 2: self.add_to_reduction(u)

    cpdef void revert_update_constraint(self, list cover_vertices):
        cdef int v
        for v in cover_vertices:
            for u in self.adj[v]:
                if self.neighborhood[u]>0:
                    if self.neighborhood[u] == 2: self.remove_from_reduction(u)
                    self.neighborhood[u] +=1
        self.is_violated = 0