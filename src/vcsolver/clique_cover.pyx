# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
import numpy as np

cimport numpy as np

DTYPE = int
ctypedef np.int_t DTYPE_t

cdef class CliqueCover:
    cdef np.ndarray clique
    cdef np.ndarray clique_size
    cdef np.ndarray tmp_size
    cdef list adj_list
    cdef np.ndarray vc_solution
    cdef int c_lb
    cdef list new_in_cover
    cdef set new_update_verices
    cdef list revert_stack

    def __cinit__(self, int V, list adj_list, np.ndarray vc_solution):
        self.clique = np.arange(V, dtype=DTYPE) # everyone is its own clique
        self.clique_size = np.ones(V, dtype=DTYPE) 
        self.tmp_size = np.zeros(V, dtype=DTYPE)
        self.adj_list = adj_list
        self.vc_solution = vc_solution
        self.c_lb = -1
        self.new_in_cover = list()
        self.new_update_verices = set()
        self.revert_stack = list()

    cpdef int get_c_lb(self):
        return self.c_lb

    cpdef void update_c_lb(self, int x, int y):
        cdef int c_x, c_y
        c_x, c_y = self.clique[x], self.clique[y]
        self.clique_size[c_x] -= 1
        self.clique_size[c_y] -= 1
        if self.clique_size[c_x]>0: self.c_lb -= 1
        if self.clique_size[c_y]>0: self.c_lb -= 1
        #print("# merge: ", x, c_x, self.clique_size[c_x], y, c_y, self.clique_size[c_y], v, self.clique[v], self.clique_size[self.clique[v]])
        
        #return c_x, c_y

    cpdef void undo_update_c_lb(self, int x, int y):
        cdef int c_x, c_y
        c_x, c_y = self.clique[x], self.clique[y]
        #print("# merge_revert: ", x, c_x, self.clique_size[c_x], y, c_y, self.clique_size[c_y])
        if self.clique_size[c_x]>0: self.c_lb += 1
        if self.clique_size[c_y]>0: self.c_lb += 1
        self.clique_size[c_x] += 1
        self.clique_size[c_y] += 1

    cpdef void reset_records(self):
        self.new_in_cover, self.new_update_verices = list(), set()

    cpdef void record_vertices(self, list in_cover, list neighborhood):
        self.new_in_cover.extend(in_cover)
        self.new_update_verices.update(neighborhood) # its a set

    cpdef int update(self):
        self.handle_in_cover_update(self.new_in_cover)
        cdef list mapping = self.update_neighbors(list(self.new_update_verices))
        self.revert_stack.append((self.new_in_cover, mapping))
        self.reset_records()
        return self.c_lb

    cpdef void revert_update(self):
        cdef list in_cover, mapping
        in_cover, mapping = self.revert_stack.pop()
        self.undo_update_neighbors(mapping)
        self.handle_in_cover_update_invert(in_cover)

    cdef void handle_in_cover_update(self, list vertices):
        cdef int v, c
        for v in vertices:
            c = self.clique[v]
            self.clique_size[c] -=1
            if self.clique_size[c]==0: continue # -1 +1 removed clique 
            self.c_lb -=1

    cdef void handle_in_cover_update_invert(self, list vertices):
        cdef int v, c
        for v in vertices:
            c = self.clique[v]
            self.clique_size[c] += 1
            if self.clique_size[c]==1: continue #-1 +1 add vertex and clique 
            self.c_lb += 1

    cdef void undo_update_neighbors(self, list mapping):
        cdef int v, old_c, imp
        for v, old_c, imp in mapping:
            if self.clique[v] != -1:
                c = self.clique[v]
                self.clique_size[c] -=1
            self.c_lb -=imp # undo posible improvement 
            self.clique[v]=old_c
            self.clique_size[old_c] +=1

    cdef list update_neighbors(self, list vertices):
        cdef list changes = list()
        cdef int v, old_c, imp
        for v in vertices:
            if self.vc_solution[v]: continue
            elif len(self.adj_list[v])==0: #if vertex is isolated
                #if self.clique[v] == -1: print("Error")
                imp, old_c = 0, self.clique[v]
                self.clique_size[old_c] -= 1
                self.clique[v] = -1
            else:   # if vertex has update neighborhood check again for impovement 
                imp, old_c = self.clique_inspect(v)
            
            if old_c != -1:
                self.c_lb += imp
                changes.append((v, old_c, imp))
        return changes

    def clique_lb(self, int init, list deg_bags):
        if init==1:
            self.clique[:], self.clique_size[:], self.c_lb = np.arange(len(self.clique)), 1, 0
        #if c_lb == -1: # during banching TODO 
        #    for v in self.inspect_for_clique_lb:
        #        c_lb += self.clique_inspect(v)
        #    return c_lb
        cdef set bag 
        cdef int v, imp, old_c
        for bag in deg_bags[1:]: # init before branching 
            for v in bag:
                imp, old_c = self.clique_inspect(v)
                self.c_lb += imp

        if init==1: return self.clique_lb(0, deg_bags) # update to find improvement 
        return self.c_lb

    cdef (int, int) clique_inspect(self, int v):
        cdef np.ndarray clique, clique_size, tmp_size
        clique, clique_size, tmp_size = self.clique, self.clique_size, self.tmp_size
        cdef int improvement = 0
        cdef int new_c, crt_size 
        new_c, crt_size = clique[v], clique_size[clique[v]]-1 # remove your self from you own clique to compair to others 
        #cdef set neighbors = self.adj_list[v]
        cdef list neighbors = list(self.adj_list[v])
        tmp_size[clique[neighbors]] = 0
        cdef int u, c, old_c
        #for u in neighbors:
        #    tmp_size[clique[u]] = 0
        for u in neighbors:
            c = clique[u]
            tmp_size[c] += 1
            if tmp_size[c] == clique_size[c] and crt_size < clique_size[c]:
                new_c = c
                crt_size = clique_size[c]

        if new_c != clique[v]:
            old_c = clique[v]
            clique_size[old_c] -= 1
            if clique_size[old_c] == 0:
                improvement = 1
            clique[v] = new_c
            clique_size[new_c] += 1
            return improvement, old_c
        return improvement, -1
