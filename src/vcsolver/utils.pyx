# cython: language_level=3

from cython.view cimport array as cvarray
import numpy as np
import cython 

@cython.boundscheck(False)
@cython.wraparound(False)
def set_remove_for_set(list adj, set  vertices, int v):
    cdef int i
    for i in vertices:
        adj[i].remove(v)

cpdef void set_add_for_set(list adj, set vertices, int v):
    cdef int i
    for i in vertices: 
        adj[i].add(v)

cpdef void set_add_for_array(list adj, int[:] degrees, int[:] vertices):
    cdef int i
    for (d, v) in zip(degrees, vertices):
        adj[d].add(v)

cpdef void degrees_update_for_array(list adj, list deg_bags, int[:] degrees, list vertices):
    cdef int d
    cdef int v
    for v in vertices:
        d = degrees[v]
        deg_bags[d].discard(v)
        d = len(adj[v])
        degrees[v]=d
        deg_bags[d].add(v)

def remove_edges_for_array(list adj, list vertices):
    cdef int v,u
    cdef list old = list()
    cdef set s
    num_edges = 0
    for v in vertices:
        s = adj[v]
        old.append(s)
        num_edges += len(s)
        for u in s:
            adj[u].remove(v)
        adj[v]=set()
    return old, num_edges

cpdef void add_edges_for_array(list adj, list vertices, list old):
    cdef int v, u
    cdef int i = 0
    for v in vertices:
        adj[v].update(old[i])
        for u in adj[v]:
            adj[u].add(v)
        i += 1
        
