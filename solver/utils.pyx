from cython.view cimport array as cvarray
import numpy as np

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

def degrees_update_for_array(list adj, list deg_bags, int[:] degrees, list vertices):
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
    cpdef list old = list()
    cpdef set s
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
        

cpdef list dom_helper(int v, list adj_list):
    cdef set v_neighbors, u_neighbors
    v_neighbors = adj_list[v]
    v_neighbors.add(v)
    cdef list tmp_vc_addition = list()
    for u in v_neighbors: 
        if u == v: continue
        u_neighbors = adj_list[u]
        u_neighbors.add(u)
        if v_neighbors.issubset(u_neighbors):
            tmp_vc_addition.append(u)
        u_neighbors.remove(u)

    v_neighbors.remove(v)
    return tmp_vc_addition


