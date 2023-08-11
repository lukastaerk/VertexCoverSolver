from collections import deque
from collections.abc import Iterable


def vectorized_func(func, num_in: int, num_out: int):
    if num_in == 2 and num_out == 0:

        def f(array, value_s):
            if isinstance(value_s, Iterable):
                deque(map(func, array, value_s))
            else:
                deque(map(lambda e: func(e, value_s), array))

        return f
    elif num_in == 2 and num_out == 1:

        def f(arr1, arr2):
            return list(map(func, arr1, arr2))

        return f
    elif num_in == num_out == 1:

        def f(array):
            return list(map(func, array))

        return f
    elif num_in == 1 and num_out == 0:

        def f(array):
            return deque(map(func, array))

        return f


def graph_set_func(func, num_in: int, num_out: int):
    if num_in == 2 and num_out == 0:

        def f(adj_list, array, value_s):
            if isinstance(value_s, Iterable):
                deque(map(lambda e, v: func(adj_list[e], v), array, value_s))
            else:
                deque(map(lambda e: func(adj_list[e], value_s), array))

        return f
    elif num_in == 2 and num_out == 1:

        def f(adj_list, arr1, arr2):
            return list(map(lambda e, v: func(adj_list[e], v), arr1, arr2))

        return f
    elif num_in == num_out == 1:

        def f(adj_list, array):
            return list(map(lambda e: func(adj_list[e]), array))

        return f
    elif num_in == 1 and num_out == 0:

        def f(adj_list, array):
            return deque(map(lambda e: func(adj_list[e]), array))

        return f


set_len = graph_set_func(len, 1, 1)
set_remove = graph_set_func(set.remove, 2, 0)
set_discard = graph_set_func(set.discard, 2, 0)
set_clear = graph_set_func(set.clear, 1, 0)
set_add = graph_set_func(set.add, 2, 0)
set_update = graph_set_func(set.update, 2, 0)
set_pop = graph_set_func(set.pop, 1, 1)
execute_all = vectorized_func(lambda func: func(), 1, 0)
