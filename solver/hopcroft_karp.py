import numpy as np

# BFS: Finds the paths of minimal length that connect currently unmatched vertices in L to
# currently unmatched vertices in R
#
# The breadth-first search partitions the vertices of the graph into layers.
#
# The free vertices in L are used as the starting vertices of this search and form the first
# layer of the partitioning. At the first level of the search, there are only unmatched edges,
# since the free vertices in L are (by definition) not incident to any matched edges.
#
# At subsequent levels of the search, the traversed edges are required to alternate between
# matched and unmatched.
# When searching for successors from a vertex in L, only unmatched edges may be traversed,
# while from a vertex in R only matched edges may be traversed.
#
# The search terminates at the first layer k where one or more free vertices in R are reached.
#

# Finds all available shortest augmenting paths in the graph
# returns True if at least one augmented path exists, otherwise False
def hopcroft_karp_bfs(L: set, R: set, adjlist: dict, pair_L: dict, pair_R: dict, dist: dict) -> bool:
    queue = list()
    r_dummy = None

    # For each vertex in L, initialize distance to 0 if it was unmatched so far
    # We will find shortest augmented paths from each unmatched vertex
    for li in L:
        if pair_L[li] == None:  # li is a free vertex in L
            dist[li] = 0
            queue.append(li)  # enqueue the unmatched vertex in queue
        else:  # li is left endpoint of a matched edge
            dist[li] = float("inf")

    # Set distance of dummy node r_dummy to infinite. When we find the shortest path,
    # we would end at dummy node and BFS would set its distance to the length of shortest path.
    # We can use this length to eliminate any paths that are longer than that shortest path.
    # If more than one vertex has same length shortest path then for both of them we can follow the
    # dist array all the way to dummy node and when we get there we can check that length of path
    # so far is the same value as in dist[r_dummy].
    dist[r_dummy] = float("inf")  # r_dummy [sic!] is dummy vertex in L, connected to all unmatched vertices of R

    while len(queue) != 0:
        li = queue.pop(0)
        # If length of path to this node has exceeded the shortest path
        # we already found then ignore this path and move on to next node
        if dist[li] < dist[r_dummy]:
            # Go through our adjacency list of this node and
            # see if any of our neighbours for this node in R has an existing match in L
            # OR is unmatched.
            # If it has existing matching then we will travel to next node in L and repeat the process.
            # If there was no existing matching then next_L is dummy node.
            # In that case, if we are getting there for the first time then we have found
            # the shortest path and we mark dummy nodes dist as length of
            # shortest path we found. If it isn't our first time then either we have another
            # shortest path of same length or it's a longer path.
            # In either case we can just ignore because if it was another shortest path
            # then dist of dummy node is already marked correctly.
            # If it was not then dist of dummy node will be less than
            # dist[li] and it would ignored by DFS.
            for r in adjlist[li]:
                next_L = pair_R[r]
                if dist[next_L] == float("inf"):
                    dist[next_L] = dist[li] + 1

                    # Note that queue will always contain vertices from L.
                    # We don't need dist array for R.
                    # This is because from vertex in L we always go to next vertex in L
                    # (or dummy node).
                    queue.append(next_L)
    # If we found a shortest path then dist[r_dummy] would contain length of the shortest path
    return dist[r_dummy] != float("inf")


# Uses values set in dist array to traverse the path from each unmatched node in L.
# If we arrive at r_dummy and if our path length matches dist[r_dummy] then we have one of the
# shortest paths. In that case, we do symmetric difference on existing matching with shortest path.
#
# Returns True if shortest augmented path was found from li.
def hopcroft_karp_dfs(li: int, adjlist: dict, pair_L: dict, pair_R: dict, dist: dict) -> bool:
    # Recursion termination: If we arrive at dummy node during traversal then
    # we have shortest augmented path and so terminate.
    if li is not None:
        # For each neighbours of li, see if it's next node in the possible path
        for r in adjlist[li]:
            next_L = pair_R[r]
            # The neigbour is next node in path if it's matching node is our distance + 1
            if dist[next_L] == dist[li] + 1:
                # Recursively see if for this next node in path, we have shortest augmented path available
                dfs_ret = hopcroft_karp_dfs(pair_R[r], adjlist, pair_L, pair_R, dist)
                if dfs_ret == True:
                    # If so then time to do symmetric difference!
                    # Note that pair_L[li] either has previous matching or is unmatched.
                    # If it had previous matching then setting pair_L[li] to new value
                    # removes that matching and then adds a new matching
                    # which in essence is symmetric difference.
                    pair_R[r] = li
                    pair_L[li] = r
                    return True
        # Mark our node unusable for getting included in any other paths
        # so that all shortest augmented paths are vertex disjoint.
        # For bipartite case (but not for general graphs) it can be proved that
        # doing this simple vertex elimination results in maximal set of vertex disjoint paths.
        # Why not for general graphs? Imagine 5 paths horizontally and one path
        # vertical that cuts across the 5. If we choose vertices of vertical path then
        # we need to eliminate all horizontal paths resulting in non-maximal set.
        dist[li] = float("inf")
        return False
    return True


def hopcroft_karp(L: set, R: set, adjlist: dict) -> (dict, dict):
    # partial matching .. will be the final matching in the end
    # dictionary with a vertex (int) in L (or R) as a key.. the value is the matched vertex in R (or L)
    pair_L = dict()
    pair_R = dict()
    dist = dict()  # Distances used in the BFS

    for li in L:
        pair_L[li] = None
    for r in R:
        pair_R[r] = None
    #    matching = 0

    # The modified version of BFS marks shortest augmented paths available using dist array.
    # When no more paths are available it returns false.
    # Note that we need to find augmented paths starting from either only L or R because
    # these paths end on other side anyway.
    bfs_ret = hopcroft_karp_bfs(L, R, adjlist, pair_L, pair_R, dist)
    while bfs_ret == True:
        for li in L:  # for each unmatched vertex in L
            if pair_L[li] is None:
                # See if we have shortest augmented path. If we do, then
                # DFS will also mark vertices along that path unusable for
                # any other paths in next calls so that paths remain vertex disjoint.
                # Also if DFS does find shortest augmented path then it will do symmetric
                # difference along that path, setting proper pairs in pair_L and pair_R
                dfs_ret = hopcroft_karp_dfs(li, adjlist, pair_L, pair_R, dist)
        #                if dfs_ret == True:
        #                    matching += 1
        bfs_ret = hopcroft_karp_bfs(L, R, adjlist, pair_L, pair_R, dist)
    return pair_L, pair_R


def main():
    #    adjlist = [ set(), {4,5}, {4,5,6}, {5,6}, {1,2}, {1,2,3}, {2,3} ]
    #    edges_AB = [ (1,4), (1,5), (2,4), (2,5), (2,6), (3,5), (3,6) ]
    #    adjlist = [ [5,6], [5,9], [7,8], [5,9], [6,8], [0,1,3], [0,4], [2], [2,4], [1,3] ]
    adjlist = {
        0: [5, 6],
        1: [5, 9],
        2: [7, 8],
        3: [5, 9],
        4: [6, 8],
        5: [0, 1, 3],
        6: [0, 4],
        7: [2],
        8: [2, 4],
        9: [1, 3],
    }
    edges_AB = [
        (0, 5),
        (0, 6),
        (1, 5),
        (1, 9),
        (2, 7),
        (2, 8),
        (3, 5),
        (3, 9),
        (4, 6),
        (4, 8),
    ]
    edges_AB = np.array(edges_AB)
    A = set(edges_AB[:, 0])
    B = set(edges_AB[:, 1])
    print("A: ", A)
    print("B: ", B)

    pair_A, pair_B = hopcroft_karp(A, B, adjlist)
    print(pair_A, pair_B)


if __name__ == "__main__":
    main()
