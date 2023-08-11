from pysat.examples.hitman import Hitman


def min_sat_solver(sets):
    with Hitman(bootstrap_with=sets, solver="g3", htype="sorted") as hitman:
        for hs in hitman.enumerate():
            return hs
