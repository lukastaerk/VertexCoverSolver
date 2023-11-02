
def verifier(vc, edges):
    for edge in edges:
        if edge[0] not in vc and edge[1] not in vc:
            print("Edge not covered: ", edge)
            return False
    return True