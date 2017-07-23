import math

def nth_root(n, k):
    assert(k > 0 and n > 0)
    r = 0
    while math.pow(r, k) < n:
        r+= 1
    return r

class Level:
    def __init__(self, num_buckets):
        self.__buckets = [[] for _ in range(num_buckets)]
        self.__base_distance = 0
        self.__first_bucket = 0
        self.__last_bucket = -1

# do not insert anything that is larger than the current base distance + max_edge into the queue
#   instead put this into a todo list
# when the last element is extracted from the queue
#   it provides a new base distance
#   this distance ensures that all todo items will fit into the buckets
# note that this scheme also works for priority queues and reduces the logarithmic factor from n*log(n) to n*log(max_edge)
class MLB:

    NumLevels = 4

    def __init__(self, max_edge):
        self.__todo_list = []
        self.__max_edge = max_edge
        self.__num_buckets = nth_root(self.__max_edge, MLB.NumLevels)
        self.__levels = [Level(self.__num_buckets) for _ in range(MLB.NumLevels)]

        print self.__num_buckets

    def extract(self):
        # - if the bottom level is empty
        #   create a new one starting to expand from the lowest non-empty level
        # - extract node from the bottom level
        #   (it is always there)
        # - if the buckets are empty now insert the todo items
        #   w.r.t. the new base distance given by the extracted node
        pass

    def insert(self, n):
        # - if the path length does not fit into the top level put it into the todo list
        #   otherwise put the node into the lowest level that fits




q = MLB(20)


