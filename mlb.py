import math

def nth_root(n, k):
    assert(k > 0 and n > 0)
    r = 0
    while math.pow(r, k) <= n:
        r+= 1
    return r

class Level:
    def __init__(self, base_offset, bucket_width, num_buckets):
        self.__buckets = [[] for _ in range(num_buckets)]
        self.__min_offsets = [0 for _ in range(num_buckets)]
        self.__bucket_width = bucket_width
        self.__num_buckets = num_buckets
        self.__base_offset = base_offset
        self.__first_bucket = num_buckets
        self.__last_bucket = 0

    def empty(self):
        return self.__last_bucket < self.__first_bucket

    def minimum(self):
        return self.__buckets[self.__first_bucket][-1]

    def set_base_weight(self, weight):
        self.__base_offset = weight / self.__bucket_width
        self.__first_bucket = self.__num_buckets
        self.__last_bucket = 0

    def insert(self, weight):
        i = weight / self.__bucket_width - self.__base_offset
        if 0 <= i and i < self.__num_buckets:
            self.__first_bucket = min(self.__first_bucket, i)
            self.__last_bucket = max(self.__last_bucket, i)
            #print ("inserted {} in bucket {} of width {}".format(weight, i, self.__bucket_width))
            bucket = self.__buckets[i]
            bucket.append(weight)
            # keep all the minima at the end
            if len(bucket) == 1 or bucket[-1] < bucket[self.__min_offsets[i]]:
                self.__min_offsets[i] = len(bucket) - 1
            elif bucket[self.__min_offsets[i]] < bucket[-1]:
                bucket[self.__min_offsets[i]], bucket[-1] = bucket[-1], bucket[self.__min_offsets[i]]
                self.__min_offsets[i] += 1
            return True
        return False

    def __update_first(self):
        while self.__first_bucket <= self.__last_bucket and len(self.__buckets[self.__first_bucket]) == 0:
            self.__first_bucket += 1

    def migrate(self, level):
        for weight in self.__buckets[self.__first_bucket]:
            ret = level.insert(weight)
            assert(ret)
        self.__buckets[self.__first_bucket] = []
        self.__min_offsets[self.__first_bucket] = 0
        self.__update_first()

    def extract(self):
        # NOTE: extract is only called on the bottom level where the weights of all elements in a bucket are equal
        assert(self.__min_offsets[self.__first_bucket] == 0)
        # NOTE: keeping all the minima at the end is only interesting for the c++ implementation
        #       because here paths are are additionally partitioned into relevant/irrelevant paths
        #       also, all of this could be avoided by making the relevant bit part of the weight (which is difficult in the current implementation)
        ret = self.__buckets[self.__first_bucket].pop()
        if self.__min_offsets[self.__first_bucket] == len(self.__buckets):
            self.__min_offsets[self.__first_bucket] = 0
        self.__update_first()
        return ret

class MLB:

    NumLevels = 4

    def __init__(self, max_edge):
        self.__base_distance = 0
        self.__last_level = MLB.NumLevels - 1
        self.__todo_list = []
        self.__max_edge = max_edge
        self.__num_buckets = nth_root(self.__max_edge, MLB.NumLevels)
        self.__levels = [Level(0, pow(self.__num_buckets, i), self.__num_buckets) for i in range(MLB.NumLevels)]

    def __update_level(self):
        while self.__levels[self.__last_level].empty() and self.__last_level < MLB.NumLevels - 1:
            self.__last_level += 1
        return self.__last_level

    def empty(self):
        return self.__levels[self.__last_level].empty()

    def extract(self):
        # if the bottom level is empty
        # create a new one starting to expand from the lowest non-empty level
        if self.__levels[0].empty():
            assert(not self.empty())
            current_level = self.__levels[self.__last_level]
            min_weight = current_level.minimum()
            while self.__last_level > 0:
                self.__last_level -= 1
                next_level = self.__levels[self.__last_level]
                next_level.set_base_weight(min_weight)
                current_level.migrate(next_level)
                current_level = next_level
        assert(not self.__levels[0].empty())

        # extract node from the bottom level
        ret = self.__levels[0].extract() + self.__base_distance
        self.__update_level()

        # if the buckets are empty now, insert the todo items
        # w.r.t. the new base distance given by the extracted node
        if self.__levels[self.__last_level].empty():
            self.__base_distance = ret
            self.__levels[-1].set_base_weight(0)
            for weight in self.__todo_list:
                assert(weight > ret)
                self.__insert(weight)
            self.__todo_list = []
        return ret

    def decrease(self, old_length, new_length):
        # it does not make much sense to implement this in the current python version
        # better do this in the c++ version where additional properties for each element are stored
        #
        # if the element is in the todo queue:
        #   (the old length can be used to determine if it is in the queue and an offset to get its position)
        #   adjust its weight
        #   if it fits into the mlb now, remove it from the queue and use __insert
        # else
        #   the offset member of the element can be used to find the element in the mlb without bucket scans
        #   determine the new insert level
        #   find the elements current level
        #   if the levels differ remove the elemnt from its current level and use Level.insert on the insert level
        #   otherwise move the element in its current bucket
        #   (note that moving is cheaper than removing because no bucket scans are necessary)
        pass

    def __insert(self, path_length):
        assert(path_length <= self.__base_distance + self.__max_edge)
        weight = path_length - self.__base_distance
        while not self.__levels[self.__last_level].insert(weight):
            self.__last_level += 1
            assert(self.__last_level < MLB.NumLevels)

    def insert(self, path_length):
        assert(path_length >= self.__base_distance)
        if path_length > self.__base_distance + self.__max_edge:
            self.__todo_list.append(path_length)
        else:
            self.__insert(path_length)

import random, timeit, heapq

def test_mlb(m, k):
    max_edge = random.randint(1, m)

    q = MLB(max_edge)
    a = []
    b = []

    for _ in range(k):
        x = random.randint(0, max_edge)
        y = x + random.randint(0, max_edge)
        q.insert(x)
        q.insert(y)
        a.append(x)
        a.append(y)

    while not q.empty():
        b.append(q.extract())

    assert(sorted(a) == b)

def test_heap(m, k):
    max_edge = random.randint(1, m)

    q = []
    a = []
    b = []

    for _ in range(k):
        x = random.randint(0, max_edge)
        y = x + random.randint(0, max_edge)
        heapq.heappush(q, x)
        heapq.heappush(q, y)
        a.append(x)
        a.append(y)

    while len(q) > 0:
        b.append(heapq.heappop(q))

    assert(sorted(a) == b)

if __name__ == '__main__':
    state = random.getstate()

    print "multi level bucket implementation: {}".format(timeit.timeit("test_mlb(10000, 1000)", "from __main__ import test_mlb", number=1000))
    random.setstate(state)

    print

    print "heap based implementation:         {}".format(timeit.timeit("test_heap(10000, 1000)", "from __main__ import test_heap", number=1000))
    random.setstate(state)

