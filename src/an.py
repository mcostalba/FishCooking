#!/usr/bin/python
import struct
from collections import defaultdict

S_BETA = 0
S_DEPTH = 1
S_FLAGS = 2
S_EVAL = 3
S_EVAL_MARGIN = 4
S_EXTRA1 = 5
S_EXTRA2 = 6
S_EXIT_TYPE = 7

S_FLAG_IN_CHECK = 1
S_FLAG_ALL_NODE = 2
S_FLAG_FROM_NULL = 4

def parse_stats():
  f = open('stats', 'rb')
  file_data = f.read()
  stat_struct = struct.Struct('<iBBiiiBc')
  stat_struct_size = stat_struct.size

  data = []
  for i in xrange(len(file_data) / stat_struct_size):
    d = stat_struct.unpack_from(file_data, i * stat_struct_size)
    data.append(d)

  return data

def filter_depth(stats):
  results = defaultdict(list)
  for s in stats:
    results[s[S_DEPTH]].append(s)
  return results

def main():
  stats = parse_stats()
  d = filter_depth(stats)
  for k,v in d.iteritems():
    print k, len(v)

  return

  all = []
  cut = []
  for s in stats:
    if s[S_FLAGS] & S_FLAG_ALL_NODE:
      all.append(s)
    else:
      cut.append(s) 

  print len(all), len(cut)

if __name__ == '__main__':
  main()
