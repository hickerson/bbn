from __future__ import division
def run(args):
  import time
  show = args + [str(time.time())]
  print " ".join(show)

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
