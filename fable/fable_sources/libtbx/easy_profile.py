from __future__ import division
class easy_profile(object):
  """ cProfile.Profile easy-to-use wrapper """

  def __init__(self, func, file_name, func_name, line, runs=1):
    """ Profiling of the callable func.
    file_name, func_name, line shall tell where func is defined.
    runs is the number of calls which will be performed
    """
    import cProfile
    self.prof = cProfile.Profile()
    self.func = func
    self.file_name, self.func_name, self.line = file_name, func_name, line
    self.runs = runs

  def time(self, *args, **kwds):
    """ Time spent per-call in self.func(*args, **kwds) """
    for i in xrange(self.runs):
      self.prof.runcall(self.func, *args, **kwds)
    self.prof.create_stats()
    for (file_name, line, func), data in self.prof.stats.iteritems():
      if self.file_name is not None:
        if not file_name.endswith(self.file_name): continue
      if self.func_name is not None:
        if func != self.func_name: continue
      if self.line is not None:
        if line != self.line: continue
      break
    else:
      return None
    calls = data[0]
    cumulative = data[3]
    t = cumulative/calls
    return t
