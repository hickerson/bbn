from __future__ import division
from libtbx.queuing_system_utils import sge_utils, pbs_utils
from libtbx.str_utils import show_string
try: import gzip
except ImportError: gzip = None
try: import bz2
except ImportError: bz2 = None
try:
  import hashlib
  hashlib_md5 = hashlib.md5
except ImportError:
  import md5
  hashlib_md5 = md5.new
from stdlib import math
import warnings
import shutil
import glob
import time
import atexit
import traceback
import re
import sys, os
op = os.path

windows_device_names = """\
CON PRN AUX NUL COM1 COM2 COM3 COM4 COM5 COM6 COM7 COM8 COM9
LPT1 LPT2 LPT3 LPT4 LPT5 LPT6 LPT7 LPT8 LPT9""".split()

def xfrange(start, stop=None, step=None, tolerance=None):
  """A float range generator."""

  if stop is None:
    stop = start + 0.0
    start = 0.0
  else:
    start += 0.0 # force it to be a float
  if step is None:
    step = 1.0
  else:
    assert step != 0.0
  count = int(math.ceil((stop - start) / step))
  if (    tolerance is not None
      and abs(start + count * step - stop) < abs(step * tolerance)):
    count += 1
  for i in xrange(count):
    yield start + i * step

def frange(start, stop=None, step=None):
  return list(xfrange(start, stop=stop, step=step))

def xsamples(start, stop=None, step=None, tolerance=1e-6):
  return xfrange(start, stop, step, tolerance)

def samples(start, stop=None, step=None, tolerance=1e-6):
  return list(xsamples(start, stop, step, tolerance))

def escape_sh_double_quoted(s):
  "the result is supposed to be double-quoted when passed to sh"
  if (s is None): return None
  return s.replace('\\','\\\\').replace('"','\\"')

def xlen(seq):
  if (seq is None): return seq
  return len(seq)

def product(seq):
  result = None
  for val in seq:
    if (result is None):
      result = val
    else:
      result *= val
  return result

def sequence_index_dict(seq, must_be_unique=True):
  result = {}
  for i,elem in enumerate(seq):
    if (must_be_unique): assert elem not in result
    result[elem] = i
  return result

def number_from_string(string):
  # similar to libtbx.phil.number_from_value_string
  # (please review if making changes here)
  if (string.lower() in ["true", "false"]):
    raise ValueError(
      'Error interpreting "%s" as a numeric expression.' % string)
  try: return int(string)
  except KeyboardInterrupt: raise
  except Exception: pass
  try: return eval(string, math.__dict__, {})
  except KeyboardInterrupt: raise
  except Exception:
    raise ValueError(
      'Error interpreting "%s" as a numeric expression: %s' % (
        string, format_exception()))

def gzip_open(file_name, mode):
  assert mode in ["r", "rb", "w", "wb", "a", "ab"]
  if (gzip is None):
    un = ""
    if (mode[0] == "r"): un = "un"
    raise RuntimeError(
      "gzip module not available: cannot %scompress file %s"
        % (un, show_string(file_name)))
  return gzip.open(file_name, mode)

def bz2_open(file_name, mode):
  assert mode in ('r', 'w')
  if bz2 is None:
    raise RuntimeError('bz2 module not available: cannot %compress file %s'
                       % ({'r':'un', 'w':''}[mode], file_name))
  return bz2.BZ2File(file_name, mode)

def warn_if_unexpected_md5_hexdigest(
      path,
      expected_md5_hexdigests,
      hints=[],
      out=None):
  m = hashlib_md5()
  m.update("\n".join(open(path).read().splitlines()))
  current_md5_hexdigest = m.hexdigest()
  if (m.hexdigest() in expected_md5_hexdigests): return False
  warning = "Warning: unexpected md5 hexdigest:"
  file_name = "  File: %s" % show_string(path)
  new_hexdigest = "  New md5 hexdigest: %s" % m.hexdigest()
  width = max([len(s) for s in [warning, file_name, new_hexdigest]])
  if (out is None): out = sys.stdout
  print >> out, "*"*width
  print >> out, warning
  print >> out, file_name
  print >> out, new_hexdigest
  for hint in hints:
    print >> out, hint
  print >> out, "*"*width
  return True

def get_memory_from_string(mem_str):
  if type(mem_str)==type(1): return mem_str
  if type(mem_str)==type(1.): return mem_str
  mem_str = mem_str.replace(" ","").strip().upper()
  if mem_str == "": return 0
  factor=1024
  for i, greek in enumerate(["K","M","G","T","E","Z","Y"]):
    num_str=None
    if mem_str[-1]==greek:
      num_str = mem_str[:-1]
    if mem_str.find("%sB" % greek)==len(mem_str)-2:
      num_str = mem_str[:-2]
    if num_str is not None:
      try:
        num = float(num_str)
      except ValueError, e:
        raise RuntimeError("""
   The numerical portion of %s is not a valid float
""" % mem_str)
      break
    factor*=1024
  else:
    try:
      num = int(mem_str)
    except ValueError, e:
      raise RuntimeError("""
   There is no memory unit or valid float in %s
""" % mem_str)
    factor=1
  return num*factor

def getenv_bool(variable_name, default=False):
  value = os.environ.get(variable_name, None)
  if (value is None): return default
  value_lower = value.lower()
  if (value_lower not in ["false", "true", "0", "1"]):
    raise Sorry(
      'Environment variable %s must be "True", "False", "0", or "1"'
      ' (current value: "%s").' % (variable_name, value))
  return (value_lower in ["true", "1"])

def file_size(file_name):
  return os.stat(file_name).st_size

def copy_file(source, target, compress=None):
  assert op.isfile(source)
  if (op.isdir(target)):
    target = op.join(target, op.basename(source))
  if (compress is None):
    t = open(target, "wb")
  else:
    assert compress == ".gz"
    t = gzip_open(file_name=target+compress, mode="wb")
  t.write(open(source, "rb").read())
  del t

def remove_files(pattern=None, paths=None, ensure_success=True):
  assert [pattern, paths].count(None) == 1
  if (paths is None):
    paths = glob.glob(pattern)
  for path in paths:
    if (ensure_success):
      if (op.exists(path)):
        os.remove(path)
        if (op.exists(path)):
          raise RuntimeError("Cannot remove file: %s" % show_string(path))
    else:
      if (op.isfile(path)):
        os.remove(path)

def find_files (dir_name, pattern="*", files_only=True) :
  assert os.path.isdir(dir_name) and (pattern is not None)
  regex = re.compile(pattern)
  files = os.listdir(dir_name)
  matching_files = []
  for file_name in files :
    full_path = os.path.join(dir_name, file_name)
    if (files_only) and (not os.path.isfile(full_path)) :
      continue
    if (regex.search(file_name) is not None) :
      matching_files.append(full_path)
  return matching_files

def sort_files_by_mtime (file_names=None, dir_name=None, reverse=False) :
  assert ([file_names, dir_name].count(None) == 1)
  if (dir_name is not None) :
    assert os.path.isdir(dir_name)
    file_names = [ os.path.join(dir_name, fn) for fn in os.listdir(dir_name) ]
  files_and_mtimes = []
  for file_name in file_names :
    files_and_mtimes.append((file_name, os.path.getmtime(file_name)))
  files_and_mtimes.sort(lambda x,y: cmp(x[1], y[1]))
  if (reverse) :
    files_and_mtimes.reverse()
  return [ file_name for file_name, mtime in files_and_mtimes ]

def tupleize(x):
  try:
    return tuple(x)
  except KeyboardInterrupt: raise
  except Exception:
    return (x,)

def plural_s(n, suffix="s"):
  if (n == 1): return n, ""
  return n, suffix

def n_dim_index_from_one_dim(i1d, sizes):
  assert len(sizes) > 0
  result = []
  for sz in reversed(sizes):
    assert sz > 0
    result.append(i1d % sz)
    i1d //= sz
  result.reverse()
  return result

def flat_list(nested_list):
  result = []
  if (hasattr(nested_list, "__len__")):
    for sub_list in nested_list:
      result.extend(flat_list(sub_list))
  else:
    result.append(nested_list)
  return result

def select_matching(key, choices, default=None):
  for key_pattern, value in choices:
    m = re.search(key_pattern, key)
    if m is not None: return value
  return default


class Keep: pass

class Sorry(Exception):
  """
  Basic exception type for user errors; the traceback will be suppressed.
  """
  __orig_module__ = __module__
  # trick to get just "Sorry" instead of "libtbx.utils.Sorry"
  __module__ = Exception.__module__

  def reset_module (self) :
    self.__class__.__module__ = self.__class__.__orig_module__

disable_tracebacklimit = "LIBTBX_DISABLE_TRACEBACKLIMIT" in os.environ

__prev_excepthook = sys.excepthook

def sorry_excepthook(type, value, traceback):
  tb_off = (not disable_tracebacklimit and isinstance(value, Sorry))
  if (tb_off):
    class __not_set(object): pass
    prev_tracebacklimit = getattr(sys, "tracebacklimit", __not_set)
    sys.tracebacklimit = 0
  result = __prev_excepthook(type, value, traceback)
  if (tb_off):
    if (prev_tracebacklimit is __not_set):
      del sys.tracebacklimit
    else:
      sys.tracebacklimit = prev_tracebacklimit
  return result

sys.excepthook = sorry_excepthook

class Usage(Sorry):
  """
  Subclass of Sorry, for printing out usage instructions upon program
  invocation without arguments (or --help, etc.).
  """
  __module__ = Exception.__module__

class Abort(Sorry) :
  """
  Subclass of Sorry, primarily used in the Phenix GUI in response to user
  input.
  """
  __module__ = Exception.__module__

class Failure(Sorry) :
  __module__ = Exception.__module__

def detect_multiprocessing_problem():
  vers_info = sys.version_info[:2]
  if (vers_info < (2,6)):
    return "multiprocessing module not available:" \
      " Python 2.6 or higher is required" \
      " (version currently in use: %d.%d)" % vers_info
  import libtbx.load_env
  if (libtbx.env.has_module("omptbx")) :
    import omptbx
    if (omptbx.omp_version is not None) :
      return "multiprocessing is not compatible with OpenMP"
  sem_open_msg = "This platform lacks a functioning sem_open implementation"
  pool = None
  try:
    try:
      import multiprocessing
      pool = multiprocessing.Pool(processes=2)
      pool.map(func=abs, iterable=range(2), chunksize=1)
    except ImportError, e:
      if (not str(e).startswith(sem_open_msg)):
        raise
      return "multiprocessing import error: " + sem_open_msg
  finally:
    if (pool is not None):
      pool.close()
      pool.join()
  return None

def if_none(value, default):
  if (value is None): return default
  return value

def format_exception():
  ei = sys.exc_info()
  type_ = ei[0].__name__
  value = str(ei[1])
  if (value != ""):
    value = value.replace(" (<string>, line ", " (line ")
  else:
    file_name, line = traceback.extract_tb(sys.exc_info()[2], 1)[0][:2]
    if (file_name is not None):
      value = file_name+" "
    if (line is not None):
      value += "line %d" % line
  return ("%s: %s" % (type_, value)).rstrip()

def show_exception_info_if_full_testing(prefix="EXCEPTION_INFO: "):
  import libtbx.load_env
  if (    not libtbx.env.full_testing
      and not disable_tracebacklimit):
    return
  from libtbx import introspection
  from cStringIO import StringIO
  sio = StringIO()
  introspection.show_stack(out=sio)
  traceback.print_exc(file=sio)
  msg = "\n".join([prefix+line for line in sio.getvalue().splitlines()]) + "\n"
  del sio
  done = []
  for out in [sys.stdout, sys.stderr, sys.__stdout__, sys.__stderr__]:
    def is_done():
      for o in done:
        if (o is out): return True
      return False
    if (is_done()): continue
    out.write(msg)
    flush = getattr(out, "flush", None)
    if (flush is not None): flush()
    done.append(out)
  return msg

def base36_encode(integer, width=None):
  digit_set = "0123456789abcdefghijklmnopqrstuvwxyz"
  digits = []
  while (integer != 0):
    integer, i = divmod(integer, 36)
    digits.append(digit_set[i])
  if (width is not None):
    while (len(digits) < width):
      digits.append("0")
  digits.reverse()
  return "".join(digits)

def base36_timestamp(seconds_since_epoch=None, multiplier=1000, width=10):
  s = seconds_since_epoch
  if (s is None):
    s = time.time()
  return base36_encode(integer=int(s * multiplier + 0.5), width=width)

def date_and_time():
  seconds_since_epoch = time.time()
  localtime = time.localtime(seconds_since_epoch)
  if (time.daylight and localtime[8] != 0):
    tzname = time.tzname[1]
    offs = -time.altzone
  else:
    tzname = time.tzname[0]
    offs = -time.timezone
  return time.strftime("Date %Y-%m-%d Time %H:%M:%S", localtime) \
       + " %s %+03d%02d (%.2f s)" % (
           tzname, offs//3600, offs//60%60, seconds_since_epoch)


class timer_base(object):

  def __init__(self):
    self.t = self.get()

  def elapsed(self):
    t = self.get()
    d = t - self.t
    return d

  def delta(self):
    t = self.get()
    d = t - self.t
    self.t = t
    return d

  def show_elapsed(self, prefix="", out=None):
    if (out == None): out = sys.stdout
    print >> out, prefix+"%.2f s" % self.elapsed()

  def show_delta(self, prefix="", out=None):
    if (out == None): out = sys.stdout
    print >> out, prefix+"%.2f s" % self.delta()


class user_plus_sys_time(timer_base):

  def get(self):
    t = os.times()
    return t[0] + t[1]


class wall_clock_time(timer_base):
  """ motivation: when running multithreaded code, user_plus_sys_time
  would report the cumulated times for all threads: not very useful
  to analyse the scaling with the number of threads! Wall clock time, although
  it is less reliable is the only solution in that case """

  def get(self):
    return time.time()


class time_log(object):

  def __init__(self, label, use_wall_clock=False):
    self.label = label
    self.use_wall_clock = use_wall_clock
    self.accumulation = 0
    self.n = 0
    self.delta = 0
    self.timer = None

  def start(self):
    if (self.use_wall_clock):
      self.timer = wall_clock_time()
    else:
      self.timer = user_plus_sys_time()
    return self

  def stop(self):
    self.delta = self.timer.delta()
    self.timer = None
    self.accumulation += self.delta
    self.n += 1

  def average(self):
    return self.accumulation / max(1,self.n)

  def log(self):
    self.stop()
    return self.report()

  def log_elapsed(self, local_label):
    return "time_log: %s: %.2f elapsed %s" % (
      self.label, self.timer.elapsed(), local_label)

  legend = "time_log: label: n accumulation delta average"

  def report(self):
    assert self.timer is None
    return "time_log: %s: %d %.2f %.3g %.3g" % (
      self.label, self.n, self.accumulation,
      self.delta, self.average())

def human_readable_time(time_in_seconds):
  time_units = time_in_seconds
  time_unit = "seconds"
  if (time_units > 120):
    time_units /= 60
    time_unit = "minutes"
    if (time_units > 120):
      time_units /= 60
      time_unit = "hours"
      if (time_units > 48):
        time_units /= 24
        time_unit = "days"
  return time_units, time_unit

def human_readable_time_as_seconds(time_units, time_unit):
  if (isinstance(time_units, str)): time_units = float(time_units)
  if (time_unit == "seconds"): return time_units
  if (time_unit == "minutes"): return time_units*60
  if (time_unit == "hours"): return time_units*60*60
  if (time_unit == "days"): return time_units*60*60*24
  raise RuntimeError("Unknown time_unit: %s" % time_unit)

def format_timestamp_12_hour (unix_time, short=False, replace_with="unknown") :
  if unix_time is None :
    return replace_with
  elif short :
    return time.strftime("%d-%m-%y %I:%M %p", time.localtime(float(unix_time)))
  else :
    return time.strftime("%b %d %Y %I:%M %p", time.localtime(float(unix_time)))

def format_timestamp_24_hour (unix_time, short=False, replace_with="unknown") :
  if unix_time is None :
    return "unknown"
  elif short :
    return time.strftime("%d-%m-%y %H:%M", time.localtime(float(unix_time)))
  else :
    return time.strftime("%b %d %Y %H:%M", time.localtime(float(unix_time)))

format_timestamp = format_timestamp_12_hour

def format_cpu_times(show_micro_seconds_per_tick=True):
  t = os.times()
  result = "u+s,u,s: %.2f %.2f %.2f" % (t[0] + t[1], t[0], t[1])
  if (show_micro_seconds_per_tick):
    try: python_ticker = sys.gettickeraccumulation()
    except AttributeError: pass
    else:
      result += " micro-seconds/tick: %.3f" % ((t[0]+t[1])/python_ticker*1.e6)
  return result

def show_total_time(
      out=None,
      show_micro_seconds_per_bytecode_instruction=True):
  if (out == None): out = sys.stdout
  total_time = user_plus_sys_time().get()
  try: python_ticker = sys.gettickeraccumulation()
  except AttributeError: pass
  else:
    print >> out, "Time per interpreted Python bytecode instruction:",
    print >> out, "%.3f micro seconds" % (total_time / python_ticker * 1.e6)
  print >> out, "Total CPU time: %.2f %s" % human_readable_time(total_time)

def show_wall_clock_time(seconds, out=None):
  if (out is None): out = sys.stdout
  print >> out, "wall clock time:",
  if (seconds < 120):
    print >> out, "%.2f seconds" % seconds
  else:
    m = int(seconds / 60 + 1.e-6)
    s = seconds - m * 60
    print >> out, "%d minutes %.2f seconds (%.2f seconds total)" % (
      m, s, seconds)
  out_flush = getattr(out, "flush", None)
  if (out_flush is not None):
    out_flush()

class show_times:

  def __init__(self, time_start=None, out=None):
    if (time_start is None):
      t = os.times()
      self.time_start = time.time() - (t[0] + t[1])
    elif (time_start == "now"):
      self.time_start = time.time()
    else:
      self.time_start = -(0-time_start) # be sure time_start is a number
    self.out = out

  def __call__(self):
    out = self.out
    if (out is None): out = sys.stdout
    t = os.times()
    usr_plus_sys = t[0] + t[1]
    try: ticks = sys.gettickeraccumulation()
    except AttributeError: ticks = None
    s = "usr+sys time: %.2f seconds" % usr_plus_sys
    if (ticks is not None):
      s += ", ticks: %d" % ticks
      if (ticks != 0):
        s += ", micro-seconds/tick: %.3f" % (usr_plus_sys*1.e6/ticks)
    print >> out, s
    show_wall_clock_time(seconds=time.time()-self.time_start, out=out)

def show_times_at_exit(time_start=None, out=None):
  atexit.register(show_times(time_start=time_start, out=out))

class host_and_user:

  def __init__(self):
    self.host = os.environ.get("HOST")
    self.hostname = os.environ.get("HOSTNAME")
    self.computername = os.environ.get("COMPUTERNAME")
    self.hosttype = os.environ.get("HOSTTYPE")
    self.processor_architecture = os.environ.get("PROCESSOR_ARCHITECTURE")
    self.machtype = os.environ.get("MACHTYPE")
    self.ostype = os.environ.get("OSTYPE")
    self.vendor = os.environ.get("VENDOR")
    self.user = os.environ.get("USER")
    self.username = os.environ.get("USERNAME")
    self.homedir = None
    if (os.name == "nt") :
      homedrive = os.environ.get("HOMEDRIVE")
      homepath = os.environ.get("HOMEPATH")
      if (not None in [homedrive, homepath]) :
        self.homedir = os.path.join(homedrive, homepath)
    else :
      self.homedir = os.environ.get("HOME")
    getpid = getattr(os, "getpid", None)
    if (getpid is None):
      self.pid = None
    else:
      self.pid = getpid()
    self.sge_info = sge_utils.info()
    self.pbs_info = pbs_utils.chunk_info()

  def get_user_name (self) :
    if (self.user is not None) :
      return self.user
    else :
      return self.username

  def get_host_name (self) :
    if (self.host is not None) :
      return self.host
    elif (self.hostname is not None) :
      return self.hostname
    elif (self.computername is not None) :
      return self.computername
    return None

  def show(self, out=None, prefix=""):
    if (out is None): out = sys.stdout
    if (self.host is not None):
      print >> out, prefix + "HOST =", self.host
    if (    self.hostname is not None
        and self.hostname != self.host):
      print >> out, prefix + "HOSTNAME =", self.hostname
    if (    self.computername is not None
        and self.computername != self.host):
      print >> out, prefix + "COMPUTERNAME =", self.computername
    if (self.hosttype is not None):
      print >> out, prefix + "HOSTTYPE =", self.hosttype
    if (self.processor_architecture is not None):
      print >> out, prefix + "PROCESSOR_ARCHITECTURE =", \
        self.processor_architecture
    if (   self.hosttype is None
        or self.machtype is None
        or self.ostype is None
        or "-".join([self.machtype, self.ostype]) != self.hosttype):
      if (self.machtype is not None):
        print >> out, prefix + "MACHTYPE =", \
          self.machtype
      if (self.ostype is not None):
        print >> out, prefix + "OSTYPE =", \
          self.ostype
    if (self.vendor is not None and self.vendor != "unknown"):
      print >> out, prefix + "VENDOR =", \
        self.vendor
    if (self.user is not None):
      print >> out, prefix + "USER =", self.user
    if (    self.username is not None
        and self.username != self.user):
      print >> out, prefix + "USERNAME =", self.username
    if (self.pid is not None):
      print >> out, prefix + "PID =", self.pid
    self.sge_info.show(out=out, prefix=prefix)
    self.pbs_info.show(out=out, prefix=prefix)

def allow_delete_directory (target_dir) :
  """
  Check for specified reserved directories which are standard on many systems;
  these should never be deleted as part of any program.
  """
  homedir = host_and_user().homedir
  safe_dirs = [
    homedir,
    os.path.join(homedir, "Documents"),
    os.path.join(homedir, "Desktop"),
    os.path.join(homedir, "Downloads"),
    os.path.join(homedir, "Library"),
    os.path.join(homedir, "Movies"),
    os.path.join(homedir, "data"),
    "/",
    "/home",
    "/Users",
  ]
  target_dir = os.path.abspath(target_dir)
  for safe_dir in safe_dirs :
    if (target_dir == safe_dir) :
      return False
  return True

def _indentor_write_loop(write_method, indent, incomplete_line, lines):
  for line in lines:
    if (len(line) == 0):
      incomplete_line = False
    elif (incomplete_line):
      write_method(line)
      incomplete_line = False
    else:
      write_method(indent)
      write_method(line)
    write_method("\n")

class indentor(object):

  def __init__(self, file_object=None, indent="", parent=None):
    if (file_object is None):
      if (parent is None):
        file_object = sys.stdout
      else:
        file_object = parent.file_object
    self.file_object = file_object
    if (hasattr(self.file_object, "flush")):
      self.flush = self._flush
    self.indent = indent
    self.parent = parent
    self.incomplete_line = False

  def write(self, block):
    if (len(block) == 0): return
    if (block.endswith("\n")):
      _indentor_write_loop(
        write_method=self.file_object.write,
        indent=self.indent,
        incomplete_line=self.incomplete_line,
        lines=block.splitlines())
      self.incomplete_line = False
    else:
      lines = block.splitlines()
      if (len(lines) == 1):
        if (self.incomplete_line):
          self.file_object.write(lines[-1])
        else:
          self.file_object.write(self.indent + lines[-1])
      else:
        _indentor_write_loop(
          write_method=self.file_object.write,
          indent=self.indent,
          incomplete_line=self.incomplete_line,
          lines=lines[:-1])
        self.file_object.write(self.indent + lines[-1])
      self.incomplete_line = True

  def _flush(self):
    self.file_object.flush()

  def shift_right(self, indent="  "):
    return self.__class__(indent=self.indent+indent, parent=self)

class buffered_indentor(indentor):

  def __init__(self, file_object=None, indent="", parent=None):
    indentor.__init__(self, file_object, indent, parent)
    self.buffer = []

  def write(self, block):
    self.buffer.append(block)

  def write_buffer(self):
    if (self.parent is not None):
      self.parent.write_buffer()
    for block in self.buffer:
      indentor.write(self, block)
    self.buffer = []

class null_out(object):
  """Pseudo-filehandle for suppressing printed output."""

  def isatty(self): return False
  def close(self): pass
  def flush(self): pass
  def write(self, str): pass
  def writelines(self, sequence): pass

class raise_if_output(object):
  "example use: sys.stdout = raise_if_output()"

  def isatty(self): return False
  def close(self): pass
  def flush(self): pass
  def write(self, str): raise RuntimeError
  def writelines(self, sequence): raise RuntimeError

class multi_out(object):
  """
  Multiplexing output stream, e.g. for simultaneously printing to stdout
  and a logfile.
  """

  def __init__(self):
    self.labels = []
    self.file_objects = []
    self.atexit_send_to = []
    self.closed = False
    self.softspace = 0
    atexit.register(self._atexit)

  def _atexit(self):
    if (not self.closed):
      for f,a in zip(self.file_objects, self.atexit_send_to):
        if (a is not None): a.write(f.getvalue())

  def register(self, label, file_object, atexit_send_to=None):
    """Adds an output stream to the list."""
    assert not self.closed
    self.labels.append(label)
    self.file_objects.append(file_object)
    self.atexit_send_to.append(atexit_send_to)
    return self

  def replace_stringio(self,
        old_label,
        new_label,
        new_file_object,
        new_atexit_send_to=None):
    i = self.labels.index(old_label)
    old_file_object = self.file_objects[i]
    new_file_object.write(old_file_object.getvalue())
    old_file_object.close()
    self.labels[i] = new_label
    self.file_objects[i] = new_file_object
    self.atexit_send_to[i] = new_atexit_send_to

  def isatty(self):
    return False

  def close(self):
    for file_object in self.file_objects:
      if (file_object is sys.__stdout__): continue
      if (file_object is sys.__stderr__): continue
      file_object.close()
    self.closed = True

  def flush(self):
    for file_object in self.file_objects:
      flush = getattr(file_object, "flush", None)
      if (flush is not None): flush()

  def write(self, str):
    for file_object in self.file_objects:
      file_object.write(str)

  def writelines(self, sequence):
    for file_object in self.file_objects:
      file_object.writelines(sequence)

def write_this_is_auto_generated(f, file_name_generator):
  print >> f, """\
/* *****************************************************
   THIS IS AN AUTOMATICALLY GENERATED FILE. DO NOT EDIT.
   *****************************************************

   Generated by:
     %s
 */
""" % file_name_generator

class import_python_object:

  def __init__(self, import_path, error_prefix, target_must_be, where_str):
    path_elements = import_path.split(".")
    if (len(path_elements) < 2):
      raise ValueError(
        '%simport path "%s" is too short%s%s' % (
          error_prefix, import_path, target_must_be, where_str))
    module_path = ".".join(path_elements[:-1])
    try: module = __import__(module_path)
    except ImportError:
      raise ImportError("%sno module %s%s" % (
        error_prefix, module_path, where_str))
    for attr in path_elements[1:-1]:
      module = getattr(module, attr)
    try: self.object = getattr(module, path_elements[-1])
    except AttributeError:
      raise AttributeError(
        '%sobject "%s" not found in module "%s"%s' % (
          error_prefix, path_elements[-1], module_path, where_str))
    self.path_elements = path_elements
    self.module_path = module_path
    self.module = module

class input_with_prompt(object):

  def __init__(self, prompt, tracebacklimit=0):
    try: import readline
    except Exception: pass
    try: self.previous_tracebacklimit = sys.tracebacklimit
    except Exception: self.previous_tracebacklimit = None
    if (tracebacklimit is not None):
      sys.tracebacklimit = tracebacklimit
    self.input = raw_input(prompt)

  def __del__(self):
    if (self.previous_tracebacklimit is None):
      del sys.tracebacklimit
    else:
      sys.tracebacklimit = self.previous_tracebacklimit

def count_max(assert_less_than):
  i = 0
  while True:
    yield None
    i += 1
    assert i < assert_less_than

class detect_binary_file(object):

  def __init__(self, monitor_initial=None, max_fraction_non_ascii=None):
    if (monitor_initial is None):
      self.monitor_initial = 1000
    else:
      self.monitor_initial = monitor_initial
    if (max_fraction_non_ascii is None):
      self.max_fraction_non_ascii = 0.05
    else:
      self.max_fraction_non_ascii = max_fraction_non_ascii
    self.n_ascii_characters = 0
    self.n_non_ascii_characters = 0
    self.status = None

  def is_binary_file(self, block):
    if (self.monitor_initial > 0):
      for c in block:
        if (1 < ord(c) < 128):
          self.n_ascii_characters += 1
        else:
          self.n_non_ascii_characters += 1
        self.monitor_initial -= 1
        if (self.monitor_initial == 0):
          if (  self.n_non_ascii_characters
              > self.n_ascii_characters * self.max_fraction_non_ascii):
            self.status = True
          else:
            self.status = False
          break
    return self.status

  def from_initial_block(
        file_name,
        monitor_initial=None,
        max_fraction_non_ascii=None):
    detector = detect_binary_file(
      monitor_initial=monitor_initial,
      max_fraction_non_ascii=max_fraction_non_ascii)
    block = open(file_name, "rb").read(detector.monitor_initial)
    if (len(block) == 0): return False
    detector.monitor_initial = min(len(block), detector.monitor_initial)
    return detector.is_binary_file(block=block)
  from_initial_block = staticmethod(from_initial_block)

def search_for(
      pattern,
      mode,
      re_flags=0,
      lines=None,
      file_name=None):
  assert mode in ["==", "find", "startswith", "endswith", "re.search", "re.match"]
  assert [lines, file_name].count(None) == 1
  if (lines is None):
    lines = open(file_name).read().splitlines()
  result = []
  a = result.append
  if (mode == "=="):
    for l in lines:
      if (l == pattern): a(l)
  elif (mode == "startswith"):
    for l in lines:
      if (l.startswith(pattern)): a(l)
  elif (mode == "endswith"):
    for l in lines:
      if (l.endswith(pattern)): a(l)
  elif (mode == "find"):
    for l in lines:
      if (l.find(pattern) >= 0): a(l)
  elif (mode == "re.search"):
    import re
    for l in lines:
      if (re.search(pattern=pattern, string=l, flags=re_flags) is not None):
        a(l)
  else:
    import re
    for l in lines:
      if (re.match(pattern=pattern, string=l, flags=re_flags) is not None):
        a(l)
  return result


class progress_displayed_as_fraction(object):

  def __init__(self, n):
    self.n = n
    self.i = 0
    if self.n == 1: self.advance = lambda: None
    self.advance()

  def advance(self):
    if self.i > 0: sys.stdout.write('\r')
    sys.stdout.write("%i / %i" % (self.i, self.n))
    sys.stdout.flush()
    self.i += 1

  def done(self):
    if self.n == 1: return
    sys.stdout.write("\n")
    sys.stdout.flush()


class progress_bar(progress_displayed_as_fraction):

  def advance(self):
    characters = ['|']
    if self.i > 0:
      characters.extend(['=']*(self.i-1))
      characters.append('>')
    characters.extend(' '*(self.n - self.i))
    characters.append('|\r')
    sys.stdout.write(''.join(characters))
    sys.stdout.flush()
    self.i += 1

def format_float_with_standard_uncertainty(value, standard_uncertainty):
  if standard_uncertainty < 1e-16: return str(value)
  precision = -int(round(math.log10(standard_uncertainty)))
  if precision > -1:
    su = standard_uncertainty * math.pow(10, precision)
    if round(su,1) < 2:
      su *= 10
      precision += 1
    fmt_str = "%%.%if(%%i)" %precision
    return fmt_str %(value, round(su))
  else:
    precision += 1
    su = int(round(standard_uncertainty, precision))
    fmt_str = "%.0f(%i)"
    return fmt_str %(round(value, precision), su)

def random_hex_code(number_of_digits):
  import random
  digits = []
  for i_digit in xrange(number_of_digits):
    i = random.randrange(16)
    digits.append("0123456789abcdef"[i])
  return "".join(digits)

def get_svn_revision(path=None):
  # adapted from:
  #   http://code.djangoproject.com/browser/django/trunk/django/utils/version.py
  rev = None
  if path is None:
    import libtbx.load_env
    path = op.dirname(libtbx.env.dist_path(module_name="libtbx"))
  entries_path = '%s/.svn/entries' % path
  try:
    entries = open(entries_path, 'r').read()
  except IOError:
    pass
  else:
    # Versions >= 7 of the entries file are flat text.  The first line is
    # the version number. The next set of digits after 'dir' is the revision.
    if re.match('(\d+)', entries):
      rev_match = re.search('\d+\s+dir\s+(\d+)', entries)
      if rev_match:
        rev = int(rev_match.groups()[0])
  return rev

def get_build_tag(path=None):
  tag = None
  if path is None:
    import libtbx.load_env
    path = op.dirname(libtbx.env.dist_path(module_name="libtbx"))
  tag_file_path = "%s/TAG" %path
  if op.exists(tag_file_path):
    tag = open(tag_file_path).readline().strip()
  return tag

def getcwd_safe () :
  try :
    cwd = os.getcwd()
  except OSError, e :
    if (e.errno == 2) :
      raise Sorry("Could not determine the current working directory because "+
        "it has been deleted or unmounted.")
    else :
      raise e
  return cwd

def getcwd_or_default (default=None) :
  if (default is None) :
    if (os.name == "nt") :
      home_drive = os.environ.get("HOMEDRIVE", "C:")
      home_dir = os.environ.get("HOMEPATH", "\\")
      default = home_drive + home_dir
    else :
      default = os.environ.get("HOME", "/")
  try :
    cwd = os.getcwd()
  except OSError, e :
    if (e.errno == 2) :
      cwd = default
    else :
      raise e
  return cwd

def create_run_directory (prefix, default_directory_number=None) :
  """
  Create a program output directory using sequential numbering, picking the
  highest run ID.  In other words, if the prefix is 'Refine' and the current
  directory contains subdirectories named Refine_2 and Refine_9, the new
  directory will be Refine_10.
  """
  dir_number = default_directory_number
  if (dir_number is None) :
    dir_ids = []
    for file_name in os.listdir(os.getcwd()) :
      if (os.path.isdir(file_name)) and (file_name.startswith(prefix)) :
        dir_id = file_name.split("_")[-1]
        if (dir_id.isdigit()) :
          dir_ids.append(int(dir_id))
    if (len(dir_ids) > 0) :
      dir_number = max(max(dir_ids) + 1, 1)
    else :
      dir_number = 1
  dir_name = prefix + "_" + str(dir_number)
  if (os.path.isdir(dir_name)) :
    raise OSError("The directory %s already exists."%os.path.abspath(dir_name))
  else :
    os.makedirs(dir_name)
  return os.path.abspath(dir_name)

class tmp_dir_wrapper (object) :
  """
  Convenience methods for running in a (presumably empty) temporary directory
  and copying all files to another directory.  Can be used whether or not the
  temporary directory is actually defined; if None, no action will be taken.
  Otherwise, both tmp_dir and dest_dir (default is current directory) must be
  existing paths.
  """
  def __init__ (self, tmp_dir, dest_dir=None, out=sys.stdout) :
    if (dest_dir is None) :
      dest_dir = os.getcwd()
    self.tmp_dir = tmp_dir
    self.dest_dir = dest_dir
    if (tmp_dir is None) :
      pass
    elif (not os.path.isdir(tmp_dir)) :
      raise Sorry("The temporary directory %s does not exist." % tmp_dir)
    else :
      if (not os.path.isdir(dest_dir)) :
        raise Sorry("The destination directory %s does not exist." % dest_dir)
      print >> out, "Changing working directory to %s" % tmp_dir
      print >> out, "Ultimate destination is %s" % dest_dir
      os.chdir(tmp_dir)

  def transfer_files (self, out=sys.stdout) :
    if (self.tmp_dir is None) : return False
    assert os.path.isdir(self.dest_dir)
    files = os.listdir(self.tmp_dir)
    print >> out, "Copying all output files to %s" % self.dest_dir
    for file_name in files :
      print >> out, "  ... %s" % file_name
      shutil.copy(os.path.join(self.tmp_dir, file_name), self.dest_dir)
    print >> out, ""
    return True

def show_development_warning (out=sys.stdout) :
  print >> out, """
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                  WARNING - EXPERIMENTAL PROGRAM                        !!
  !!                                                                        !!
  !! This program is still in development - some functionality may be       !!
  !! missing and/or untested.  Use at your own risk!  For bug reports, etc. !!
  !! email bugs@phenix-online.org.                                          !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
"""

def check_if_output_directory_exists (file_name=None, dir_name=None) :
  if (file_name is not None) :
    assert (dir_name is None)
    dir_name = os.path.dirname(file_name)
  if (dir_name == "") : return
  if (dir_name is None) :
    raise Sorry("No output directory specified.")
  if (not op.isdir(dir_name)) :
    raise Sorry(("The specified output directory (%s) does not exist or "+
      "is not a directory.") % dir_name)
  else :
    # XXX writing to Dropbox folders is generally not a good idea
    head, tail = os.path.split(dir_name)
    while tail != "" :
      if (tail == "Dropbox") :
        warnings.warn("You are directing output to a Dropbox directory.  "+
          "Please note that this is not guaranteed to work in all cases; "+
          "use at your own risk.", UserWarning)
      head, tail = os.path.split(head)

def concatenate_python_script (out, file_name) :
  """
  Insert a Python script into an existing file, removing any __future__
  import to prevent syntax errors.  (This could be dangerous in most contexts
  but is required for some of our Coot-related scripts to work.)
  """
  data = open(file_name, "r").read()
  print >> out, ""
  print >> out, "#--- script copied from %s" % os.path.basename(file_name)
  for line in data.splitlines() :
    if line.startswith("from __future__") :
      continue
    else :
      print >> out, line
  print >> out, "#--- end"
  print >> out, ""

def greek_time(secs):
  for greek in ["","milli", "micro", "nano"]:
    if secs>1:
      break
    secs*=1000
  return secs, greek

###########################
# URL retrieval functions #
###########################

libtbx_urllib_proxy = None

def install_urllib_http_proxy (server, port=80, user=None, password=None) :
  global libtbx_urllib_proxy
  import urllib2
  if (user is None) :
    proxy = urllib2.ProxyHandler({'http': '%s:%d' % (server, port) })
    opener = urllib2.build_opener(proxy)
  else :
    proxy = urllib2.ProxyHandler({
      'http': 'http://%s:%s@%s:%s' % (user, password, server, port),
    })
    auth = urllib2.HTTPBasicAuthHandler()
    opener = urllib2.build_opener(proxy, auth, urllib2.HTTPHandler)
  libtbx_urllib_proxy = proxy
  urllib2.install_opener(opener)
  print "Installed urllib2 proxy at %s:%d" % (server, port)
  return proxy

def urlopen (*args, **kwds) :
  """
  Substitute for urllib2.urlopen, with automatic HTTP proxy configuration
  if specific environment variables are defined.
  """
  if ("CCTBX_HTTP_PROXY" in os.environ) and (libtbx_urllib_proxy is None) :
    server = os.environ["CCTBX_HTTP_PROXY_SERVER"]
    port = os.environ.get("CCTBX_HTTP_PROXY_PORT", 80)
    user = os.environ.get("CCTBX_HTTP_PROXY_USER", None)
    passwd = os.environ.get("CCTBX_HTTP_PROXY_PASSWORD", None)
    if (user is not None) and (password is None) :
      raise Sorry("You have defined a user name for the HTTP proxy, but "+
        "no password was specified.  Please set the environment variable "+
        "CCTBX_HTTP_PROXY_PASSWORD.")
    install_urllib_http_proxy(
      server=server,
      port=port,
      user=user,
      password=password)
  import urllib2
  return urllib2.urlopen(*args, **kwds)

class download_progress (object) :
  """
  Simple proxy for displaying download status - here with methods for
  writing to the console, but can be subclassed and used for graphical display.
  """
  def __init__ (self, log=None, n_kb_total=None) :
    if (log is None) :
      log = null_out()
    self.log = log
    self.n_kb_total = n_kb_total
    self.n_kb_elapsed = 0

  def set_total_size (self, n_kb_total) :
    self.n_kb_total = n_kb_total
    self.n_kb_elapsed = 0

  def increment (self, n_kb) :
    assert (self.n_kb_total is not None)
    self.n_kb_elapsed += n_kb
    return self.show_progress()

  def show_progress (self) :
    self.log.write("\r%d/%d KB downloaded" % (self.n_kb_elapsed,
      self.n_kb_total))
    self.log.flush()
    return True

  def percent_finished (self) :
    assert (self.n_kb_total is not None)
    return 100 * min(1.0, self.n_kb_elapsed / self.n_kb_total)

  def complete (self) :
    self.log.write("\rDownload complete")

  def run_continuously (self) :
    """
    Placeholder for cases where the download is not being run asynchronously.
    """
    pass

class download_target (object) :
  """
  Flexible callable object for retrieving a file from a URL, with optional
  HTTPS authentication.  Designed to be runnable in a separate thread with
  graphical progress update.

  Note that in some circumstances SSL support may be missing from the socket
  module, in which case we use 'curl' to download securely.  (This will not
  work on Windows, obviously.)
  """
  def __init__ (self,
      url,
      file_name,
      use_curl=None, # SSL only
      user=None, # SSL only
      password=None, # SSL only
      base_url=None) :  # SSL only
    self.url = url
    self.file_name = file_name
    self.use_curl = use_curl
    self.user = user
    self.password = password
    self.base_url = base_url
    if (not None in [self.user, self.password]) :
      assert (self.base_url is not None)
      import socket
      if ((not self.use_curl) and (hasattr(socket, "ssl")) and
          (hasattr(socket.ssl, "__call__"))) :
        self.use_curl = False
      else :
        self.use_curl = True

  def __call__ (self, log=None, progress_meter=None) :
    if (log is None) :
      log = null_out()
    if (progress_meter is None) :
      progress_meter = download_progress(log=log)
    from libtbx import easy_run
    import urllib2
    file_name = self.file_name # return value
    if (not self.use_curl) :
      if (not None in [self.user, self.password]) :
        passman = urllib2.HTTPPasswordMgrWithDefaultRealm()
        passman.add_password(None, self.base_url, params.user, params.password)
        authhandler = urllib2.HTTPBasicAuthHandler(passman)
        opener = urllib2.build_opener(authhandler)
        urllib2.install_opener(opener)
      req = urllib2.urlopen(self.url)
      info = req.info()
      n_kb_total = int(info['Content-length']) / 1024
      progress_meter.set_total_size(n_kb_total)
      # TODO adjust chunk size automatically based on download speed
      n_kb_chunk = getattr(self, "n_kb_chunk", 512)
      chunksize = n_kb_chunk * 1024
      fp = open(self.file_name, 'wb')
      while True:
        chunk = req.read(chunksize)
        if not chunk: break
        if not progress_meter.increment(n_kb_chunk) :
          file_name = None
          break
        fp.write(chunk)
      fp.close()
      progress_meter.complete()
    else :
      progress_meter.run_continuously()
      if (not None in [self.user, self.password]) :
        curl_args = "--user %s:%s" % (self.user, self.password)
      rc = easy_run.call("curl %s \"%s\" -o %s" % (curl_args, self.url,
        self.file_name))
      progress_meter.complete()
      if (rc != 0) :
        raise RuntimeError("curl exited with code %d" % rc)
    if (file_name is None) :
      return None
    return op.abspath(self.file_name)
