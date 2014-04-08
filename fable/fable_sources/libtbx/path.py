from __future__ import division
import locale
import shutil
import os
import sys
op = os.path

def norm_join(*args):
  return op.normpath(op.join(*args))

def abs_real_norm(path):
  return op.normpath(op.realpath(op.abspath(path)))

def abs_norm(path):
  return op.normpath(op.abspath(path))

def posix_relpath(
      path,
      start=".",
      enable_abspath_if_through_root=False):
  # based on relpath() in Python-2.7.2/Lib/posixpath.py
  if not path:
    raise ValueError("no path specified")
  def as_list(ap):
    result = []
    for _ in ap.split("/"):
      if (_): result.append(_)
    return result
  start_list = as_list(op.abspath(start))
  path_abs = op.abspath(path)
  path_list = as_list(path_abs)
  i = len(op.commonprefix([start_list, path_list]))
  if (i == 0 and enable_abspath_if_through_root):
    return path_abs
  rel_list = [".."] * (len(start_list)-i) + path_list[i:]
  if not rel_list:
    return "."
  return op.join(*rel_list)

def nt_is_true_abs(path):
  prefix, _ = op.splitunc(path)
  if (bool(prefix)): return True
  prefix, _ = op.splitdrive(path)
  return bool(prefix)

def nt_relpath(
      path,
      start=".",
      enable_abspath_if_through_root=False):
  # based on relpath() in Python-2.7.2/Lib/ntpath.py
  if not path:
    raise ValueError("no path specified")
  if (enable_abspath_if_through_root
        and not (    nt_is_true_abs(start)
                 and nt_is_true_abs(path))):
    raise RuntimeError(
      "nt_relpath(): both path and start must be absolute paths.")
  start_abs = abs_norm(start)
  path_abs = abs_norm(path)
  def _abspath_split(abs):
    prefix, rest = op.splitunc(abs)
    is_unc = bool(prefix)
    if not is_unc:
      prefix, rest = op.splitdrive(abs)
    rest_list = []
    for _ in rest.split("\\"):
      if _: rest_list.append(_)
    return is_unc, prefix, rest_list
  start_is_unc, start_prefix, start_list = _abspath_split(start_abs)
  path_is_unc, path_prefix, path_list = _abspath_split(path_abs)
  if path_is_unc ^ start_is_unc:
    if (enable_abspath_if_through_root):
      return path_abs
    raise ValueError("Cannot mix UNC and non-UNC paths (%s and %s)"
                              % (path, start))
  if path_prefix.lower() != start_prefix.lower():
    if (enable_abspath_if_through_root):
      return path_abs
    if path_is_unc:
      raise ValueError("path is on UNC root %s, start on UNC root %s"
                        % (path_prefix, start_prefix))
    else:
      raise ValueError("path is on drive %s, start on drive %s"
                        % (path_prefix, start_prefix))
  i = 0
  for e1, e2 in zip(start_list, path_list):
    if e1.lower() != e2.lower():
      break
    i += 1
  if (i == 0 and enable_abspath_if_through_root):
    return path_prefix + "\\" + "\\".join(path_list)
  rel_list = [".."] * (len(start_list)-i) + path_list[i:]
  if not rel_list:
    return "."
  return op.join(*rel_list)

if (os.name == "nt"):
  relpath = nt_relpath
else:
  relpath = posix_relpath

def tail_levels(path, number_of_levels):
  return op.join(*path.split(op.sep)[-number_of_levels:])

def create_target_dir(target_file):
  target_dir = op.split(target_file)[0]
  if (not op.isdir(target_dir)):
    os.makedirs(target_dir)

def move_old(path, serial_sep="_", serial_fmt="%03d"):
  if (not op.exists(path)): return
  bns = op.basename(path) + serial_sep
  dn = op.dirname(op.abspath(path))
  max_i = 0
  for ex in os.listdir(dn):
    if (ex.startswith(bns)):
      s = ex[len(bns):]
      try: i = int(s)
      except ValueError: pass
      else: max_i = max(max_i, i)
  nn = op.join(dn, bns + serial_fmt % (max_i+1))
  os.rename(path, nn)

def move_old_create_new_directory(path, serial_sep="_", serial_fmt="%03d"):
  move_old(path=path, serial_sep=serial_sep, serial_fmt=serial_fmt)
  os.makedirs(path)

def canonical_path(file_name, effective_current_working_directory=None):
  if not isinstance(file_name, (str, unicode)):
    file_name = abs(file_name)
  if (not op.isabs(file_name)):
    if (effective_current_working_directory is None):
      effective_current_working_directory = os.getcwd()
    file_name = op.join(effective_current_working_directory, file_name)
  return op.normpath(file_name)

def is_same_canoncial_file(file_names):
  assert len(file_names) == 2
  if (file_names[0] == file_names[1]): return True
  if (hasattr(op, "samefile")):
    return op.samefile(file_names[0], file_names[1])
  return False

def is_same_file(file_names, effective_current_working_directory=None):
  return is_same_canoncial_file(
    [canonical_path(file_name, effective_current_working_directory)
      for file_name in file_names])

def full_command_path(command, search_first=[], search_last=[]):
  dirs = search_first + os.environ["PATH"].split(os.pathsep) + search_last
  for path in dirs:
    path_command = op.join(path, command)
    if (op.exists(path_command)):
      return abs_norm(path=path_command)
  return None

class directory(object):

  def __init__(self, path):
    self.path = path

  def get(self, name, must_exist=True):
    assert name is not None
    result = op.join(self.path, name)
    if (must_exist and not op.exists(result)):
      from libtbx.str_utils import show_string
      raise RuntimeError("No such file or directory: %s" % show_string(result))
    return result

  def sub_directory(self, name, must_exist=True):
    result = directory(self.get(name))
    if (must_exist and not op.isdir(result.path)):
      from libtbx.str_utils import show_string
      raise RuntimeError("Not a directory: %s" % show_string(result.path))
    return result

def walk_source_tree(top, arg=None):
  def visitor(result, dirname, names):
    names_keep = []
    for name in names:
      path = op.join(dirname, name)
      if (not op.isdir(path)):
        if (not name.endswith(".pyc")):
          result.append(path)
        continue
      def is_file_in_subdir(name):
        return op.isfile(op.join(path, name))
      if (   (name == "CVS" and is_file_in_subdir("Entries"))
          or (name == ".svn"
                and (   is_file_in_subdir("README.txt")
                     or is_file_in_subdir("entries")))):
        continue
      names_keep.append(name)
    if (len(names_keep) != len(names)):
      del names[:]
      names.extend(names_keep)
  result = []
  op.walk(top, visitor, result)
  return result

def random_new_directory_name(prefix="tmp_dir_", number_of_hex_code_digits=8):
  from libtbx.utils import random_hex_code
  for i_trial in xrange(10**6):
    name = prefix + random_hex_code(number_of_digits=number_of_hex_code_digits)
    if (not op.exists(name)):
      return name
  else:
    raise AssertionError

def makedirs_race(
      path=None,
      max_trials=10,
      delay_if_exists=0.001,
      delay_after_exception=0.5):
  if (path is None):
    path = random_new_directory_name()
  import time
  for i_trial in xrange(max_trials):
    if (op.exists(path)):
      if (delay_if_exists is not None):
        # in case the OS needs time to finalize makedirs from another process
        time.sleep(delay_if_exists)
      break
    try:
      os.makedirs(path)
    except Exception:
      if (delay_after_exception is not None):
        time.sleep(delay_after_exception)
  if (not op.isdir(path)):
    raise RuntimeError("makedirs_race(%s) failure." % path)
  return path


class path_mixin(object):

  def __truediv__(self, path):
    return self.__div__(path)

  def isdir(self):
    return op.isdir(abs(self))

  def isfile(self):
    return op.isfile(abs(self))

  def exists(self):
    return op.exists(abs(self))

  def open(self, *args, **kwds):
    return open(abs(self), *args, **kwds)

  def makedirs(self):
    os.makedirs(abs(self))

  def remove(self):
    if self.exists(): os.remove(abs(self))
    assert not self.exists()

  def remove_tree(self):
    from distutils.dir_util import remove_tree
    if self.isdir():
      remove_tree(abs(self))
    else:
      self.remove()

  def listdir(self):
    return os.listdir(abs(self))

  def chmod(self, *args, **kwds):
    return os.chmod(abs(self), *args, **kwds)

  def access(self, *args, **kwds):
    return os.access(abs(self), *args, **kwds)

  def basename(self):
    return op.basename(abs(self))

  def ext(self):
    return op.splitext(self.basename())[1]

  def split(self):
    return (self.dirname(), self.basename())

  def samefile(self, other):
    if isinstance(other, str) or isinstance(other, unicode):
      return op.samefile(abs(self), other)
    else:
      return op.samefile(abs(self), abs(other))

  def is_relocatable(self):
    return isinstance(self, relocatable_path)

  def sh_value(self, anchor_var="LIBTBX_BUILD"):
    if (self.is_relocatable()):
      return op.join("$%s" % anchor_var, self.relocatable)
    return abs(self)

  def bat_value(self, anchor_var="LIBTBX_BUILD"):
    if (self.is_relocatable()):
      return op.join("%%%s%%" % anchor_var, self.relocatable)
    return abs(self)


class absolute_path(path_mixin):

  def __init__(self, path, case_sensitive=False):
    assert op.isabs(path)
    if not case_sensitive:
      path = op.normcase(path)
    path = op.normpath(path)
    self._path = path

  def reset(self, new_path, case_sensitive=False):
    self._path = absolute_path(new_path, case_sensitive)._path

  def __div__(self, other):
    return absolute_path(op.join(self._path, other))

  def __abs__(self):
    return self._path

  def __add__(self, ext):
    return absolute_path(self._path + ext)

  def __repr__(self):
    return 'absolute_path("%s")' % self._path

  def dirname(self):
    return absolute_path(op.dirname(self._path))


class relocatable_path(path_mixin):

  def __init__(self, anchor, relocatable):
    assert isinstance(anchor, absolute_path)
    self._anchor = anchor
    if op.isabs(relocatable):
      relocatable = relpath(
        path=op.realpath(abs(absolute_path(relocatable))),
        start=op.realpath(abs(self._anchor)),
        enable_abspath_if_through_root=True)
    self.relocatable = relocatable

  def anchor(self):
    return self._anchor
  anchor = property(anchor)

  def __div__(self, path):
    return relocatable_path(self._anchor, op.join(self.relocatable, path))

  def __idiv__(self, path):
    self.relocatable = op.join(self.relocatable, path)
    return self

  def __add__(self, ext):
    return relocatable_path(self._anchor, self.relocatable + ext)

  def self_or_abs_if(self, flag):
    if flag:
      return abs(self)
    else:
      return self

  def __abs__(self):
    return op.abspath(op.join(abs(self.anchor), self.relocatable))

  def __repr__(self):
    return 'relocatable_path(anchor="%s", relocatable="%s")' % (
      abs(self.anchor), self.relocatable)

  def dirname(self):
    assert self.relocatable
    return relocatable_path(self._anchor, op.dirname(self.relocatable))

  def basename(self):
    return op.basename(self.relocatable)

  def normcase(self):
    return relocatable_path(self._anchor, op.normcase(self.relocatable))

  def __eq__(self, other):
    return (    self._anchor == other._anchor
            and self.relocatable == other.relocatable)

class clean_out_directory (object) :
  """
  Utility for cleaning out Phenix (etc.) project folders, which tend to
  accumulate large amounts of temporary files and other large objects (CCP4
  maps, .geo files, etc.) which we may not want to keep around forever.
  """
  def __init__ (self,
      path_name,
      delete_kin_files=True,
      delete_geo_files=True,
      delete_map_files=True,
      delete_temp_dirs=True,
      delete_probe_files=True) :
    self.path_name = path_name
    self.n_files = self.n_dirs = self.n_bytes = 0
    self.file_paths = []
    self.dir_paths = []
    for dirname, dirnames, filenames in os.walk(path_name) :
      base_dir = os.path.basename(dirname)
      base_dir_name = base_dir.split("_")[0]
      # phenix-specific stuff
      if (base_dir_name in ["AutoSol","AutoBuild","AutoMR","LigandFit",
                            "StructureComparison",]) :
        if ("TEMP0" in dirnames) and delete_temp_dirs :
          dirnames.remove("TEMP0")
          full_path = os.path.join(dirname, "TEMP0")
          self.delete_directory(full_path)
      elif (base_dir_name in ["Refine"]) :
        if (".comm" in dirnames) and delete_temp_dirs and delete_probe_files :
          dirnames.remove(".comm")
          full_path = os.path.join(dirname, ".comm")
          self.delete_directory(full_path)
      elif (base_dir_name in ["FFT", "SuperposeMaps"]) :
        continue
      for file_name in filenames :
        full_path = os.path.join(dirname, file_name)
        if ((file_name.endswith(".kin") and delete_kin_files) or
            (file_name.endswith(".geo") and delete_geo_files) or
            (file_name.endswith(".ccp4") and delete_map_files) or
            (file_name.endswith(".xplor") and delete_map_files) or
            (file_name == "probe.txt" and delete_probe_files)) :
          self.delete_file(full_path)

  def delete_file (self, file_name) :
    if not os.path.isfile(file_name) : return
    self.n_bytes += os.path.getsize(file_name)
    self.file_paths.append(file_name)
    self.n_files += 1

  def delete_directory (self, dir_name) :
    self.n_bytes += directory_size(dir_name)
    self.dir_paths.append(dir_name)
    self.n_dirs += 1

  @property
  def n_total (self) :
    return self.n_dirs + self.n_files

  def run (self, out=sys.stdout) :
    self.show(out=out)
    print >> out, "Deleting all selected files and directories..."
    for file_name in self.file_paths :
      os.remove(file_name)
    for dir_name in self.dir_paths :
      shutil.rmtree(dir_name)

  def get_freed_space (self) :
    locale.setlocale(locale.LC_ALL, 'en_US')
    n_kb = self.n_bytes / 1000
    if (n_kb > 10000) :
      return locale.format("%.1f", n_kb / 1000, grouping=True) + " MB"
    else :
      return locale.format("%.1f", n_kb, grouping=True) + " KB"

  def show (self, out=sys.stdout) :
    if (self.n_dirs > 0) :
      print >> out, "The following %d directories will deleted:" % \
        self.n_dirs
      for dir_name in sorted(self.dir_paths) :
        print >> out, "  %s" % dir_name
    if (self.n_files > 0) :
      print >> out, "The following %d files will be deleted:" % \
        self.n_files
      for file_name in sorted(self.file_paths) :
        print >> out, "  %s" % file_name
    if (self.n_bytes > 0) :
      print >> out, "%s of disk space will be freed." % self.get_freed_space()

# http://stackoverflow.com/questions/1392413/calculating-a-directory-size-using-python
def directory_size(path):
  total_size = 0
  seen = set()
  for dirpath, dirnames, filenames in os.walk(path):
    for f in filenames:
      fp = os.path.join(dirpath, f)
      try:
        stat = os.stat(fp)
      except OSError:
        continue
      if stat.st_ino in seen:
        continue
      seen.add(stat.st_ino)
      total_size += stat.st_size
  return total_size  # size in bytes
