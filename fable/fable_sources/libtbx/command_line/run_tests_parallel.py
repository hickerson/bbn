from __future__ import division

import libtbx.test_utils.parallel
from libtbx.utils import Sorry, Usage
import libtbx.phil
import random
import os
import sys

master_phil = libtbx.phil.parse("""
directory = None
  .type = path
  .multiple = True
module = None
  .type = str
  .multiple = True
nproc = 1
  .type=  int
shuffle = False
  .type = bool
quiet = False
  .type = bool
run_in_tmp_dir = False
  .type = bool
output_junit_xml = False
  .type = bool
  .help = "Create junit-style xml output"
          "Requires junit_xml module:"
          "  https://pypi.python.org/pypi/junit-xml"
""")

def run (args) :
  if (len(args) == 0) :
    raise Usage("""libtbx.run_tests_parallel [module=NAME] [directory=path]""")
  user_phil = []
  for arg in args :
    if os.path.isdir(arg) :
      user_phil.append(libtbx.phil.parse("directory=%s" % arg))
    else :
      try :
        arg_phil = libtbx.phil.parse(arg)
      except RuntimeError :
        raise Sorry("Unrecognized argument '%s'" % arg)
      else :
        user_phil.append(arg_phil)

  params = master_phil.fetch(sources=user_phil).extract()

  if params.run_in_tmp_dir:
    import tempfile
    run_dir = tempfile.mkdtemp(suffix='', prefix='cctbxtst')
    print 'Running tests in %s' % run_dir
    os.chdir(run_dir)
  else:
    cwd = os.getcwd()
    cwd_files = os.listdir(cwd)
    if (len(cwd_files) > 0) :
      raise Sorry("Please run this program in an empty directory.")

  if params.output_junit_xml:
    try:
      import junit_xml
    except ImportError, e:
      raise Sorry(
        "Cannot import junit_xml. Try running with output_junit_xml=False")

  if (len(params.directory) == 0) and (len(params.module) == 0) :
    raise Sorry("Please specify modules and/or directories to test.")
  all_tests = []
  for dir_name in params.directory :
    if os.path.split(dir_name)[-1].find("cctbx_project")>-1:
      print 'DANGER '*10
      print 'Using the directory option in cctbx_project can be very time consuming'
      print 'DANGER '*10
    dir_tests = libtbx.test_utils.parallel.find_tests(dir_name)
    all_tests.extend(libtbx.test_utils.parallel.make_commands(dir_tests))
  for module_name in params.module :
    module_tests = libtbx.test_utils.parallel.get_module_tests(module_name)
    all_tests.extend(module_tests)
  if (len(all_tests) == 0) :
    raise Sorry("No test scripts found in %s." % params.directory)
  if (params.shuffle) :
    random.shuffle(all_tests)
  if (not params.quiet) :
    print "Running the following %d tests on %d processors:" % (len(all_tests),
      params.nproc)
    for test in all_tests :
      print "  " + test
  log = open("zlog", "wb")
  libtbx.test_utils.parallel.run_command_list(
    cmd_list=all_tests,
    nprocs=params.nproc,
    log=log,
    quiet=params.quiet,
    output_junit_xml=params.output_junit_xml)
  log.close()
  print """
============================================================================
Reminder: Please do not forget: libtbx.find_clutter
          See also: cctbx_project/libtbx/development/dev_guidelines.txt
============================================================================
"""

if (__name__ == "__main__") :
  run(sys.argv[1:])
