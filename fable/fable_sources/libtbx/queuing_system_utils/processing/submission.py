from __future__ import division
from __future__ import with_statement

import subprocess

class Submission(object):
  """
  Handles job submissions
  """

  def __init__(self, cmds, name, out, err):

    self.cmds = cmds
    self.name = name
    self.out = out
    self.err = err


  def create(self, name):

    outfile = "%s.out" % name
    errfile = "%s.err" % name
    commands = (
      self.cmds + [ self.name, name, self.out, outfile, self.err, errfile ]
      )

    return ( execute( args = commands ), outfile, errfile )


class Synchronous(Submission):
  """
  Submits jobs synchronously
  """

  SCRIPT = \
"""\
source %s
%s << EOF
%s
EOF
"""

  def __call__(self, name, executable, script, include, cleanup):

    ( process, outfile, errfile ) = self.create( name = name )
    process.stdin.write( self.SCRIPT % ( include, executable, script ) )
    process.stdin.close()

    from libtbx.queuing_system_utils.processing import status

    return status.Synchronous(
      outfile = outfile,
      errfile = errfile,
      additional = cleanup,
      process = process,
      )


  @classmethod
  def SGE(cls, command = [ "qsub" ]):

    return cls(
      cmds = command + [ "-S", "/bin/sh", "-cwd", "-sync", "y" ],
      name = "-N",
      out = "-o",
      err = "-e",
      )


  @classmethod
  def LSF(cls, command = [ "bsub" ]):

    return cls(
      cmds = command + [ "-K" ],
      name = "-J",
      out = "-o",
      err = "-e",
      )


class SlurmStream(object):
  """
  Submits jobs synchronously to Slurm srun, using redirected file handles for
  standard streams
  """

  def __init__(self, command):

    self.cmds = command + [ "-Q" ]


  def __call__(self, name, executable, script, include, cleanup):

    process = execute( args = self.cmds + [ executable ] )
    process.stdin.write( script )
    process.stdin.close()

    from libtbx.queuing_system_utils.processing import status

    return status.SlurmStream( process = process, additional = cleanup )


class AsynchronousCmdLine(Submission):
  """
  Submits jobs asynchronously
  """

  def __init__(self, cmds, qdel, name, out, err, extract, poller, handler):

    super( AsynchronousCmdLine, self ).__init__(
      cmds = cmds,
      name = name,
      out = out,
      err = err,
      )
    self.qdel = qdel
    self.extract = extract
    self.poller = poller
    self.handler = handler


  def __call__(self, name, executable, script, include, cleanup):

    ( process, outfile, errfile ) = self.create( name = name )

    ( out, err ) = process.communicate(
      input = self.handler.script( include = include, executable = executable, script = script )
      )

    if err:
      raise RuntimeError, err

    return self.handler(
      jobid = self.extract( output = out ),
      poller = self.poller,
      outfile = outfile,
      errfile = errfile,
      additional = cleanup,
      qdel = self.qdel,
      )


  @classmethod
  def SGE(cls, poller, handler, command = [ "qsub" ]):

    return cls(
      cmds = command + [ "-S", "/bin/sh", "-cwd", "-terse" ],
      qdel = [ "qdel" ],
      name = "-N",
      out = "-o",
      err = "-e",
      extract = generic_jobid_extract,
      poller = poller,
      handler = handler,
      )


  @classmethod
  def LSF(cls, poller, command = [ "bsub" ]):

    from libtbx.queuing_system_utils.processing import status

    return cls(
      cmds = command,
      qdel = [ "bkill" ],
      name = "-J",
      out = "-o",
      err = "-e",
      extract = lsf_jobid_extract,
      poller = poller,
      handler = status.StdStreamStrategy,
      )


  @classmethod
  def PBS(cls, poller, command = [ "qsub" ]):

    from libtbx.queuing_system_utils.processing import status

    return cls(
      cmds = command + [ "-d", "." ],
      qdel = [ "qdel" ],
      name = "-N",
      out = "-o",
      err = "-e",
      extract = generic_jobid_extract,
      poller = poller,
      handler = status.StdStreamStrategy,
      )


  @classmethod
  def Slurm(cls, poller, command = [ "sbatch" ]):

    from libtbx.queuing_system_utils.processing import status

    return cls(
      cmds = command,
      qdel = [ "scancel" ],
      name = "-J",
      out = "-o",
      err = "-e",
      extract = slurm_jobid_extract,
      poller = poller,
      handler = status.SlurmStdStreamStrategy,
      )


class AsynchronousScript(object):
  """
  Submits jobs asynchronously. Options are passed using a script
  """

  CONDOR_SCRIPT = \
"""
Executable     = %s
Universe       = vanilla
initialdir = .

Output  = %s
Error   = %s
Log = %s
Log_xml = True
Notification = Never
Queue
"""

  def __init__(self, cmds, qdel, script, extract, poller, handler):

    self.cmds = cmds
    self.qdel = qdel
    self.script = script
    self.extract = extract
    self.poller = poller
    self.handler = handler


  def create(self):

    return execute( args = self.cmds )


  def __call__(self, name, executable, script, include, cleanup):

    outfile = "%s.out" % name
    errfile = "%s.err" % name
    scriptfile = "%s.condor.script" % name
    logfile = "%s.condor.xml" % name

    with open( scriptfile, "w" ) as ofile:
      ofile.write(
        self.handler.script( include = include, executable = executable, script = script )
        )

    import os
    import stat
    os.chmod( scriptfile, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR )

    process = self.create()

    ( out, err ) = process.communicate(
      input = self.script % ( scriptfile, outfile, errfile, logfile )
      )

    if err:
      raise RuntimeError, err

    return self.handler(
      jobid = self.extract( output = out ),
      poller = self.poller,
      outfile = outfile,
      errfile = errfile,
      logfile = logfile,
      additional = cleanup + [ scriptfile ],
      qdel = self.qdel,
      )

  @classmethod
  def Condor(cls, poller, command = [ "condor_submit" ], script = CONDOR_SCRIPT):

    from libtbx.queuing_system_utils.processing import status

    return cls(
      cmds = command,
      qdel = [ "condor_rm" ],
      script = script,
      extract = condor_jobid_extract,
      poller = poller,
      handler = status.LogfileStrategy,
      )


# Helpers
def execute(args):

  try:
    process = subprocess.Popen(
      args,
      stdin = subprocess.PIPE,
      stdout = subprocess.PIPE,
      stderr = subprocess.PIPE,
      )

  except OSError, e:
    raise RuntimeError, "Error: '%s': %s" % ( " ".join( args ), e )

  return process


def lsf_jobid_extract(output):

  match = LSF_JOBID_EXTRACT_REGEX().search( output )

  if not match:
    raise RuntimeError, "Unexpected response from queuing system"

  return match.group(1)


def generic_jobid_extract(output):

  return output.strip()


def condor_jobid_extract(output):

  match = CONDOR_JOBID_EXTRACT_REGEX().search( output )

  if not match:
    raise RuntimeError, "Unexpected response from queuing system"

  return match.group(1)


def slurm_jobid_extract(output):

  match = SLURM_JOBID_EXTRACT_REGEX().search( output )

  if not match:
    raise RuntimeError, "Unexpected response from queuing system"

  return match.group(1)


# Regex caching
from libtbx.queuing_system_utils.processing import util

LSF_JOBID_EXTRACT_REGEX = util.get_lazy_initialized_regex(
  pattern = r"Job <(\d+)> is submitted",
  )

CONDOR_JOBID_EXTRACT_REGEX = util.get_lazy_initialized_regex(
  pattern = r"\d+ job\(s\) submitted to cluster (\d+)\.",
  )

SLURM_JOBID_EXTRACT_REGEX = util.get_lazy_initialized_regex(
  pattern = r"Submitted batch job (\d+)"
  )
