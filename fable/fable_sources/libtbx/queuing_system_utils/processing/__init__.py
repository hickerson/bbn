"""
Generic module to provide parallel job execution on queuing systems

Provides drop-in replacement classes to those defined in the multiprocessing
module (Queue and Process), with certain restrictions placed by the pickle
module
"""
from __future__ import division
from __future__ import with_statement

import cPickle as pickle
import os
import time
import itertools
import glob
from Queue import Empty as QueueEmptyException

class InstantTimeout(object):
  """
  Timeout immediately
  """

  def delay(self, waittime):

    raise QueueEmptyException, "No data found in queue"


class TimedTimeout(object):
  """
  Timeout after given time
  """

  def __init__(self, max_delay):

      self.max_delay = max_delay


  def delay(self, waittime):

    if 0 < self.max_delay:
      waittime = min( self.max_delay, waittime )
      self.max_delay -= waittime
      time.sleep( waittime )

    else:
      raise QueueEmptyException, "No data found in queue within timeout"


class NoTimeout(object):
  """
  No timeout
  """

  def delay(self, waittime):

    time.sleep( waittime )


class Queue(object):
  """
  Queue object to receive data from jobs running on remote machines

  Data transfer is achieved via files. It is safe to use any number of
  Queue objects in the same directory, even with a matching identifier
  """

  def __init__(self, identifier, waittime = 1):

    self.waittime = waittime
    self.root = "%s_%d_%d" % ( identifier, os.getpid(), id( self ) )
    self.count = itertools.count()


  def put(self, obj):

    index = self.count.next()
    # Writing a tempfile and renaming it may prevent reading incomplete files
    tmp_name = "tmp_%s.%d" % ( self.root, index )
    assert not os.path.exists( tmp_name )

    with open( tmp_name, "wb" ) as ofile:
      pickle.dump( obj, ofile )

    target_name = "%s.%d" % ( self.root, index )
    assert not os.path.exists( target_name )
    os.rename( tmp_name, target_name )


  def get(self, block = True, timeout = None):

    if not block:
      predicate = InstantTimeout()

    else:
      if timeout is not None:
        predicate = TimedTimeout( max_delay = timeout )

      else:
        predicate = NoTimeout()

    while True:
      try:
        data = self.next()

      except StopIteration:
        predicate.delay( waittime = self.waittime )
        continue

      break

    return data


  def next(self):

    waiting = []

    for fname in glob.glob( "%s.*" % self.root ):
      ( root, ext ) = os.path.splitext( fname )

      try:
        number = int( ext[1:] )

      except ValueError:
        continue

      waiting.append( ( number, fname ) )

    if not waiting:
      raise StopIteration

    selected = min( waiting, key = lambda p: p[0] )[1]
    data = pickle.load( open( selected, "rb" ) )
    os.remove( selected )
    return data


class Job(object):
  """
  Job object to execute function calls on remote machines accessible via
  queuing systems

  Restrictions: target, args and kwargs has to be pickleable
  """

  def __init__(self, qinterface, target, args = (), kwargs = {}):

    self.qinterface = qinterface

    self.target = target
    self.args = args
    self.kwargs = kwargs

    self.status = None


  @property
  def name(self):

    return "%s_%d" % ( self.qinterface.root, id( self ) )


  @property
  def jobid (self):

    return getattr(self.status, "jobid", None)

  def start(self):

    if self.status is not None:
      raise RuntimeError, "start called second time"

    data = self.qinterface.input(
      name = self.name,
      target = self.target,
      args = self.args,
      kwargs = self.kwargs,
      )

    self.status = self.qinterface.submitter(
      name = self.name,
      executable = self.qinterface.executable,
      script = self.qinterface.script % data.script(),
      include = self.qinterface.include,
      cleanup = data.files(),
      )


  def is_alive(self):

    if self.status is None:
      raise RuntimeError, "job has not been submitted yet"

    return not self.status.is_finished()


  def join(self):

    while self.is_alive():
      time.sleep( 0.1 )

    ( stdout, stderr, self.exitcode ) = self.status.results()

    if stdout:
      print stdout

    if stderr and self.qinterface.display_stderr:
      print stderr

    if self.qinterface.save_error:
      if self.exitcode != 0:
        self.err = RuntimeError( stderr )

      else:
        self.err = None

    self.status.cleanup()


  def terminate(self):

    self.status.terminate()


  def __str__(self):

    return "%s(name = '%s')" % ( self.__class__.__name__, self.name )


class QueueHandler(object):
  """
  Handles submission requests for common queuing systems
  """
  SCRIPT = \
"""\
import cPickle as pickle
%s
target( *args, **kwargs )
"""

  def __init__(
    self,
    submitter,
    input,
    include,
    root,
    executable = "libtbx.python",
    script = SCRIPT,
    save_error = False,
    display_stderr = True,
    ):

    self.submitter = submitter
    self.root = "%s%s" % ( root, os.getpid() )
    self.input = input
    self.include = include

    self.executable = executable
    self.script = script

    self.save_error = save_error
    self.display_stderr = display_stderr


  def Job(self, target, args = (), kwargs = {}):

    return Job(
      qinterface = self,
      target = target,
      args = args,
      kwargs = kwargs,
      )


def get_libtbx_env_setpaths():

  import libtbx.load_env
  return libtbx.env.under_build( "setpaths.sh" )


def SGE(
  name = "libtbx_python",
  command = None,
  asynchronous = True,
  input = None,
  include = None,
  poller = None,
  handler = None,
  save_error = False,
  display_stderr = True,
  ):

  command = process_command_line( command = command, default = [ "qsub" ] )

  from libtbx.queuing_system_utils.processing import submission

  if asynchronous:
    if handler is None:
      from libtbx.queuing_system_utils.processing import status
      handler = status.StdStreamStrategy

    if poller is None:
      from libtbx.queuing_system_utils.processing import polling
      poller = polling.CentralPoller.SGE()

    submitter = submission.AsynchronousCmdLine.SGE(
      poller = poller,
      handler = handler,
      command = command,
      )

  else:
    submitter = submission.Synchronous.SGE( command = command )

  if input is None:
    from libtbx.queuing_system_utils.processing import transfer
    input = transfer.TemporaryFile

  if include is None:
    include = get_libtbx_env_setpaths()

  return QueueHandler(
    submitter = submitter,
    input = input,
    include = include,
    root = name,
    save_error = save_error,
    display_stderr = display_stderr,
    )


def LSF(
  name = "libtbx_python",
  command = None,
  asynchronous = True,
  input = None,
  include = None,
  poller = None,
  handler = None,
  save_error = False,
  display_stderr = True,
  ):

  command = process_command_line( command = command, default = [ "bsub" ] )

  from libtbx.queuing_system_utils.processing import submission

  if asynchronous:
    if poller is None:
      from libtbx.queuing_system_utils.processing import polling
      poller = polling.CentralPoller.LSF()

    submitter = submission.AsynchronousCmdLine.LSF(
      poller = poller,
      command = command,
      )

  else:
    submitter = submission.Synchronous.LSF( command = command )

  if input is None:
    from libtbx.queuing_system_utils.processing import transfer
    input = transfer.TemporaryFile

  if include is None:
    include = get_libtbx_env_setpaths()

  return QueueHandler(
    submitter = submitter,
    input = input,
    include = include,
    root = name,
    save_error = save_error,
    display_stderr = display_stderr,
    )


def PBS(
  name = "libtbx_python",
  command = None,
  asynchronous = True,
  input = None,
  include = None,
  poller = None,
  handler = None,
  save_error = False,
  display_stderr = True,
  ):

  command = process_command_line( command = command, default = [ "qsub" ] )

  from libtbx.queuing_system_utils.processing import submission

  if asynchronous:
    if poller is None:
      from libtbx.queuing_system_utils.processing import polling
      poller = polling.CentralPoller.PBS()

    submitter = submission.AsynchronousCmdLine.PBS(
      poller = poller,
      command = command,
      )

  else:
    raise RuntimeError, "PBS does not support synchronous submission"

  if input is None:
    from libtbx.queuing_system_utils.processing import transfer
    input = transfer.TemporaryFile

  if include is None:
    include = get_libtbx_env_setpaths()

  return QueueHandler(
    submitter = submitter,
    input = input,
    include = include,
    root = name,
    save_error = save_error,
    display_stderr = display_stderr,
    )


def Condor(
  name = "libtbx_python",
  command = None,
  asynchronous = True,
  input = None,
  include = None,
  poller = None,
  handler = None,
  save_error = False,
  display_stderr = True,
  ):

  command = process_command_line( command = command, default = [ "condor_submit" ] )

  from libtbx.queuing_system_utils.processing import submission

  if asynchronous:
    if poller is None:
      from libtbx.queuing_system_utils.processing import polling
      poller = polling.CentralPoller.Condor()

    submitter = submission.AsynchronousScript.Condor(
      poller = poller,
      command = command,
      )

  else:
    raise RuntimeError, "Condor does not support synchronous submission"

  if input is None:
    from libtbx.queuing_system_utils.processing import transfer
    input = transfer.TemporaryFile

  if include is None:
    include = get_libtbx_env_setpaths()

  return QueueHandler(
    submitter = submitter,
    input = input,
    include = include,
    root = name,
    save_error = save_error,
    display_stderr = display_stderr,
    )


def Slurm(
  name = "libtbx_python",
  command = None,
  asynchronous = True,
  input = None,
  include = None,
  poller = None,
  handler = None,
  save_error = False,
  display_stderr = True,
  ):

  from libtbx.queuing_system_utils.processing import submission

  if asynchronous:
    command = process_command_line( command = command, default = [ "sbatch" ] )

    if poller is None:
      from libtbx.queuing_system_utils.processing import polling
      poller = polling.CentralPoller.Slurm()

    submitter = submission.AsynchronousCmdLine.Slurm(
      poller = poller,
      command = command,
      )

  else:
    command = process_command_line( command = command, default = [ "srun" ] )
    submitter = submission.SlurmStream( command = command )

  if input is None:
    from libtbx.queuing_system_utils.processing import transfer
    input = transfer.TemporaryFile

  if include is None:
    include = get_libtbx_env_setpaths()

  return QueueHandler(
    submitter = submitter,
    input = input,
    include = include,
    root = name,
    save_error = save_error,
    display_stderr = display_stderr,
    )


def process_command_line(command, default):

  if command is None:
    return default

  else:
    import shlex
    return shlex.split( command )


def sge_evaluate(command):

  from libtbx.queuing_system_utils.processing import polling

  return polling.SinglePoller(
    command = process_command_line( command = command, default = [ "qstat", "-j" ] ),
    evaluate = polling.sge_single_evaluate,
    )


def lsf_evaluate(command):

  from libtbx.queuing_system_utils.processing import polling

  return polling.SinglePoller(
    command = process_command_line( command = command, default = [ "bjobs" ] ),
    evaluate = polling.lsf_single_evaluate,
    )


def pbs_evaluate(command):

  from libtbx.queuing_system_utils.processing import polling

  return polling.SinglePoller(
    command = process_command_line( command = command, default = [ "qstat", "-f" ] ),
    evaluate = polling.pbs_single_evaluate,
    )


def condor_evaluate(command):

  from libtbx.queuing_system_utils.processing import polling

  return polling.CentralPoller(
    command = process_command_line(
      command = command,
      default = [ "condor_q", "-xml" ],
      ),
    evaluate = polling.condor_xml_evaluate,
    )


def slurm_evaluate(command):

  from libtbx.queuing_system_utils.processing import polling

  return polling.SinglePoller(
    command = process_command_line(
      command = command,
      default = [ "squeue", "-o", "%t", "-h", "-j" ],
      ),
    evaluate = polling.slurm_single_evaluate,
    )


INTERFACE_FOR = {
  "sge": ( SGE, sge_evaluate ),
  "lsf": ( LSF, lsf_evaluate ),
  "pbs": ( PBS, pbs_evaluate ),
  "condor": ( Condor, condor_evaluate ),
  "slurm": ( Slurm, slurm_evaluate )
  }

def qsub (
  target,
  name="libtbx_python",
  platform="sge",
  command=None,
  polling_command=None,
  ):

  assert hasattr(target, "__call__")

  if platform not in INTERFACE_FOR:
    raise RuntimeError, "Unknown platform: %s" % platform

  ( factory, poller_factory ) = INTERFACE_FOR[ platform ]

  qinterface = factory(
      name = name,
      command = command,
      poller = poller_factory( command = polling_command ),
      asynchronous = True,
      )

  return qinterface.Job( target = target )

