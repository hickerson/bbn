from __future__ import division

# see also xmlrpc_server_example.py

# FIXME: rewrite ServerProxy with built-in threading for handling failed
# requests

# XXX: The ServerProxy here is based on the xmlrpclib ServerProxy
# object, but rewritten from scratch to cache requests which failed due
# to a connection error and retry them later.  I'm not sure why I can't
# just subclass it - I suspect it's a "feature" of old-style classes.
#
# XXX: Note that using the hacked ServerProxy violates the intended behavior
# of the XML-RPC protocol.  Therefore, this module allows either the original
# or modified version to be used - the original is left as the default.
#
#
# Original copyright information:
#
# XML-RPC CLIENT LIBRARY
# $Id: xmlrpclib.py 65467 2008-08-04 00:50:11Z brett.cannon $
#
# an XML-RPC client interface for Python.
#
# the marshalling and response parser code can also be used to
# implement XML-RPC servers.
#
# --------------------------------------------------------------------
# The XML-RPC client interface is
#
# Copyright (c) 1999-2002 by Secret Labs AB
# Copyright (c) 1999-2002 by Fredrik Lundh
#
# By obtaining, using, and/or copying this software and/or its
# associated documentation, you agree that you have read, understood,
# and will comply with the following terms and conditions:
#
# Permission to use, copy, modify, and distribute this software and
# its associated documentation for any purpose and without fee is
# hereby granted, provided that the above copyright notice appears in
# all copies, and that both that copyright notice and this permission
# notice appear in supporting documentation, and that the name of
# Secret Labs AB or the author not be used in advertising or publicity
# pertaining to distribution of the software without specific, written
# prior permission.
#
# SECRET LABS AB AND THE AUTHOR DISCLAIMS ALL WARRANTIES WITH REGARD
# TO THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANT-
# ABILITY AND FITNESS.  IN NO EVENT SHALL SECRET LABS AB OR THE AUTHOR
# BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR ANY
# DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
# WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
# ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE
# OF THIS SOFTWARE.
# --------------------------------------------------------------------

from libtbx import adopt_init_args
import xmlrpclib
import httplib
import socket
import subprocess
import threading
import time
import string
import random
import os
import sys

# http://stackoverflow.com/questions/372365/set-timeout-for-xmlrpclib-serverproxy
class TimeoutTransport(xmlrpclib.Transport):
  def __init__(self, timeout=socket._GLOBAL_DEFAULT_TIMEOUT, *args, **kwargs):
    xmlrpclib.Transport.__init__(self, *args, **kwargs)
    self.timeout = timeout

  def make_connection (self, host) :
    if self._connection and host == self._connection[0]:
      return self._connection[1]
    # create a HTTP connection object from a host descriptor
    chost, self._extra_headers, x509 = self.get_host_info(host)
    #store the host argument along with the connection object
    self._connection = host, httplib.HTTPConnection(chost,
      timeout=self.timeout)
    return self._connection[1]

default_timeout = 5.0
if (sys.platform == "win32") :
  default_timeout = 1.0
class ServerProxy (object) :
  def __init__(self, uri, transport=None, encoding=None, verbose=0,
               allow_none=0, use_datetime=0, timeout=default_timeout,
               raise_errors=True):
    self._pending = []
    # establish a "logical" server connection

    # get the url
    import urllib
    type, uri = urllib.splittype(uri)
    if type not in ("http", "https"):
      raise IOError, "unsupported XML-RPC protocol"
    self.__host, self.__handler = urllib.splithost(uri)
    if not self.__handler:
      self.__handler = "/RPC2"

    if transport is None:
      if type == "https":
        transport = xmlrpclib.SafeTransport(use_datetime=use_datetime)
      else:
        transport = TimeoutTransport(timeout=timeout,
          use_datetime=use_datetime)
    self.__transport = transport

    self.__encoding = encoding
    self.__verbose = verbose
    self.__allow_none = allow_none
    self._timeouts = 0
    self._errors = []
    self.raise_errors = raise_errors

  def __request(self, methodname, params):
    self._pending.append((methodname, params))
    return self.flush_requests()

  def flush_requests (self, show_timings=False) :
    t1 = time.time()
    result = None
    while len(self._pending) > 0 :
      (methodname, params) = self._pending.pop(0)
      # call a method on the remote server
      try :
        request = xmlrpclib.dumps(params, methodname,
                                  encoding=self.__encoding,
                                  allow_none=self.__allow_none)

        response = self.__transport.request(
            self.__host,
            self.__handler,
            request,
            verbose=self.__verbose
        )

        if len(response) == 1 :
          result = response[0]
      except KeyboardInterrupt :
        raise
      except Exception, e :
        msg = str(e)
        if (hasattr(e, "errno")) :
          if (e.errno in [32,54,61,104,111,10054,10061]) :
            self._pending.insert(0, (methodname, params))
            t = time.strftime("%H:%M:%S", time.localtime())
            self._errors.append("%s -- %s" % (t, msg))
            break
        if ("timed out" in msg) :
          print "XMLRPC timeout, ignoring request"
          self._timeouts += 1
        elif str(e).startswith("<ProtocolError ") :
          self._pending = []
          break
        elif ("exceptions.SystemExit" in str(e)) :
          self._pending = []
          break
        else :
          msg = "XMLRPC error: %s\nMethod: %s\nParams: %s\n" % \
            (str(e), str(methodname), ", ".join([ str(p) for p in params ]))
          if (not self.raise_errors) :
            print >> sys.stderr, msg
          else :
            raise RuntimeError(msg)
    t2 = time.time()
    if (show_timings) :
      sys.stderr.write("flush_requests: %.3fs\n" % (t2-t1))
    return result


  def __repr__(self):
    return (
          "<ServerProxy for %s%s>" %
          (self.__host, self.__handler)
          )

  __str__ = __repr__

  # note: to call a remote object with an non-standard name, use
  # result getattr(server, "strange-python-name")(args)

  def __getattr__(self, name):
    # magic method dispatcher
    return xmlrpclib._Method(self.__request, name)

  def number_of_timeout_errors (self) :
    return self._timeouts

  def get_error_messages (self) :
    return self._errors

#-----------------------------------------------------------------------
class external_program_thread (threading.Thread) :
  def __init__ (self, command_args, program_id, log=None,
      intercept_output=True) :
    adopt_init_args(self, locals())
    if self.log is None :
      self.log = sys.stdout
    threading.Thread.__init__(self)
    self._alive = True

  def run (self) :
    if self.intercept_output :
      p = subprocess.Popen(args=self.command_args, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, shell=True)
    else :
      p = subprocess.Popen(args=self.command_args, shell=True)
    while True :
      if p.poll() is not None :
        break
      else :
        time.sleep(0.5)
      if self.intercept_output :
        output = p.stdout.readline()
        if output is not None and output != "" :
          self.log.write(output)
          self.log.flush()
    self._alive = False

  # XXX: this is probably a bad idea
  def is_alive (self) :
    return self._alive

class external_program_server (object) :
  port_ranges = [ (40001, 40840),
                  (46000, 46999) ]
  def __init__ (self, command_args, program_id, timeout, cache_requests=False,
                local_port=None, log=None, intercept_output=False) :
    adopt_init_args(self, locals())
    assert isinstance(command_args, list) or isinstance(command_args, tuple)
    self._process = None
    self._server = None
    self.initialize_server()

  def initialize_server (self) :
    if self._process is None and self._server is None :
      valid_ports = []
      for (start, end) in self.port_ranges :
        valid_ports.extend([ n for n in range(start, end) ])
      i = int(random.random() * (len(valid_ports) - 1))
      self._port = valid_ports[i]
      prog_port_env = "CCTBX_%s_PORT" % string.upper(self.program_id)
      os.environ[prog_port_env] = str(self._port)
      if self.timeout is not None :
        os.environ["CCTBX_XMLRPC_TIMEOUT"] = str(self.timeout)
      if self.local_port is not None :
        os.environ["CCTBX_XMLRPC_PORT"] = str(self.local_port)
      self._process = external_program_thread(
        command_args=self.command_args,
        program_id=self.program_id,
        log=self.log,
        intercept_output=self.intercept_output)
      self._process.start()
      if self.cache_requests :
        proxy_class = ServerProxy
      else :
        proxy_class = xmlrpclib.ServerProxy
      self._server = proxy_class(uri="http://127.0.0.1:%d/RPC2" %
                                             self._port)

  def flush_requests (self, *args, **kwds) :
    if not self.cache_requests :
      return False
    elif self._server is not None :
      return self._server.flush_requests(*args, **kwds)

  def restart (self) :
    self._process = None
    self._server = None
    self.initialize_server()

  def is_alive (self) :
    if self._process is None or self._server is None :
      return False
    try :
      status = self._server.is_alive()
    except KeyboardInterrupt :
      raise
    except Exception :
      return False
    else :
      if status is None :
        return False
      else :
        return True

  def get_port (self) :
    return self._port

  def _ignore (self, *args, **kwds) :
    return True

  def __getattr__ (self, name) :
    if self._process is None or self._server is None :
      return self._ignore
    else :
      return getattr(self._server, name)

#---end
