from __future__ import division

class MultiQueue(object):

  def __init__(self):

    self.queue_for = {}


  def create(self, name):

    import Queue
    self.queue_for[ name ] = Queue.Queue()


  def remove(self, name):

    del self.queue_for[ name ]


  def put(self, name, message):

    self.queue_for[ name ].put( message )


  def get(self, name, block = True, timeout = None):

    return self.queue_for[ name ].get( block, timeout )


  def get_nowait(self, name):

    return self.queue_for[ name ].get_nowait()


class Queue(object):

  def __init__(self, server):

    import socket
    import os
    self.identifier = "%s-%s-%s" % ( socket.getfqdn(), os.getpid(), id( self ) )

    self.server = server
    self.server.multiqueue.create( self.identifier )


  def shutdown(self):

    self.server.multiqueue.remove( self.identifier )
    self.server = None


  def put(self, value):

    self.server.multiqueue.put( self.identifier, value )


  def get(self, block = True, timeout = None):

    return self.server.multiqueue.get(
      self.identifier,
      block = block,
      timeout = timeout,
      )


  def get_nowait(self):

    return self.server.multiqueue.get_nowait( self.identifier )


class Manager(object):

  def __init__(self, manager):

    self.manager = manager
    self.multiqueue = manager.get() # caching


  def Queue(self):

    return Queue( server = self )


  def __getstate__(self):

    result = self.__dict__.copy()
    result[ "address" ] = self.manager.address
    result[ "authkey" ] = str( self.manager._authkey )
    del result[ "manager" ]
    del result[ "multiqueue" ]
    return result


  def __setstate__(self, result):

    manager = self.get_client_manager(
      address = result[ "address" ],
      authkey = result[ "authkey" ],
      )
    result[ "manager" ] = manager
    result[ "multiqueue" ] = manager.get()
    del result[ "address" ]
    del result[ "authkey" ]
    self.__dict__ = result


  @classmethod
  def Server(cls, port = 0, keylength = 16):

    from multiprocessing.managers import BaseManager

    class QManager(BaseManager):

      pass

    multiqueue = MultiQueue()
    QManager.register( "get", lambda: multiqueue )

    import socket
    import string
    import random

    manager = QManager(
      address = ( socket.getfqdn(), port ),
      authkey = "".join(
        random.choice( string.ascii_letters ) for i in range( keylength )
        ),
      )
    manager.start()

    return cls( manager = manager )


  @classmethod
  def Client(cls, server, port, authkey):

    return cls(
      manager = cls.get_client_manager(
        address = (server, port ),
        autkey = authkey,
        )
      )


  @staticmethod
  def get_client_manager(address, authkey):

    from multiprocessing.managers import BaseManager

    class QManager(BaseManager):

      pass

    QManager.register( "get" )

    manager = QManager( address = address, authkey = authkey )
    manager.connect()

    return manager

