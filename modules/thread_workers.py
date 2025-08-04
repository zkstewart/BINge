from multiprocessing import Process, Pipe, Queue

class BasicProcess(Process):
    def __init__(self, *args, **kwargs):
        Process.__init__(self)
        self.args = args
        self.kwargs = kwargs
        self._exception_receiver, self._exception_sender = Pipe()
        self.exception = None
    
    def run(self):
        try:
            self.task(*self.args, **self.kwargs)
        except Exception as e:
            self._exception_sender.send(e)
    
    def task(self, *args, **kwargs):
        # Override this method in a subclass
        raise NotImplementedError
    
    def check_errors(self):
        if self._exception_receiver.poll():
            self.exception = self._exception_receiver.recv()
        if self.exception:
            raise self.exception

class ReturningProcess(BasicProcess):
    def __init__(self, *args, **kwargs):
        BasicProcess.__init__(self)
        self.args = args
        self.kwargs = kwargs
        self.queue = Queue()
        self._exception_receiver, self._exception_sender = Pipe()
        self.exception = None
    
    def run(self):
        try:
            result = self.task(*self.args, **self.kwargs)
            self.queue.put(result)
        except Exception as e:
            self._exception_sender.send(e)
            self.queue.put(None)
    
    def task(self, *args, **kwargs):
        # Override this method in a subclass
        raise NotImplementedError
    
    def get_result(self, timeout=None):
        return self.queue.get(timeout=timeout)
