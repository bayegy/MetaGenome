import os
import sys


class OSEnv:
    """docstring for OSEnv"""

    def __init__(self, **kwargs):
        self.kwargs = {k.upper(): v for k, v in kwargs.items()}

    def __enter__(self):
        for env, val in self.kwargs.items():
            origin_env = self.get_env(env)
            setattr(self, env, origin_env)
            if env == "PYTHONPATH":
                pp = val.split(':')
                self.len_pp = len(pp)
                sys.path = pp + sys.path
            os.environ[env] = ':'.join([val, origin_env]).strip(":")
            print("Current {} is {}".format(env, self.get_env(env)))
        return self

    def __exit__(self, type, value, trace):
        for env, val in self.kwargs.items():
            if env == "PYTHONPATH":
                del sys.path[:self.len_pp]
            origin_env = getattr(self, env, "")
            if origin_env:
                os.environ[env] = origin_env
            else:
                os.environ.pop(env)
        if trace:
            print(trace)

    def get_env(self, key):
        try:
            return os.environ[key]
        except KeyError:
            return ""
