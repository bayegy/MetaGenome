import os


class SGEQueue(object):
    """docstring for SGEQueue"""

    def __init__(self, name):
        self.name = name

    @property
    def info_dict(self):
        tmp_dict = {}
        for line in list(os.popen("qconf -sq {}".format(self.name))):
            key, *values = line.split()
            tmp_dict[key] = values
        return tmp_dict

    @property
    def this_host(self):
        return list(os.popen("hostname"))[0].strip()

    @property
    def this_resources(self):
        return dict(
            h_vmem=self.h_vmem.get(self.this_host),
            slots=self.slots.get(self.this_host),
            processors=self.processors.get(self.this_host),
        )

    def __getattr__(self, name):
        values = self.info_dict.get(name)
        if not values:
            raise Exception("SGE queue does not have this attribute")
        default, *specials = values[0].split(',')
        special_dict = {}
        for special in specials:
            key, value = special.split('=')
            special_dict[key.lstrip('[')] = value.rstrip(']')
        return {
            host: (special_dict.get(host) or default) for host in self.info_dict['hostlist']
        }
