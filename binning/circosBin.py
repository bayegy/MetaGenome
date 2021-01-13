from Bayegy.ampliconLibs.systemMixin import SystemMixin
from MetaGenome.pyutils.read import iter_fa
import pandas as pd


class CircosBin(SystemMixin):
    """docstring for CircosBin"""

    def __init__(self, bin_fa, top=20):
        kwargs = locals()
        kwargs = {k: v for k, v in kwargs.items() if k not in ["self"]}
        self.set_attr(**kwargs)

    def parse_bin_fa(self):
        seqs_dict = {}
        kar_ids = []
        lens = []
        for header, seq in iter_fa(bp, trim_line_break=True):
            kar_id = header.split()[0].lstrip(">")
            seqs_dict[kar_id] = seq
            kar_ids.append(kar_id)
            lens.append(len(seq))

        df = pd.DataFrame({"Length": lens}, index=kar_ids)

    def plot_circos(self):
