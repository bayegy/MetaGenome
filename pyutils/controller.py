import time
import os
from .read import trunc_fa
from MetaGenome.pipconfig import settings
from multiprocessing import cpu_count
from concurrent.futures import ThreadPoolExecutor


def wait_sge(first_check=1):
    time.sleep(first_check * 60)
    while True:
        if not list(os.popen("qstat -q {} 2>&1".format(settings.sge_queue))):
            break
        time.sleep(60)


def all_path_exists(paths):
    for path in paths:
        if not os.path.exists(path):
            return False
    return True


def keep_only(path):
    directory, file = os.path.split(path)
    for fl in os.listdir(directory):
        if not file == fl:
            cmd = "rm -r {}".format(os.path.join(directory, fl))
            print(cmd)
            os.system(cmd)


def sge_parallel(func):
    def wfunc(self, fa, out_dir, cat_file=None, first_check=10, splits=10, tmp_dir=None, remove_tmp=False, **kwargs):
        if splits == 1:
            kwargs.update(threads=int(cpu_count() * 0.9))
            func(self, fa, out_dir, escape_sge=True, **kwargs)
            # wait_sge(first_check)
            return
        fa = os.path.abspath(fa)
        out_dir = os.path.abspath(out_dir)
        tmp_dir = tmp_dir or os.path.join(out_dir, func.__name__ + '_tmpdir')
        # dirname = os.path.dirname(fa)
        basename = os.path.basename(fa)
        if os.path.exists(tmp_dir):
            print("Tmp_dir already exists, presume you have trunk the fasta and skip trunk step.")
        else:
            trunc_fa(fa, splits, out_dir=tmp_dir)
        print("######################Running " + str(func))
        start_time = time.time()

        if not len(os.listdir(tmp_dir)) == splits:
            raise Exception(
                "Number of trunk is not equal to the splits,\
                 change splits or remove tmp_dir and try again.")

        for i in range(splits):
            trunk_out = os.path.join(tmp_dir, "{basename}.trunk{i}".format(**locals()))
            fa_trunk = os.path.join(trunk_out, basename)
            if not os.path.exists(fa_trunk):
                raise Exception("Looks like the tmp_dir is not correct, remove tmp_dir and try again.")
            if not os.path.exists(os.path.join(trunk_out, "done")):
                keep_only(fa_trunk)
                func(self, fa=fa_trunk, out_dir=trunk_out, **kwargs)

        wait_sge(first_check)
        number_done = len(list(os.popen("ls {tmp_dir}/*/done".format(tmp_dir=tmp_dir))))
        if not number_done == splits:
            raise Exception("Not all workers exit correctly, simply run again will solve this problem.")
        args = locals()
        del args['self']
        if cat_file:
            self.system("cat {tmp_dir}/{basename}.trunk*/{cat_file} > {out_dir}/{cat_file}", **args)
        if remove_tmp and os.path.exists(os.path.join(out_dir, cat_file)):
            self.system("rm -r {tmp_dir}", **args)
        end_time = time.time()
        time_used = (end_time - start_time) / 60
        print("######################" + str(func) + " done; time used: {} min".format(time_used))
    return wfunc


def default_parallel(func):
    def wfunc(self, fa, out_dir, cat_file=None, first_check=10, splits=10, tmp_dir=None, remove_tmp=True, **kwargs):
        # print("decorator did nothing")
        kwargs.update(threads=int(cpu_count() * 0.9))
        return func(self, fa, out_dir, escape_sge=True, **kwargs)
    return wfunc


def sge_decorator(func):
    """
    对于被修饰的函数：fq_list参数值应是单个fq文件的路径（或路径对）
    对于最终函数：fq_list参数是多个fq文件路径的二(至少)维数组(list)
    """

    def wfunc(self, fq_list, first_check=1, pass_if_exists=[], clean_before=[], **kwargs):
        print("######################Running " + str(func))
        start_time = time.time()
        for fq in fq_list:
            parsed_fqs = self.parse_fq_list(fq)
            if pass_if_exists:
                paths = [p.format(**parsed_fqs, **self.context) for p in pass_if_exists]
                if all_path_exists(paths):
                    continue
            if clean_before:
                for path in clean_before:
                    path = path.format(**parsed_fqs, **self.context)
                    if os.path.exists(path):
                        self.system("rm -r {}".format(path))
            func(self, fq_list=fq, **kwargs)
        wait_sge(first_check)
        end_time = time.time()
        time_used = (end_time - start_time) / 60
        print("######################" + str(func) + " done; time used: {} min".format(time_used))
    return wfunc


def pool_decorator(func):
    """
    对于被修饰的函数：fq_list参数值应是单个fq文件的路径（或路径对）
    对于最终函数：fq_list参数是多个fq文件路径的二(至少)维数组(list)
    """

    def wfunc(self, fq_list, first_check=1, pass_if_exists=[], clean_before=[], **kwargs):
        """
        pass_if_exists: 路径列表， 如果列表中的路径都存在，的跳过相应样本的处理，路径中的{sample}及MetagenomePipline的context将会被相应属性替代
        """
        print("######################Running " + str(func))
        max_workers = int(cpu_count() // int(kwargs.get('threads')))
        start_time = time.time()
        executor = ThreadPoolExecutor(max_workers=max_workers)
        for fq in fq_list:
            parsed_fqs = self.parse_fq_list(fq)
            if pass_if_exists:
                paths = [p.format(**parsed_fqs, **self.context) for p in pass_if_exists]
                if all_path_exists(paths):
                    continue
            if clean_before:
                for path in clean_before:
                    path = path.format(**parsed_fqs, **self.context)
                    if os.path.exists(path):
                        self.system("rm -r {}".format(path))
            executor.submit(func, self, fq_list=fq, **kwargs)
        executor.shutdown(True)
        end_time = time.time()
        time_used = (end_time - start_time) / 60
        print("######################" + str(func) + " done; time used: {} min".format(time_used))
    return wfunc
