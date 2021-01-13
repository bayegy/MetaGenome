from Bayegy.ampliconLibs.osEnv import OSEnv
from Bayegy.ampliconLibs.libShell import solve_conda_env

class metabat2(OSEnv):

    def __init__(self):
        super().__init__(**solve_conda_env('metabat2_home'))


class prokka(OSEnv):

    def __init__(self):
        super().__init__(**solve_conda_env('prokka_home'))
