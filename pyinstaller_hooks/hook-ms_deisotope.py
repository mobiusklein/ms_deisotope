
from PyInstaller.utils.hooks import collect_submodules

hiddenimports = collect_submodules("ms_deisotope._c")
