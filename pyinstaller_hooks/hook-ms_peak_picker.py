
from PyInstaller.utils.hooks import collect_submodules

hiddenimports = collect_submodules("ms_peak_picker._c")
