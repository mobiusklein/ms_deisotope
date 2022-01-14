
from PyInstaller.utils.hooks import collect_submodules, collect_data_files

hiddenimports = collect_submodules("ms_deisotope._c")

datas = list(filter(lambda x: "test_data" not in x[1] and not x[0].endswith(
    '.c') and not x[0].endswith('.html') and not x[0].endswith('.pyx'), collect_data_files("ms_deisotope")))
