# @ File path
# @ Integer dummy

import os
symlink_path = path.getAbsolutePath() + "LINK"
file_path = path
os.system('cmd /c "mklink %s %s"' % (symlink_path, file_path))