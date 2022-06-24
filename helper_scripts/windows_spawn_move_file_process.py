#@ File file
#@ File file2
#@ File(style='directory') destination_dir
import subprocess

def start_moving_files(file_paths, destination_dir):
	for path in file_paths:
		print("Started move process of the file %s to %s" % (path, destination_dir))
		move_command = ["move", path, destination_dir]
		print(move_command)
		p = subprocess.Popen(move_command, shell=True)

file_paths = [file.getAbsolutePath(), file2.getAbsolutePath()]
destination_dir = destination_dir.getAbsolutePath()
start_moving_files(file_paths, destination_dir)