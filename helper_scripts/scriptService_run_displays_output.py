#@ ScriptService scriptService
#@ CommandService command
#@ ModuleService module
from ij import IJ, ImagePlus, ImageStack, WindowManager
from java.io import File

Arguments = ["imp", WindowManager.getCurrentImage(), "projected_dimension", "Z", "projection_type", "Max"]
impInput = scriptService.run(
    File("/home/tema/work/for_Mashunja/fiji_scripts/helper_scripts/project_imagej2.py"),
    True, Arguments)
