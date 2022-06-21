from ij import IJ, WindowManager

def get_transformed_PSF_from_bigstitcher_dataset(dataset_xml_path, timepoint, view):
    IJ.run("View PSFs", "select=%s process_angle=[Single angle (Select from List)] process_channel=[All channels] process_illumination=[All illuminations] process_tile=[All tiles] process_timepoint=[Single Timepoint (Select from List)] processing_angle=[angle %s] processing_timepoint=[Timepoint %s] display=[Averaged transformed PSF]" % (dataset_xml_path, view, timepoint))
    psf = WindowManager.getImage("Averaged transformed PSF")
    psf.hide()
    return psf



dataset_xml_path = "/media/tema/big_storage/work/goethe/fiji/test_dataset/DS0016_MEME6/downsampled_23_timepoint/dataset.xml"

psf = get_transformed_PSF_from_bigstitcher_dataset(dataset_xml_path, 0, 2)
psf.show()