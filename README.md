# TribFuse: an image processing pipeline for multiview light-sheet microsopy images of Tribolium castaneum development

A Fiji plugin for processing of light-sheet microscopy images of Tribolium development. The plugin has two branches of processing which can be run separately:
- **B-branch** - processes each view(direction) of the multiview dataset independetly. Does segmenting the embryo, cropping around the embryo, rotating images to align embryo axis with rotation axis and making maximum intensity projections along with projection motange for easy overview of the dataset. B-branch does not require beads(fiducial markers) to be present in the image, since it does not require any fusion.
- **Fusion-branch** - registeres and fuses multiple views(directions) of each timepoint with multiple fusion/deconvolution options. Embryo is also segmented and the dataste is cropped, so in the fused iamge there is minimal space around the embryo and the anterior–posterior embryonic axis is aligned with the Y-image axis.

The plugin is designed to work with raw TIFF images coming from the microscope and processing several datasets at a time.

## Installation 
If you do not have it yet, install Fiji (as shown [here](https://imagej.net/software/fiji/downloads)).
Make sure you have Fiji version at least >=2.3.0/1.53q (this is the oldest the script was tested to be working for).

Install plugin dependencies in Fiji:
- Bigstitcher
- (optional, required for GPU accelerated deconvolution) clij, clij2, clijx-assisstant, clijx-assisstant

Download this repository as a .zip archive (press on the green "Code" button, as shown [here](https://sites.northwestern.edu/researchcomputing/resources/downloading-from-github/)). Copy `tribolium_processing` folder from the archive to the `plugins` folder in your Fiji installation folder (`<your Fiji installation folder>/Fiji.app/plugins`).

### Optional Telegram notifications: 
If you would like to receive notifications via Telegram bot when the pipeline fineshes and crashes you need to do the following steps:
- copy `telegram_bot_notifications` folder somewhere on the PC you will be running the pipeline on. 
- Install `telegram-send` and `pyinstaller` via [pip](https://pypi.org/project/pip/). It is advisable to do this in a fresh [Anaconda](https://docs.anaconda.com/anaconda/user-guide/getting-started/) environment.
- Create a bot in Telegram and configure `telegram-send` to it as shown [here](https://medium.com/@robertbracco1/how-to-write-a-telegram-bot-to-send-messages-with-python-bcdf45d0a580)
- in the `telegram_bot_notifications` folder, run `pyinstaller -F send_message_to_bot.py`
- copy the resulting executable `telegram_bot_notifications/dist/send_message_to_bot.exe` to `<your Fiji installation folder>/Fiji.app/plugins/tribolium_processing/`

Example commands with Anaconda on Windows to do this:
(steps after copying `telegram_bot_notifications` folder from the archive and creating telegram bot)
Open `Anaconda Powershell Prompt` from windows menu, then
``` 
cd <path to the folder>\telegram_notifications_app\
conda create -n telegramNotif pip
pip install telegram-send
pip install -U pyinstaller
telegram-send --configure
    <Bot token that you got from Bot-father>
```
Put the password that you recieved from `telegram-send --configure` into your Telegram bot.
```
pyinstaller -F send_message_to_bot.py
copy .\dist\send_message_to_bot.exe <your Fiji installation folder>\Fiji.app\plugins\tribolium_processing\
```


### Updating the plugin
Delete `tribolium_processing` folder from `<your Fiji installation folder>/Fiji.app/plugins` folder. Download this repository as a .zip archive (press on the green "Code" button, as shown [here](https://sites.northwestern.edu/researchcomputing/resources/downloading-from-github/)). Copy `tribolium_processing` folder from the archive to the `plugins` folder in your Fiji installation folder (`<your Fiji installation folder>/Fiji.app/plugins`).



## Inputs
The plugin requires a folder with the structure as show below and a JSON file with metadata for the datasets to be processed.
Input folder strucutre:
```
├── DS0001_name
│   └── (P0)-ZStacks-Raw
│   │   ├── <Dataset specific name>TM13155919SPC000TL0000ANG000FRQ000PH00CM0CHN00.tif
│   │   ├── <Dataset specific name>TM13155919SPC000TL0001ANG000FRQ000PH00CM0CHN00.tif
│   │   ├── <Dataset specific name>TM13155919SPC001TL0000ANG000FRQ000PH00CM0CHN00.tif
│   │   ├── <Dataset specific name>TM13155919SPC001TL0001ANG000FRQ000PH00CM0CHN00.tif
|   |   ...
├── DS0002_another_name
│   └── (P0)-ZStacks-Raw
└── DS0003_another_name
    └── (P0)-ZStacks-Raw
```

Tiff files in the *(P0)-ZStacks-Raw* are separate for each timepoint and view. Speciments and views of each specimen are encoded in "SPCXXXX" part of the file name, hence for correspondace between which SPC values are which views for which sample you need to specify metadata in the JSON (see below). Timepoints are specified in "TLXXXX" part of the file name. The plugin currently is written for stacks with dimensions 1392x1040x179 and will probably not work with significatly different ones.

An example JSON metadata file can be found in the repository: metadata_TEMPLATE.json
For each folder with datasets you need to write such JSON file, where:
- **"ID"** - is the dataset ID, that should be at the beginning of the folder name (DSXXXX_...)
- **"specimens_for_directions_1234"** - is a list with correspondance of specimen numbers to directions of imaging the embryo. Direction is the index of the element in the list.
- **"head_direction"** - is the direction that the Tribolium embryo is looking at.
- **"use_manual_bounding_box"** - is wheather to use a manually specified box around an embryo for cropping the images
- **"planes_to_keep_per_direction"** - is a list where indexes are directions for which you specify manually which planes will be kept in the image. Currently you have to specify which planes to keep for all directions if you are using this parameter.
Consider using a text editor such as VSCode that supports good syntax hilighting for JSON files as well as syntax error checking.

If you have already organised files in *(P0)-ZStacks-Raw* folder from the previous run, you can you this as well. Then your files in *(P0)-ZStacks-Raw* folder should look like it is shown in the outputs.

_Optional manual crop box input:_ you can put a file with the name *manual_crop_box.roi* that will contain a ROI object saved from Fiji RoiManager that specifies a crop box around an embryo that will be used for this direction if the parameter "use_manual_bounding_box" is set to TRUE in the JSON file.
There is a small plugin to help you generate this manual crop box: *create_manual_crop_box.py*

___

# Everything below is under construction!!! (documentation only for part of the B-branch for now)

## Outputs 
The plugin will generate the following folder structure for each dataset:
```
DS0001_name
├── (B1)-Metadata
│   └── CH0001
│       ├── DR0001
│       │   ├── cropped_max_time_proj.tif
│       │   ├── crop_template.roi
│       │   ├── DatasetNamePrefix-DS0001TP(TM)DR0001CH0001PL(ZM).tif
│       │   └── Y_projected_raw_stack_for_assessment_of_plane_selection.tif
|       ...
├── (B2)-ZStacks
│   └── CH0001
│       ├── DR0001
│       │   ├── DatasetNamePrefix-DS0001TP0001DR0001CH0001PL(ZS).tif
│       │   ├── DatasetNamePrefix-DS0001TP0002DR0001CH0001PL(ZS).tif
|           ...
├── (B3)-TStacks-ZM
│   └── CH0001
│       ├── DR0001
│       │   └── DatasetNamePrefix-DS0001TP(TS)DR0001CH0001PL(ZN).tif
|       ...
├── (B4)-TStacks-ZN
│   └── CH0001
│       ├── DR0001
│       │   └── DatasetNamePrefix-DS0001TP(TS)DR0001CH0001PL(ZH).tif
|       ...
├── (B5)-TStacks-ZN-Montage
│   ├── DatasetNamePrefix-DS0001TP(TS)DR(AX)CH0001PL(ZM).tif
|   ├── DatasetNamePrefix-DS0001TP(TS)DR(AX)CH0001PL(ZH).tif
│   └── DatasetNamePrefix-DS0001TP(TS)DR(AY)CH0001PL(ZH).tif
|       ...
├── (P0)-ZStacks-Raw
|    └── CH0001
|        ├── bigstitcher_dataset
|        ├── DR0001
|        │   ├── DatasetNamePrefix-DS0001TP0001DR0001CH0001PL(ZS).tif
|        │   ├── DatasetNamePrefix-DS0001TP0002DR0001CH0001PL(ZS).tif
|            ...
└── B_BRANCH_FINISHED
```
* **_(B1)-Metadata_** folder contains files used during cropping of the image in X,Y and Z directions. 
    * DatasetNamePrefix-DS0001TP(TM)DR0001CH0001PL(ZM).tif maximum projection over Z and time dimensions simultaniously. This file is generated to get all embryo movements in one image and define a crop box around them.
    * cropped_max_time_proj.tif contains cropped maximum time projection. Used for assessment how well the plugin is cropping the images for this direction.
    * crop_template.roi ROI object from Fiji RoiManager defining the bounding box for cropping the embryo.
    * Y_projected_raw_stack_for_asessment_of_plane_selection.tif - Y projection of the stack to assess how well the plugin has cropped in Z dimension.
* **_(B2)-ZStacks_** folder contains cropped raw Z-stacks with each timepoint as a separate file.
* **_(B3)-TStacks-ZM_** folder contains maximum Z-projections of each timepoint aggregated in a stack
* **_(B4)-TStacks-ZN_** the same as _(B3)-TStacks-ZM_ but with adjusted image histogram for better viewing. Files ending with "PL(ZA)" are adjusted by thresholding the histogram and with "PL(ZH)" adjusted with histograms matched between all planes and thresholded.
* **_(B5)-TStacks-ZN-Montage_** folder contains montages adjusted and unadjusted Z-projections with directions combined in one image. Files ending with "PL(ZM)" are unadjusted, with "PL(ZA)" adjusted by thresholding the histogram and with "PL(ZH)" adjusted with histograms matched between all planes and thresholded.
* **_B_BRANCH_FINISHED_** Is generated when the plugin has finished processing this dataset. During the processing of the current dataset there will be a *B_BRANCH_ACTIVE* file, or if the plugin skipped this dataset with an error - a *B_BRANCH_ERRORED* file. If the plugin has crashed, you will see a *B_BRANCH_ACTIVE* file presist at the dataset that was being processed when the plugin has crashed.

In the directory with all datasets there will be a log file such as *2022-Mar-10-201337-b_branch.log*.

## Usage
After you click Run in a plugin editor in Fiji a window will appear. 
- You have to specify a directory with dataset folders and a JSON file with metadata for them (See **Inputs** section).
- You can specify a dataset name prefix that will be added to the beginning of all image file names.
- You have an option to copmress the images (IT IS SLOW! currently implemented as BioFormats zlib compression, which is single threaded and thus slow).
- You have an option to use previously cropped stacks if you had run the plugin on this dataset already. This will continue the plugin execution from the second to last generated cropped Z-stack.
- You can launch plugin several times with the same JSON file on the same folder with datasets and they will run in parallel. It is recommended to use separate Fiji instance for each plugin launch.

You can see output and potential errors in the output window of the script editor or in the log file in the direcotory with all datasets folders (see **Outputs**)

### Manual crop box creation
The script bounding box identification struggles with datasets where embryo flourescence is very low or only on one side. In these cases you need to create a manual crop box for each direction. You can do this as follows:
- Download create_manual_crop_box.py script the same way you did auto_B_branch_processing.py
- Open it in Fiji script editor
- Open maximum time and Z-projeciton file *DatasetNamePrefix-DS0001TP(TM)DR0001CH0001PL(ZM).tif* from the *(B1)-Metadata* folder.
- Create a bounding box using rotated rectangle tool in Fiji (can be selected by right clicking on a regular rectangle icon)
- You can run the create_manual_crop_box.py script to refine your rotated rectangle selection by rotating your bounding rectangle and assess how well the cropping worked. The script will rotate your bounding box by specified rotation angle and crop the selection for assessment
- The angle of rotation for bounding box is specified in a window that appears when you run create_manual_crop_box.py script. The rotation of the bounding box persists after the script has been run, meaning that the script rotates the selection and this rotation stays after you close the script! 
- When you are satisfied with the bounding rectangle, you need to save it as a ROI object in a file. To do this:
- Open ROI Manager in Fiji (you can type "roi manager" in Fiji search box)
    - Click *Add* 
    - Select the newly added item in the ROI Manager
    - Go to *More...* and then click *Save...*
    - Save the ROI object as a file *manual_crop_box.roi* in the *(B1)-Metadata* folder
- Repeat this for all directions

To use the manual crop boxes that you've created you need to set *"use_manual_bounding_box"* parameter to *TRUE* in the metadata JSON file as described in section **Inputs**
