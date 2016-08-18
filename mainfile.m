%-------------------------------------------------------------------------
%-------- ASSIGN DEFAULT PARAMETERS %%% ADVANCED USERS ONLY --------------
%-------------------------------------------------------------------------
clear; global param; param = struct;
param.analysisStage = 'Setting Defaults';
% ARRAY ORIENTATION
param.AO.scale_umperpx = 5; % um per pix scale (req: number, >0)
param.AO.wellSpacing_px = 80; %px (req: positive integer)
param.AO.laneSpacing_px = 200; %px (req: positive integer)
param.AO.minYCropOffset_px = 30;% px, minimum distance required on top and bottom of array (req: integer)
param.AO.multiDirectionROI = 0; % 0 = No, 1 = Yes: Create additional ROI up and down for all devices (req: 0 or 1)
param.AO.refArrayTemplate=''; % full file path to reference array template (optional)

% CURVE FITTING
param.CF.base_folder=''; %base folder for curve fitting, e.g. [param.AO.base_folder,'/prep/']
param.CF.lane2process=1; %1 or 2, sets the number of lanes to be processed
param.CF.flankSubtraciton_px=2; % pixel width along the side of each lane for the background subtraciton steps
param.CF.spacerHeight_px=0; 
param.CF.spacerWidth_px=0;
param.CF.subtraction='Average'; % subtraction mode
param.CF.loop=1; % 1 or 0, loop through all files or only run first file
param.CF.counter=1; % 1 or 0, whether to display counter during curve fitting (may slow fitting)
param.CF.overlayArray=0; %1 or 0, controls whether overlay array will be produced for images
%options related to manual selection:
param.CF.refManSelect=1; %1 or 0, controls whether or not to run manual selection in curve fitting (goodmanualselect)
param.CF.columnPerView=25; % req positive integer, controls coumns displayed in each image in goodmanualselect
param.CF.rowPerView=5; % req positive integer, controls rows displayed in each image in goodmanualselect
param.CF.defaultManualSelect='off'; % 'off' or 'on': toggles starting value in goodmanualselect func.

% DATA ANALYSIS
% not yet written...

%% -----------------------------------------------------------------------
%-------------------- BEGINNING GUI CODE ---------------------------------
%-------------------------------------------------------------------------
% Array Orientation GUI 'UI_ArrayExtraction.m'
uiwait(UI_ArrayExtraction);
%   imports defaults from global variable 'param.AO'
%   generates UI_ArrayExtraction.fig
%   users follow the UI to do the following:
%       -select files to anlayze
%       -assign array and image templates
%       -set operating parameters of array extraciton
%       -set nicknames for files (e.g. protein or condition abreviation)
%   Exports parameters in the global variable 'param'

%% Array Orientation Execution 'ArrayExtactionExecution'
ArrayExtracitonExecution
%   imports defaults from global variable 'param.AO'
%   users select 'wells features' on the slide to do the following:
%       -rotate the slide
%       -establish the ROI array
%   follows array or image templates (defined in UI_ArrayExtraction) for
%       ease of use.
%   exports rotated,croped array images and .mat files (containing array
%   details) into the subfolder of the main directory (directory/prep/)

%% Curve fit GUI
uiwait(UI_Operational);
%   imports defaults saved in global variable 'param.CF'
%   generates UI_Operational.fig
%   imports .mat files from selected (or imported 'param.filePrepPath')
%   User do the following
%       -set which files will perform ROI extract for
%       -set which files will perform CurveFit (cannot be done w/o ROI
%           subtract)
%       -Set the Ref ROI array (either loaded externally or set in table)
%       -Set curve fitting parameters
%           Subtraction mode: 'Average', 'Axial
%           FlankSubtraction, Spacer Height, Spacer Width, Loop all, Export
%           overlayed images.

%% Curve Fitting Execution
CurveFittingExecution
%   imports defaults saved in global variable 'param.CF'
%   generates CurveFittingExecution.fig
%   imports .mat files from selected path ('param.CF.base_folder')
%   set's ROI reference
%   User do the following
%       -set which files will perform ROI extract for
%       -set which files will perform CurveFit (cannot be done w/o ROI
%           subtract)
%       -Set the Ref ROI array (either loaded externally or set in table)
%       -Set curve fitting parameters
%           Subtraction mode: 'Average', 'Axial
%           FlankSubtraction, Spacer Height, Spacer Width, Loop all, Export
%           overlayed images.