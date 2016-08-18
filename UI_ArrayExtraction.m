%%% ARRAY ORIENTATION SETUP GRAPHIC USER INTERFACE %%%

%{ 
---------------------------------------------------------------------------
----------------------------GUI Insturctions-------------------------------
---------------------------------------------------------------------------
Instructions for the Array Orientation Setup GUI
Written by Todd A. Duncombe
    facilitated by MATLAB's GUIDE (gui construction interface)
Updated: 12/4/15

RECOMENDATIONS:
-We recommend analyzing a single slide at a time using this analysis code.
In which case, the best organization practice is creating a separate folder
for each unique slide to be analyzed.
    i.e. a slide folder may contain multiple proteins or antibodies on the
    same slide, either from multichannel fluorescense scan,
    or stripping and reprobing.
    
    This is not required. If analyzing multiple slides in parallel turn 
    'Array Match' off.

BEFORE USE:
-Organize the microarray images of interest into a new directory. (tif required)
    Often, different images will represent different proteins or antibodies on
    this same slide.

-Visually inspect the images (e.g. in imageJ) that the array section of 
    interest is visable and shared among all slides. Do not crop the array
    seciton 'tight', it should have >1 mm of space extra on all sides.

(1) 'SELECT DIRECTORY' (top left)
Select the directory containing the .tif files of interest. The path for
the selected directory will be displayed near the headure of the GUI. All
.tif files from the directory will be loaded in the file list box.

(2) Move files from 'FILE LIST' (left) to the 'FILES TO BE PROCESSED LIST' (right)
-Click on the files of interest to move them to files to be processed.
-Click 'ADD ALL' to quickly move all files.
-To remove files from 'files to be processed', simply click on them again 
in the files to be processed box.

(3) 'ADD FILES TO TABLE' (center of gui)
Once all files to be processed are in the correct box, click 'ADD FILES TO
TABLE'. This will update the file assignment table with these files. The
image sizes will be automatically loaded into the 'image size' column. If
buttons 'Array Match ON' or 'Auto Image Match ON', are on they will be
automatically performed (described below).

(4) Type the file nicknames into the 'NAME' column on the table
Fill in file nicknames for all images (e.g. GFP, aGFP, ab-GFP, anti-GFP). 
There must be no duplicates. This facilitates easy/intuitive referencing of the
images & corresponding data later on. It's best to choose a notation and 
stay consitent between slides. (example below)

(5) set 'ARRAY-MATCH'
Array match sets the ROI constant for all images. It is recommended when you
want to compare protein expression from different markers in the same cell.
-If all images are of the same slide, set ARRAY-MATCH to '1'.
-If there are multiple slides being analyzed, matching slides slides should
have a matching numeric value. Set non array matching images as = 0 (example below)

(6) set 'IMAGE-MATCH'
Image Match identifies images that match in image dimensions and naming, to
reduce the number of alignment steps required. Each matching image set will
share a positive integer. Non image matching images = 0(example below)

File Assignment Table Example:
FileName:                           Name:     Array Match:    Image Match:
date_slide1_scan1_gfp_bTub_488       gfp            1               1
date_slide1_scan1_gfp_bTub_555       bTUB           1               1
date_slide1_scan2_HER2_ERK_488       HER2           1               2
date_slide1_scan2_HER2_ERK_555       ERK            1               2
date_slide1_scan3_mTOR_488           mTOR           1               0
date_slide1_scan4_aACT_488           aACT           1               0

--- if multiple slides are being processed simultanously (not recommened)--
date_slide2_scan1_gfp_bTub_488       slide2_gfp     2               3
date_slide2_scan1_gfp_bTub_555       slide2_bTUB    2               3
date_slide3_scan1_gfp_488            slide3_gfp     0               0

(7) Set Image/Array Settings
'SCALE' - resolution of the the image, in microns per pixel
'WELL TO WELL LATERAL DISTANCE' - transverse pixel spacing between wells
'LANE LENGTH' - axial pixel spacing between wells

Note: 'LAST ACTION' (bottom of gui) is updated for (almost) all steps, this is the 
place to check if the operation performance of the GUI is perplexing.

OPTIONAL - ADVANCED USERS ONLY
(*) 'REFERENCE ARRAY'
Imports an array from a previously processed file for Array # '1' only.
Useful if adding data to a slide that was previously analyzed for other proteins.

(*) 'ADVANCED SETTINGS'
'MULTIDIRECTIONAL ROI GENERATION' - Creates a double ROI for each image (up
and down). Useful for multiple single cell fractionations.
'MINIMUM Y CROP OFFSET' - Minimum y crop offset around the array.

---------------------------------------------------------------------------
-----------------DATA OUTPUT: 'param.AO' (GLOBAL STRUCT VARIABLE)-------------
---------------------------------------------------------------------------
param.AO.analysisStage
    a string updated in each stage of the code (debugging tool,
    e.g. this .m file is the 'Array Orientation Setup GUI' stage)
param.AO.base_folder
    complete folder path for files to be processed  (req: string)
param.AO.file 
    filename list (req: string cell array 1xN, with file extension)
param.AO.name
    file nickname list (req: string cell array 1xN, with file extension)
param.AO.arrayTemplate
    array matching list for consistent ROI assignment for the same slide (req: positive integer cell array 1xN)
param.AO.transformTemplate
    image matching list to eliminate duplicate alignment steps (req: positive integer cell array 1xN)
param.AO.scale_umperpx
    resolution of image, um per pix scale (req: number, >0)
param.AO.wellSpacing_px
    well to well spacing, px (req: positive integer)
param.AO.laneSpacing_px
    lane to lane spacing, px (req: positive integer)
param.AO.refArrayTemplate
    complete path for a reference array, or '' for blank (req: string)
param.AO.minYcropOffset_um 
    minimum distance required on top and bottom of array, um (req: number)
    (150 um is the default, changing it is not recommended)
param.AO.multiDirectionROI 
    0 = No, 1 = Yes: Create additional ROI up and down for all devices
    (req: 0 or 1) (for performing bidirecitonal analysis,
    e.g. nuclear & cytoplasmic fracitons)

param.AO variable is printed in the command line prior to exiting the GUI
(when run is selected)
%}
%------------------------------------------------------------------------
%-------------------UNUSED FUNCTIONS - DO NOT MODIFY---------------------
%---------------REQUIRED FOR GUIDE MATLAB GUI CONSTRUCTOR----------------
function varargout = UI_ArrayExtraction(varargin)

% UI_ARRAYEXTRACTION MATLAB code for UI_ArrayExtraction.fig
%      UI_ARRAYEXTRACTION, by itself, creates a new UI_ARRAYEXTRACTION or raises the existing
%      singleton*.
%
%      H = UI_ARRAYEXTRACTION returns the handle to a new UI_ARRAYEXTRACTION or the handle to
%      the existing singleton*.
%
%      UI_ARRAYEXTRACTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UI_ARRAYEXTRACTION.M with the given input arguments.
%
%      UI_ARRAYEXTRACTION('Property','Value',...) creates a new UI_ARRAYEXTRACTION or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UI_ArrayExtraction_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UI_ArrayExtraction_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UI_ArrayExtraction

% Last Modified by GUIDE v2.5 15-Dec-2015 15:24:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UI_ArrayExtraction_OpeningFcn, ...
                   'gui_OutputFcn',  @UI_ArrayExtraction_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.

% Update handles structure
guidata(handles.figure1, handles);

% --- Outputs from this function are returned to the command line.
function varargout = UI_ArrayExtraction_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function lane_px_CreateFcn(hObject,eventdata,handles)
function well_px_CreateFcn(hObject,eventdata,handles)
function scale_CreateFcn(hObject,eventdata,handles)
function files2process_box_CreateFcn(hObject,eventdata,handles)
function filelistbox_CreateFcn(hObject,eventdata,handles)
function arrayMatchON_CreateFcn(hObject,eventdata,handles)
function tabletransfer_CreateFcn(hObject,eventdata,handles)
function assignmentTable_CreateFcn(hObject,eventdata,handles)

%------------------------------------------------------------------------
%-------------------OPENING FUNCTION - (sets defaults)-------------------
%------------------------------------------------------------------------
% Executes just before UI_ArrayExtraction is visible. Assigns defaults.
function UI_ArrayExtraction_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for UI_ArrayExtraction
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

%%% SET YOUR DEFAULT SETTINGS (if not already defined) %%%
global param
% if param.AO variable has not been created, add default settings
if ~exist('param','var') 
    param = struct;
    param.AO.scale_umperpx = 5; % um per pix scale
    param.AO.wellSpacing_px = 80; %um
    param.AO.laneSpacing_px = 200; %um
    param.AO.minYCropOffset_px = 30;% um, minimum distance required on top and bottom of array
    param.AO.multiDirectionROI = 0; % 0 = No, 1 = Yes: Create additional ROI up and down for all devices
    param.AO.refArrayTemplate=''; % full file path toe reference array template (optional)
end
param.analysisStage = 'Array Orientation Setup GUI';

%assign defaults
set(handles.scale,'String',num2str(param.AO.scale_umperpx));
set(handles.well_px,'String',num2str(param.AO.wellSpacing_px));
set(handles.lane_px,'String',num2str(param.AO.laneSpacing_px));
set(handles.wellSpacing_um,'String',num2str(param.AO.wellSpacing_px*param.AO.scale_umperpx));
set(handles.laneSpacing_um,'String',num2str(param.AO.laneSpacing_px*param.AO.scale_umperpx));
handles.param.AO.minYCropOffset_px = param.AO.minYCropOffset_px*param.AO.scale_umperpx;
handles.param.AO.multiDirectionROI = param.AO.multiDirectionROI;

% table columns (DO NOT MODIFY)
handles.table.fileNamesFull=1;
handles.table.imageSize=2;
handles.table.name=3;
handles.table.arrayMatch=4;
handles.table.imageMatch=5;
set(handles.assignmentTable,'Data',cell(6,length(fieldnames(handles.table))));

initialize_gui(hObject, handles, false);
set(hObject,'Name', 'Array Orientation Setup GUI');

%------------------------------------------------------------------------
%--------------------------CLOSING FUNCTION------------------------------
%------------------------------------------------------------------------
% Executes when 'RUN' button pressed
% -extracts & checks param.AOeters, stores them in global variable 'param.AO'
% -closes gui if param.AOeters are valid
function execute_Callback(hObject, eventdata, handles)
% hObject    handle to execute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

base_folder = get(handles.directory,'string');
if isequal(base_folder,'selected directory...')
    set(handles.lastAction, 'string', 'Cannot run, no directroy selected')
    return;
end

tabledata = get(handles.assignmentTable,'data');
%Check table has been filled out
if size(tabledata,2) < 4 || isempty(tabledata)
    set(handles.lastAction, 'string', 'Cannot run, table not completed')
    return;
end

fileNamesFull = tabledata(:,handles.table.fileNamesFull)';

fileNamesShort = tabledata(:,handles.table.name)';
%Check all names filled
if any(cellfun(@isempty,fileNamesShort))
    set(handles.lastAction, 'string', 'Cannot run, one or more names are empty');
    return;
%Check all names unique
elseif length(unique(fileNamesShort)) ~= length(fileNamesShort)
    set(handles.lastAction, 'string', 'Cannot run, names are not unique');
    return;
end

array_match = tabledata(:,handles.table.arrayMatch)';
%Check if array match cells are filled
if any(cellfun(@isempty,array_match))
    set(handles.lastAction, 'string', 'Cannot run, one or more array match cells are empty');
    return;
%Check if array match cells are numeric
elseif ~all(cellfun(@isnumeric,array_match))
    set(handles.lastAction, 'string', 'Cannot run, one or more array match cells are not numeric');
    return;
end

image_match = tabledata(:,handles.table.imageMatch)';
%Check if image_match cells are filled
if any(cellfun(@isempty,image_match))
    set(handles.lastAction, 'string', 'Cannot run, one or more image match cells are empty');
    return;
%check that image_match cells are numeric
elseif ~all(cellfun(@isnumeric,image_match))
    set(handles.lastAction, 'string', 'Cannot run, one or more array match cells are not numeric');
    return;
end

scale=str2num(handles.scale.String);
%check that the scale is numeric and not zero
if isnan(scale) || scale<=0
    set(handles.lastAction, 'string','Cannot run, invalid scale')
    return;
end

well_px=str2num(handles.well_px.String);
%check that the well_px is numeric and not zero
if isnan(well_px) || well_px<=0
    set(handles.lastAction,'string', 'Cannot run, invalid well spacing')
    return;
end

lane_px=str2num(handles.lane_px.String);
%check that lane_px is numeric and not zero
if isnan(lane_px) || lane_px<=0
    set(handles.lastAction, 'string','Cannot run, invalid lane spacing')
    return;
end

% Assign param.AO output 
global param
param.AO.base_folder = base_folder; % complete folder path (req: string)
param.AO.file = fileNamesFull; % filename list (req: string cell array 1xN, with file extension)
param.AO.name = fileNamesShort; % file nickname list (req: string cell array 1xN, with file extension)
param.AO.arrayTemplate = array_match; % array matching list for consistent ROI assignment for the same slide (req: positive integer cell array 1xN)
param.AO.transformTemplate = image_match; % image matching list for (req: positive integer cell array 1xN)
param.AO.scale_umperpx = scale; % um per pix scale (req: number, >0), resolution of image
param.AO.wellSpacing_px = well_px; % px (req: positive integer), well to well spacing
param.AO.laneSpacing_px = lane_px; % px (req: positive integer), lane to lane spacing

% Variables assigned elsewhwere 
% param.AO.minYcropOffset_um % um, minimum distance required on top and bottom of array (req: number)
% param.AO.multiDirectionROI % 0 = No, 1 = Yes: Create additional ROI up and down for all devices (req: 0 or 1)
% param.AO.refArrayTemplate % complete path for a reference array, or '' for blank (req: string)

%print param.AO struct
param.AO
    
%close GUI
close;

%------------------------------------------------------------------------
%--------------------------INSTRUCTIONS----------------------------------
%------------------------------------------------------------------------
% Executes when 'INSTRUCTIONS' button is pressed
% displays instructions in .txt file
function instructions_Callback(hObject, eventdata, handles)
% hObject    handle to instructions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open(fullfile('Instructions','UI_ArrayExtraction_instructions.txt'));

%------------------------------------------------------------------------
%-----------------FILE SELECTION ASSOCATED FUNCTIONS---------------------
%------------------------------------------------------------------------
% --- Executes when 'SELECT DIRECTORY' button pressed
% opens UI window to select directory
% verifies that tif files are found in the directory
% exports list of .tif files in folder to filelistbox
function folder_Callback(hObject, eventdata, handles)
% hObject    handle to folder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Select folder
folder = uigetdir();

% check selection
if isempty(folder) || isequal(folder,0)
    set(handles.lastAction,'string','no folder selected')
    return;
end

files = dir(fullfile(folder,'*.tif'));

% if no tif files found
if isempty(files)
    set(handles.lastAction,'string','no tif files found in the selected directory');
    errordlg('no .tif files found in the selected directory','Error');
    return;
end

% extract file list
set(handles.directory, 'String', folder);

% Display directory text in directory textbox
set(handles.filelistbox,'string',{files.name});

% Save the new values
guidata(hObject,handles)
set(handles.lastAction,'string','directory selected')

% on selection of line in filelistbox, files are moved to
% files2process_box, and removed from filelistbox
function filelistbox_Callback(hObject, eventdata, handles)
% hObject    handle to filelistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns filelistbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from filelistbox

%Extract the contents of the current listbox & fileNamesFull box
filelist = cellstr(get(hObject,'String'));
files2processFullNames = cellstr(get(handles.files2process_box,'string'));

%extract selected file in the present listbox
index = get(hObject,'value');
selected = filelist{index};

% check if default choice is made
if isequal(selected,'file list...')
    return;
end
   
% check if selected file is already in fileNamesFull, no action if present
% isempty(strfind(fileNamesFull, selected))

%Remove selected file form filelist
if index<= length(filelist)
    filelist(index)=[];
end
if isempty(filelist)
    filelist{1}='file list...';
end

%Prevents an odd bug where the listbox disapears if the 'value' is out of
%range of fileNamesFull, occurs after the last value is deleted
if index>1
    set(hObject,'value',index-1);
end

%append file to fileNamesFull
if isempty(files2processFullNames) || isequal(files2processFullNames{1},'files to be processed...')
    files2processFullNames={selected};
else
    files2processFullNames = [files2processFullNames', {selected}];
    files2processFullNames = sort(files2processFullNames);
end

%Set new fileNamesFull listbox
set(hObject,'string',filelist);
set(handles.files2process_box,'string',files2processFullNames);
set(handles.lastAction,'string','file added to be processed')

% on selection of line in files2process_box, files are moved to
% filelistbox, and removed from files2process_box
function files2process_box_Callback(hObject, eventdata, handles)
%Extract the contents of the current listbox & fileNamesFull box
files2processFullNames = cellstr(get(hObject,'String'));
filelist = cellstr(get(handles.filelistbox,'String'));

%extract selected file in the present listbox
index = get(hObject,'value');
selected = files2processFullNames{index};

% check if default choice is made
if isequal(selected,'files to be processed...')
    return;
end

% Remove selected file
if index<= length(files2processFullNames)
    files2processFullNames(index)=[];
end

if isempty(files2processFullNames)
    files2processFullNames{1}='files to be processed...';
end

%Prevents an odd bug where the listbox disapears if the 'value' is out of
%range of fileNamesFull, occurs after the last value is deleted
if index>1
    set(hObject,'value',index-1);
end

%append file to filelistbox
if isempty(filelist) || isequal(filelist{1},'file list...')
    filelist={selected};
else
    filelist = [filelist', {selected}];
    filelist = sort(filelist);
end

%Set new list
set(hObject,'string',files2processFullNames);
set(handles.filelistbox,'string',filelist);
set(handles.lastAction,'string','file removed from those to be processed')

% Executes when 'ADD ALL' button is pressed
% adds all files from filelist to files2process_box
function addAll_Callback(hObject, eventdata, handles)
filelist = get(handles.filelistbox,'string');
files2processFullNames = cellstr(get(handles.files2process_box,'string'));

  % check if default choice is made
if isa(filelist,'char') && isequal(filelist,'file list...')
    set(handles.lastAction,'string','no files to add...');
    return;
end

%append file to fileNamesFull
if isempty(files2processFullNames) || isequal(files2processFullNames{1},'files to be processed...')
    files2processFullNames=filelist;
else
    files2processFullNames = [files2processFullNames', filelist'];
    files2processFullNames = sort(files2processFullNames);
end

%reset filelist
filelist='file list...';

%save settings in figure object
set(handles.filelistbox,'string',filelist);
set(handles.filelistbox,'value',1);
set(handles.files2process_box,'string',files2processFullNames);
set(handles.lastAction,'string','all files were added to be processed');

%------------------------------------------------------------------------
%-------------- ASSIGNMENT TABLE ASSOCIATED FUNCTIONS--------------------
%------------------------------------------------------------------------

% Executes when 'ADD FILES TO TABLE' button is pressed
% adds all files from files2process_box  to the assignmentTable (clearing previous content)
% extracts tif image size, saves to table
% if on, executes arrayMatchON and/or imageMatchON
function tabletransfer_Callback(hObject, eventdata, handles)
%Load files
files2processFullNames = get(handles.files2process_box,'string');
if isequal(files2processFullNames,'files to be processed...')
    set(handles.lastAction,'string','no files to be processed');
    return;
end

%sort in alphabetical order
files2processFullNames = sort(files2processFullNames);

%insert file names into the data
data = [files2processFullNames, cell(length(files2processFullNames),4)];

%If array match is on, insert 1 into all cells
if get(handles.arrayMatchON,'Value')
    data(:,handles.table.arrayMatch)={1};
else
    data(:,handles.table.arrayMatch)={0};
end

%Loads image size for tif files into table
data = imageSizeSet(data,get(handles.directory,'string'), handles.table);

%If image match is on, match images that have the same final digits
if get(handles.imageMatchON,'Value')
    data = autoImageMatch(data, handles.table);
else
    data(:,handles.table.imageMatch)={0};
end

%set figure information
set(handles.assignmentTable,'Data',data);
set(handles.fileCount,'string',strcat('File count: ',num2str(length(files2processFullNames))));
set(handles.lastAction,'string','files to be processed added to to table')

% Called from tableTransfer
% extracts tif image info, saves image size into table
function data = imageSizeSet(data,directory,index)
    fileNamesFull = data(:,index.fileNamesFull);
    for i=1:(length(fileNamesFull))
        imageinfo=imfinfo(fullfile(directory,fileNamesFull{i}));
        width=imageinfo.Width;
        height=imageinfo.Height;

        %insert image size
        data(i,index.imageSize)={strcat(num2str(width),'x',num2str(height))};
    end

% Executes when 'ARRAY MATCH ON' button is pressed
% REQUIRED WHEN COMPARING IMAGES FROM THE SAME SLIDE, TO KEEP ROI CONSTANT
% adds '1' to all files (indicating that they have matching arrays) 
function arrayMatchON_Callback(hObject, eventdata, handles)
% hObject    handle to arrayMatchON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value')
    data = get(handles.assignmentTable,'Data');
    data(:,handles.table.arrayMatch) = {1}; 
    set(handles.assignmentTable,'Data',data);
end

% Executes when 'AUTO IMAGE MATCH ON' button is pressed
% RECOMMENDED WHEN A SLIDE IS SCANED FOR MULTIPLE WAVELENGTHS IN THE
% SAME MICROARRAY SCAN, (mitigates alignment error, reduces analysis steps)
% Adds a common number for image matched slides (stored in assignmentTable)
%   Performs image matching based on matching file names and image size.
%   File name matching is done by comparing the file before the last '_'.
%   e.g. examplefile_488.tif and examplefile_555.tif match if they
%   had the same image size
function imageMatchON_Callback(hObject, eventdata, handles)
% hObject    handle to imageMatchON (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%If image match is on, match images that have the same final digits
if get(handles.imageMatchON,'Value')
    data = get(handles.assignmentTable,'Data');
    data = autoImageMatch(data,handles.table); 
    set(handles.assignmentTable,'Data',data);
end
% Called from imageMatchON
function data = autoImageMatch(data,index)
   fileNamesFull = data(:,index.fileNamesFull);
   imsize = data(:,index.imageSize);
   %Splits files with '_' deliminted, removes final section
   for i=1:length(fileNamesFull)
        fileSplit{i} = strsplit(fileNamesFull{i},'_');
        if length(fileSplit{i}) == 2
            fileFirst(i)={(fileSplit{i}{1})};
        elseif length(fileSplit{i}) > 2
            fileFirst(i) = {strjoin({fileSplit{i}{1:(length(fileSplit{i})-1)}},'_')};
        else
            fileFirst(i) = {''};
        end
   end
   
   % Compares first section of strings in fileFirst, and imsize, to select
   % matchArray
   if ~isempty(fileFirst)
        imageMatchNum=0;
        imageMatch = zeros(length(fileFirst),1);
        for i=1:length(fileFirst)
            %intersects matches for name and size
            matchName = strmatch(fileFirst{i},fileFirst);
            matchSize = strmatch(imsize{i},imsize);
            matches = intersect(matchName,matchSize);
            
            %inserts matches into matchArray
            if length(matches)>1 && imageMatch(matches(1))==0
                imageMatchNum=imageMatchNum+1;
                imageMatch(matches)=imageMatchNum;
            end
        end
        %stores matches in the 5th column of data
        for i=1:length({data{:,1}})
            data{i,index.imageMatch}=imageMatch(i);
        end
   end

% --- Executes when entered data in editable cell(s) in assignmentTable.
% if the cell is editable, the edit is saved into the table
function assignmentTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to assignmentTable (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
set(handles.lastAction,'string','data added to table')

% Executes when 'CLEAR TABLE' is pressed, clears table...
function tableClear_Callback(hObject, eventdata, handles)
set(handles.assignmentTable,'Data',cell(6,4));
set(handles.lastAction,'string','table cleared');

%------------------------------------------------------------------------
% ----------------'REFERENCE ARRAY' ASSOCIATED FUNCTIONS-----------------
%------------------------------------------------------------------------

% --- Executes on button press in referanceArray.
function referanceArray_Callback(hObject, eventdata, handles)
% hObject    handle to referanceArray (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global param
[refArrayFile,refArrayFolder] = uigetfile('*.mat');
if isempty(refArrayFile) || isequal(refArrayFile,0)
    set(handles.lastAction,'string','no reference file selected');
    return;
end
refArrayFullFile=fullfile(refArrayFolder,refArrayFile);
param.AO.refArrayTemplate=refArrayFullFile;

% LOAD REFERENCE ARRAY DETAILS and DISPLAY:
load(refArrayFullFile,'cols','rows','name');

%Check if imageP file exists
imagePFile=fullfile(refArrayFolder,name);
imagePFile=strcat(imagePFile,'.tif');
if exist(imagePFile,'file')>0
    fileFound='yes';
else
    fileFound='Not Found';
end

%Output refArray details
refarraydetails = strcat('RefArray Details- Columns:', num2str(cols),' Rows:', num2str(rows),' Image found:', fileFound);
set(handles.refArrayDetails,'string',refarraydetails)

set(handles.refArray,'String',refArrayFile)
set(handles.lastAction,'string','reference array file added')

% --- Executes on button press in clearrefArray.
function clearrefArray_Callback(hObject, eventdata, handles)
% hObject    handle to clearrefArray (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.refArrayDetails,'string','')
set(handles.refArray,'String','no reference array selected')
set(handles.lastAction,'string','reference array file cleared')


%------------------------------------------------------------------------
% -------------ADVANCED OPTIONS POPUP-----------------
%------------------------------------------------------------------------


% --- Executes on button press in advancedSettingBTN.
function advancedSettingBTN_Callback(hObject, eventdata, handles)
% hObject    handle to advancedSettingBTN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global param
[param.AO.minYCropOffset_px, param.AO.multiDirectionROI,lastAcitonMSG] = advancedOptions(param.AO.minYCropOffset_px, param.AO.multiDirectionROI);
set(handles.lastAction,'string',lastAcitonMSG);

% Creates an Advanced Options Dialog box
% users set minycropoffset (um) & multiDirection (binary)
function [minYCropOffset, multiDirectionROI, lastActionMSG] = advancedOptions(arg1,arg2)
    global minYCropOffset multiDirectionROI lastActionMSG multiDirection_button minYCropOffset_edit minYCropOffset_default multiDirection_default
    % default settings loaded
    minYCropOffset_default=arg1;
    multiDirection_default=arg2;
    
    minYCropOffset = minYCropOffset_default;
    multiDirectionROI = multiDirection_default;
    lastActionMSG='advanced settings viewed';
    
    d = dialog('Position',[300 300 210 150],'Name','Advanced');
    
    minYCropOffset_edit = uicontrol('Parent',d,...
           'Style','edit',...
           'Position',[25 70 30 25],...
           'String', {num2str(minYCropOffset)}...
            );
    
    txt = uicontrol('Parent',d,...
           'Style','text',...
           'Position',[55 70 140 20],...
           'String','Minimum Y Offset Crop (px)');
       
        
    multiDirection_button = uicontrol('Parent',d,...
           'Style','radiobutton',...
           'Position',[25 100 200 25],...
           'String',{'Multidirection ROI generation'},...
           'Value', multiDirectionROI ...
            );
    
    %Save Button
    btn = uicontrol('Parent',d,...
           'Position',[20 20 60 25],...
           'String','Save',...
           'Callback',@callb_save);
    
    %Close Button
    btn2 = uicontrol('Parent',d,...
           'Position',[90 20 60 25],...
           'String','Close',...
           'Callback','delete(gcf)');
     
    % Wait for d to close before running to completion
    uiwait(d);
    
% Save inputs if dialog is saved
function callb_save(btn,callbackdata)
global minYCropOffset multiDirectionROI lastActionMSG multiDirection_button minYCropOffset_edit minYCropOffset_default multiDirection_default
  multiDirectionROI = multiDirection_button.Value;
  %Check that a valid value is entered
  if ~all(isstrprop(minYCropOffset_edit.String{1},'digit')) || isempty(minYCropOffset_edit.String{1})
      msgbox('minimum Y crop offset needs to be numeric');
      minYCropOffset = minYCropOffset_default;
      return;
  end
  minYCropOffset = str2num(minYCropOffset_edit.String{1});
  if ~isequal(minYCropOffset,minYCropOffset_default) || ~isequal(multiDirectionROI, multiDirection_default)
      lastActionMSG='advanced settings modified';
  end
  delete(gcf)
  return;

%------------------------------------------------------------------------
% -------------IMAGE/ARRAY SETTINGS ASSOCIATED FUNCTIONS-----------------
%------------------------------------------------------------------------

function scale_Callback(hObject, eventdata, handles)
scale = str2double(get(hObject, 'String'));
if isnan(scale)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new scale value
set(handles.lastAction,'string','scale set')
guidata(hObject, handles)

function well_px_Callback(hObject, eventdata, handles)
global param
well_px = str2num(get(hObject, 'String'));
if isnan(well_px) || well_px < 1
    set(hObject, 'String', param.AO.wellSpacing_px);
    errordlg('Input must be a positive integer','Error');
    return;
end
well_px=int16(well_px);
scale = str2num(get(handles.scale,'String'));

% Save the new well_px value
set(hObject, 'String', num2str(well_px));
set(handles.wellSpacing_um, 'String', num2str(well_px*scale));
set(handles.lastAction,'string','well to well spacing set')
guidata(hObject,handles)

function lane_px_Callback(hObject, eventdata, handles)
global param
lane_px = str2num(get(hObject, 'String'));
if isnan(lane_px) || lane_px < 1
    set(hObject, 'String', param.AO.laneSpacing_px);
    errordlg('Input must be a positive integer','Error');
    return;
end
lane_px=int16(lane_px);
scale = str2num(get(handles.scale,'String'));
% Save the new well_px value
set(hObject,'string',num2str(lane_px));
set(handles.laneSpacing_um,'string',num2str(lane_px*scale));
set(handles.lastAction,'string','lane distance set')
guidata(hObject,handles)

   
