%%% DATA EXTRACTION SETUP GRAPHIC USER INTERFACE %%%
%{ 
---------------------------------------------------------------------------
----------------------------GUI Insturctions-------------------------------
---------------------------------------------------------------------------
Instructions for the Data Extraction Setup GUI
Written by Todd A. Duncombe
    facilitated by MATLAB's GUIDE (gui construction interface)
Updated: 12/4/15

RECOMENDATIONS:
-

BEFORE USE:
-

(0) - 'SELECT DIRECTORY' (top left)
In streamlined operation, the directory will automatically be selected.
Otherwise an alternative directory can be selected. If not, select
directory that contains the cropped arrays - e.g. base_folder/prep/

(1) which files will have roi extracted, curve fit

...

(2) set curve fitting settings

'Subtraciton Mode' -
'Lanes to process' -
'Flank Subtraction (pixels)' -
'Spacer Height (pixels)' -
'Spacer Width (pixels)' -
'Loop through all selected files' -
'Export overlayed images' - 

OPTIONAL - ADVANCED USERS ONLY
(*) 'REFERENCE ARRAY'
Imports an array from a previously processed file for Array # '1' only.
Useful if adding data to a slide that was previously analyzed for other proteins.

---------------------------------------------------------------------------
-----------------DATA OUTPUT: 'param' (GLOBAL STRUCT VARIABLE)-------------
---------------------------------------------------------------------------
param.analysisStage
%}

%------------------------------------------------------------------------
%-------------------SYSTEM FUNCTIONS - DO NOT MODIFY---------------------
%---------------REQUIRED FOR 'GUIDE' - MATLAB GUI CONSTRUCTOR------------

function varargout = UI_Operational(varargin)
% UI_OPERATIONAL MATLAB code for UI_Operational.fig
%      UI_OPERATIONAL, by itself, creates a new UI_OPERATIONAL or raises the existing
%      singleton*.
%
%      H = UI_OPERATIONAL returns the handle to a new UI_OPERATIONAL or the handle to
%      the existing singleton*.
%
%      UI_OPERATIONAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in UI_OPERATIONAL.M with the given input arguments.
%
%      UI_OPERATIONAL('Property','Value',...) creates a new UI_OPERATIONAL or raises
%      the existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before UI_Operational_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to UI_Operational_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help UI_Operational

% Last Modified by GUIDE v2.5 09-Dec-2015 14:55:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @UI_Operational_OpeningFcn, ...
                   'gui_OutputFcn',  @UI_Operational_OutputFcn, ...
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

function varargout = UI_Operational_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

function initialize_gui(fig_handle, handles, isreset)
guidata(handles.figure1, handles);

function lanes2process_CreateFcn(hObject, eventdata, handles)
function subtract_CreateFcn(hObject, eventdata, handles)
function flankSubtract_CreateFcn(hObject, eventdata, handles)
function spacerHeight_CreateFcn(hObject, eventdata, handles)
function spacerWidth_CreateFcn(hObject, eventdata, handles)
function uitable2_CreateFcn(hObject, eventdata, handles)

%------------------------------------------------------------------------
%-------------------OPENING FUNCTION - (sets defaults)-------------------
%------------------------------------------------------------------------
function UI_Operational_OpeningFcn(hObject, eventdata, handles, varargin)
global param

handles.output = hObject;

% Table orientation
handles.table.name=1;
handles.table.refROI=2;
handles.table.extractROI=3;
handles.table.curveFit=4;
set(handles.uitable2,'Data',cell(4,length(fieldnames(handles.table))));

handles.externalReferenceMessage ='no external reference roi added';
set(handles.refROI_file,'string',handles.externalReferenceMessage);

% Set Defaults
subOptions ={'Average','Axial'};
set(handles.subtract,'string',subOptions);

if exist('param','var')
    if isfield(param,'AO') && isfield(param.AO,'filePrepPath') && ~isequal(param.AO.filePrepPath,'') 
        param.CF.base_folder=param.AO.filePrepPath;
        set(handles.directory,'string',param.CF.base_folder);
        addFiles(handles);
    end    
    if isfield(param,'CF')
        if isfield(param.CF,'lanes2process')
            set(handles.lanes2process,'string',num2str(param.CF.lanes2process));
        end
        if isfield(param.CF,'flankSubtract_px')
            set(handles.flankSubtract,'string',num2str(param.CF.flankSubtract_px));
        end
        if isfield(param.CF,'spacerHeight_px')
            set(handles.spacerHeight,'string',num2str(param.CF.spacerHeight_px));
        end
        if isfield(param.CF,'spacerWidth_px')
            set(handles.spacerWidth,'string',num2str(param.CF.spacerWidth_px));
        end
        if isfield(param.CF,'loop')
            set(handles.loop,'value',param.CF.loop);
        end
        if isfield(param.CF,'overlayArray')
            set(handles.overlayArray,'value',param.CF.overlayArray);
        end
        if isfield(param.CF,'refManSelect')
            set(handles.refManSelect,'value',param.CF.refManSelect);
        end
        if isfield(param.CF,'counter')
            set(handles.counter,'value',param.CF.counter);
        end
    end
else
    param = struct;
end

%Define analysis stage
param.analysisStage = 'Curve Fitting Setup GUI';

guidata(hObject, handles);
initialize_gui(hObject, handles, false);

%------------------------------------------------------------------------
%--------------------------CLOSING FUNCTION------------------------------
%------------------------------------------------------------------------

% --- Executes on button press in execute.
function execute_Callback(hObject, eventdata, handles)
% hObject    handle to execute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

base_folder = get(handles.directory,'string');
if isequal(base_folder,'selected directory...')
    set(handles.lastAction, 'string', 'Cannot run, no directroy selected')
    return;
end

tabledata = get(handles.uitable2,'data');
%Check table has been filled out
if isempty(tabledata)
    set(handles.lastAction, 'string', 'Cannot run, table not completed')
    return;
end

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

refROI = tabledata(:,handles.table.refROI)';
extractROI = tabledata(:,handles.table.extractROI)';
curveFit = tabledata(:,handles.table.curveFit)';
lanes2process=str2num(get(handles.lanes2process,'string'));
flankSubtract=str2num(get(handles.flankSubtract,'string'));
spacerHeight=str2num(get(handles.spacerHeight,'string'));
spacerWidth=str2num(get(handles.spacerWidth,'string'));
loop=get(handles.loop,'value');
overlayArray=get(handles.overlayArray,'value');
counter = get(handles.counter,'value');
refManSelect = get(handles.refManSelect,'value');

subI=get(handles.subtract,'value');
subS=get(handles.subtract,'string');
subtract=subS{subI}

% Assign param output 
global param
param.CF.base_folder = base_folder; % complete folder path (req: string)
param.CF.name = fileNamesShort; % file nickname list (req: string cell array 1xN, with file extension)
param.CF.refROI = refROI; % array matching list for consistent ROI assignment for the same slide (req: positive integer cell array 1xN)
param.CF.extractROI = extractROI; % image matching list for (req: positive integer cell array 1xN)
param.CF.curveFit = curveFit; % um per pix scale (req: number, >0), resolution of image
param.CF.subtraction = subtract; % px (req: positive integer), well to well spacing
param.CF.lanes2process = lanes2process; % px (req: positive integer), well to well spacing
param.CF.flankSubtract = flankSubtract; % px (req: positive integer), lane to lane spacing
param.CF.spacerHeight = spacerHeight; % image matching list for (req: positive integer cell array 1xN)
param.CF.spacerWidth = spacerWidth; % um per pix scale (req: number, >0), resolution of image
param.CF.loop = loop; % loop through all files, or just process first
param.CF.overlayArray = overlayArray; % px 
param.CF.counter= counter; % display counter during curve fit execution
param.CF.refManSelect= refManSelect; % display counter during curve fit execution

% Variables assigned elsewhwere 
%...

%print param struct
param.CF

close;
%mass = handles.metricdata.density * handles.metricdata.volume;
%set(handles.mass, 'String', mass);

%------------------------------------------------------------------------
%--------------------------INSTRUCTIONS----------------------------------
%------------------------------------------------------------------------

% --- Executes on button press in instructions.
function instructions_Callback(hObject, eventdata, handles)
% hObject    handle to instructions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%------------------------------------------------------------------------
%--------------CURVE FITTING FILE ASSIGNMENT IN TABLE--------------------
%------------------------------------------------------------------------
function addFiles(handles)
folder=get(handles.directory,'string');
if isempty(folder) || isequal(folder,0)
    set(handles.lastAction,'string','no folder selected')
    return;
end

files = dir(fullfile(folder,'*.mat'));

% if no .mat files found
if isempty(files)
    set(handles.lastAction,'string','no mat files were found in the selected directory');
    set(handles.directory ,'string','selected directory...');
    errordlg('no .mat files found in the selected directory','Error');
    return;
end

%removes .mat extensions
for i=1:length(files)
   [~,filenames{i},~]=fileparts(files(i).name); 
end

%load table
data=[filenames',cell(length(files),length(fieldnames(handles.table))-1)];

%insert automatic extractROI & curveFit
for i=1:size(data,1)
    data{i,handles.table.extractROI}='1';
    data{i,handles.table.curveFit}='1';
end

% Display directory text in directory textbox
set(handles.uitable2,'data',data);
% Save the new values
set(handles.lastAction,'string','directory selected')

% --- Executes when entered data in editable cell(s) in uitable2.
function uitable2_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable2 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

%load user input & table
index = eventdata.Indices;
enteredData = eventdata.EditData;
reference = get(handles.refROI_file,'string');
data = get(hObject,'data');

%Updating Reference ROI
%   must have only 1 reference ROI
if index(2)==handles.table.refROI
    for i=1:size(data,1)
        data{i,handles.table.refROI}='';
    end
    %if external reference roi is loaded, prevents table entry
    if strcmpi(reference,handles.externalReferenceMessage)
        data{index(1),handles.table.refROI}='1';
        data{index(1),handles.table.extractROI}='1';
        data{index(1),handles.table.curveFit}='1';
        set(handles.lastAction,'string','Reference ROI set');
    else
        errordlg('cannot add an internal & external reference ROI. Select "clear ref ROI" to modify.','Error');
        set(handles.lastAction,'string','Reference ROI conflict, clear external reference');
    end
%Updating CurveFit and ROI extract logic
% for CurFit ON, ROI extract must be ON
elseif index(2)==handles.table.extractROI
    if isequal(enteredData,'1')
        data{index(1),handles.table.extractROI}='1';
        set(handles.lastAction,'string','ROI extract ON');
    elseif isequal(enteredData,'0')
        data{index(1),handles.table.extractROI}='';
        data{index(1),handles.table.curveFit}='';
        set(handles.lastAction,'string','ROI extract OFF (CurveFit set OFF)');
    end
elseif index(2)==handles.table.curveFit
    if isequal(enteredData,'1')
        data{index(1),handles.table.curveFit}='1';
        data{index(1),handles.table.extractROI}='1';
        set(handles.lastAction,'string','CurveFit set ON (ROI extract ON)');
    elseif isequal(enteredData,'0')
        data{index(1),handles.table.curveFit}='';
        set(handles.lastAction,'string','CurveFit set OFF')
    end
end

%fill in remaining columns if they don't exist
if size(data,2) < 2
    data = [data,cell(length(files),length(fieldnames(handles.table))-1)];
    %automatically insert extractROI & curveFit
    for i=1:size(data,1)
        data{i,handles.table.extractROI}='1';
        data{i,handles.table.curveFit}='1';
    end
end

% Save table
set(handles.uitable2,'data',data);

%------------------------------------------------------------------------
% -------------------------SELECTING DIRECTORY/reset --------------------------
%------------------------------------------------------------------------

% --- Executes on button press in folder.
function folder_Callback(hObject, eventdata, handles)
global param
%Select folder, extract file list
folder = uigetdir();

% check selection
if isempty(folder) || isequal(folder,0)
    set(handles.lastAction,'string','no folder selected')
    return;
end

%load .mat iles in the folder
files=dir([folder,'\*.mat']);
% if no .mat files found
if size(files)==0
    set(handles.lastAction,'string','no mat files were found in the selected directory');
    errordlg('no .mat files found in the selected directory','Error');
    return;
end

%Set folder & addFiles to the table
param.CF.base_folder=folder;
set(handles.directory,'string',param.CF.base_folder);
addFiles(handles);

%------------------------------------------------------------------------
% ----------------'REFERENCE ROI' ASSOCIATED FUNCTIONS-------------------
%------------------------------------------------------------------------

% --- Executes on button press in refROI.
function refROI_Callback(hObject, eventdata, handles)
global param
%verifies that the selected file is a .mat
[refROIFile,refROIFolder] = uigetfile('*.mat');
%Checks that the file contains a '.mat'
if isempty(refROIFile) || isequal(refROIFile,0)
    set(handles.lastAction,'string','no reference file selected');
    return;
end

if ~strendswith(refROIFile,'_RoiData.mat')
    errordlg('The name of the external reference file must end with ''_RoiData.mat''','Error')
    set(handles.lastAction,'string','The name of the reference file must end with ''_RoiData.mat''');
    return;
end
refROIFullFile=fullfile(refROIFolder,refROIFile);
param.refROITemplate=refROIFullFile;

%Clear other ROI references
data = get(handles.uitable2,'data');
for i=1:size(data,1)
    data{i,2}='';
end
set(handles.uitable2,'data',data);

set(handles.refROI_file,'String',refROIFile)
set(handles.lastAction,'string','Reference ROI added')


%check whether a string ends with a given input
% from http://www.mathworks.com/matlabcentral/fileexchange/21710-string-toolkits/content/strings/strendswith.m
function b = strendswith(s, pat)
sl = length(s);
pl = length(pat);

b = (sl >= pl && strcmpi(s(sl-pl+1:sl), pat)) || isempty(pat);

% --- Executes on button press in clearRef.
function clearRef_Callback(hObject, eventdata, handles)
data = get(handles.uitable2,'data');
for i=1:size(data,1)
    data{i,2}='';
end
set(handles.uitable2,'data',data);
set(handles.refROI_file,'String',handles.externalReferenceMessage)
set(handles.lastAction,'string','Reference ROI cleared');

%------------------------------------------------------------------------
% -------------CURVE FITTING SETTINGS ASSOCIATED FUNCTIONS---------------
%------------------------------------------------------------------------

function lanes2process_Callback(hObject, eventdata, handles)
in = str2num(get(hObject, 'String'));
if isnan(in) || in<1
    set(hObject, 'String', '1');
    errordlg('Input must be a positive number','Error');
elseif in>2
    set(hObject, 'String', '1');
    errordlg('Input must be no longer than 2','Error');
end
in=int8(in);
set(hObject,'string',num2str(in));
set(handles.lastAction,'string','Lanes to process set');

function loop_Callback(hObject, eventdata, ~)
% Hint: get(hObject,'Value') returns toggle state of loop

function overlayArray_Callback(hObject, eventdata, handles)
% Hint: get(hObject,'Value') returns toggle state of overlayArray

function flankSubtract_Callback(hObject, eventdata, handles)
in = str2num(get(hObject, 'String'));
if isnan(in) || in<0
    set(hObject, 'String', '2');
    errordlg('Input must be a non negative number','Error');
end
in=int16(in);
set(hObject,'string',num2str(in));
set(handles.lastAction,'string','Flank subtraction width set');

% --- Executes on selection change in subtract.
function subtract_Callback(hObject, eventdata, handles)
set(handles.lastAction,'string','Subtraction mode set');

function spacerHeight_Callback(hObject, eventdata, handles)
in = str2num(get(hObject, 'String'));
if isnan(in) || in<0
    set(hObject, 'String', '0');
    errordlg('Input must be a non negative number','Error');
end
in=int16(in);
set(hObject,'string',num2str(in));
set(handles.lastAction,'string','Spacer height set');

function spacerWidth_Callback(hObject, eventdata, handles)
in = str2num(get(hObject, 'String'));
if isnan(in) || in < 0
    set(hObject, 'String', '0');
    errordlg('Input must be a non negative number','Error');u
end
in=int16(in);
set(hObject,'string',num2str(in));
set(handles.lastAction,'string','Spacer width set');


% --- Executes on button press in counter.
function counter_Callback(hObject, eventdata, handles)
    


% --- Executes on button press in refManSelect.
function refManSelect_Callback(hObject, eventdata, handles)
% hObject    handle to refManSelect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of refManSelect
