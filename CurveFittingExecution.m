%------------------------------------------------------------------------
%-------------------SYSTEM FUNCTIONS - DO NOT MODIFY---------------------
%---------------REQUIRED FOR GUIDE MATLAB GUI CONSTRUCTOR----------------

function varargout = CurveFittingExecution(varargin)
% CURVEFITTINGEXECUTION MATLAB code for CurveFittingExecution.fig
%      CURVEFITTINGEXECUTION, by itself, creates a new CURVEFITTINGEXECUTION or raises the existing
%      singleton*.
%
%      H = CURVEFITTINGEXECUTION returns the handle to a new CURVEFITTINGEXECUTION or the handle to
%      the existing singleton*.
%
%      CURVEFITTINGEXECUTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CURVEFITTINGEXECUTION.M with the given input arguments.
%
%      CURVEFITTINGEXECUTION('Property','Value',...) creates a new CURVEFITTINGEXECUTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CurveFittingExecution_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CurveFittingExecution_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CurveFittingExecution

% Last Modified by GUIDE v2.5 10-Dec-2015 10:58:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CurveFittingExecution_OpeningFcn, ...
                   'gui_OutputFcn',  @CurveFittingExecution_OutputFcn, ...
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
end
function initialize_gui(fig_handle, handles, isreset)
% If the metricdata field is present and the reset flag is false, it means
% we are we are just re-initializing a GUI by calling it from the cmd line
% while it is up. So, bail out as we dont want to reset the data.

% Update handles structure
guidata(handles.figure1, handles);
end

function varagout = CurveFittingExecution_OutputFcn(hObject, eventdata, handles) 
    varagout = 'Curve Fitting Execution Completed';
end

%------------------------------------------------------------------------
%-------------------OPENING FUNCTION - (sets defaults)-------------------
%------------------------------------------------------------------------

% --- Executes just before CurveFittingExecution is made visible.
function CurveFittingExecution_OpeningFcn(hObject, eventdata, handles, varargin)
    handles.output = hObject;
    guidata(hObject, handles);
    
    global resetAll slideStart param
    resetAll=1; %resets the status box (to ensure old content is cleared)
    
    %Set's analysis stage 
    param.analysisStage = 'Curve Fitting Execution';
    %Sets the first slide to start processing at
    %used for restarting on a slide
    slideStart=1;
    
    %Set default figure texts
    
    initialize_gui(hObject, handles, false);
    
    operational_curvefit(hObject, eventdata, handles,[]);
end

%------------------------------------------------------------------------
%----------------------MAIN ARRAY EXTRACITON FUNCTION--------------------
%------------------------------------------------------------------------

function operational_curvefit(hObject, eventdata, handles, globalGoodIn)
    global param slideStart fig h1 h3 globalGood
    
    %load globalGood if it has been set
    clear globalGood;
    if ~isempty(globalGoodIn)
        globalGood=globalGoodIn;
    end
    
    %set global graphics
    fig=gcf;
    set(fig,'pointer','arrow');
    set(fig,'units','normalized');
    set(fig,'outerposition',[0.05 0.1 0.9 0.85]);
    h1 = subplot(2,4,[2:4 6:8]);
    h3 = subplot(2,4,1);
    
    %Load instructions INCOMPLETE
    instructImage=struct();
    instructImage.example = imread(fullfile('Images','arrayInstructionExample.jpg'));
    imshow(instructImage.example,'parent',h3);
      
    %create folder & paths
    param.CF.fileProcPath = fullfile(param.CF.base_folder,'processing');
    if exist(param.CF.fileProcPath,'dir')~=7
        mkdir(param.CF.base_folder,'processing');
    end
    param.CF.fileOverlayPath=fullfile(param.CF.base_folder,'overlayArray');

    %Display folders in figure
    if length(param.CF.base_folder)>80
        short=param.CF.base_folder((length(param.CF.base_folder)-80):length(param.CF.base_folder));
        set(handles.inputFolder,'string',['Input: ...',short]);    
    else
        short=param.CF.base_folder(1:length(param.CF.base_folder));
        set(handles.inputFolder,'string',['Input:',short]);
    end
    
    if length(param.CF.base_folder)>50
        short=param.CF.fileProcPath((length(param.CF.fileProcPath)-80):length(param.CF.fileProcPath));
        set(handles.outputFolder,'string',['Output: ...',short]);
    else
        short=param.CF.fileProcPath(1:length(param.CF.fileProcPath));
        set(handles.outputFolder,'string',['Output:',short]);
    end
    
    % Load external reference array (if external rerence array exists
    if isfield(param.CF,'refROITemplate') && ~isequal(param.CF.refROITemplate,'')
        load(strcat(param.CF.refROITemplate,'_RoiData.mat'));
        globalGood=output{3};
    end
  
    %If an internal reference ROI exists, adjust file order to process it first,
    %saves file order in global variable 'param.CF.runOrder'
    fileOrder();
    
    % Assigns what file the curve fitting starts with, if it loops through
    % all files
    start=slideStart;
    if param.CF.loop
        last=length(param.CF.name);
    else
        last=slideStart;
    end
    
    %Main curve fitting loop
    for fileRun=start:last        
        %assign proper file index
        fileIndex=param.CF.runOrder(fileRun);
        file = param.CF.name{fileIndex};
        
        %Update statusBox
        set(handles.fileTitle,'str',['Current File:',file]);
        lineOut(strcat('-------------------file#',num2str(fileRun),'--------------------'),hObject, handles);
        statusBoxUpdate(fileIndex,hObject, handles);
        
        % Load Array information from in '/prep/' folder
        load(fullfile(param.CF.base_folder,[file,'.mat']));
        
        %insert data from loaded file
        cropOffsetX=cropOffsetX_adj;
        cropOffsetY=cropOffsetY_adj;
        devPerBlock = cols*rows;
        devPerRow = cols;
        
        % for the default array (with a single block of wells)
        % horzDevSpace and vertDevSpace are well and lane spacing
        horzDevSpace = param.AO.wellSpacing_px;
        vertDevSpace = param.AO.laneSpacing_px;
        
        %Recommend using a spacerWidth and Height equal to 0
        sepWidth=horzDevSpace-param.CF.spacerWidth_px;
        sepLength=vertDevSpace-param.CF.spacerHeight_px;
        %find upper left point in the array
        x0=cropOffsetX-round(sepWidth/2);
        y0=cropOffsetY+round(param.CF.spacerHeight_px/2);
        
        %define image filename, adding '.tif'
        name = strcat(file,'.tif'); 

        %For overrun analysis (include next lane)
        if param.CF.lanes2process==2
            sepLength = vertDevSpace*2-param.CF.spacerHeight_px;%over run correction
            if param.CF.lanes2process>(cropOffsetY)
                rows=rows-1;
            end
        end

        % Roi Analysis
            [roi, roiSub,roiSubAvg,roiData,roiDataSub,roiDataSubAvg,background] = ...
        roiTransformTD3(...
            param.CF.base_folder,name,1,devPerBlock,devPerRow,x0,y0,horzDevSpace,vertDevSpace,...
            sepLength,sepWidth,param.CF.spacerWidth_px,param.CF.spacerHeight_px,...
            0,0,param.CF.flankSubtraciton_px,param.CF.overlayArray);
        
        %name output images
        roiImgfilename=strcat('roi',name);
        roiSubImgfilename=strcat('roiSUB',name);
        roiSubAvgImgfilename=strcat('roiSUBAVG',name);
        
        %Saves ROI for all files with extractROI equal to '1'
        %   files with extractROI off only export an overlayed Array 
        %   (assuming that option is selected)
        if strcmp(param.CF.extractROI{fileIndex},'1');
            %save output images
            imwrite(roi,fullfile(param.CF.fileProcPath,roiImgfilename),'tiff');
            imwrite(roiSub,fullfile(param.CF.fileProcPath,roiSubImgfilename),'tiff');
            imwrite(roiSubAvg,fullfile(param.CF.fileProcPath,roiSubAvgImgfilename),'tiff');
        end
        
        %Set desired subtraction setting for fitting/analysis
        if strcmpi(param.CF.subtraction,'Average')
            roiSubSet=roiSubAvg;
            roiSubSetData=roiDataSubAvg;
        else
            roiSubSet=roiSub;
            roiSubSetData=roiDataSub;
        end        

        %Load reference from reference array
        %   or if globalGood does not exist, it creates ones
        if ~exist('globalGood','var')
            globalGood=ones(devPerBlock,1);
        end
        good=globalGood;
        output={};  
        output={roi,roiSubSet, good, devPerBlock,sepLength,sepWidth,param.CF.fileProcPath,file,roiData,roiSubSetData,background,param.CF.subtraction};
        
        %if curve fitting is on for the slide
        if strcmp(param.CF.curveFit{fileIndex},'1');
            %manual removal of poor separations
            good = manualBlockRemove(file,roi,roiSubSetData, devPerBlock,vertDevSpace,horzDevSpace,good, hObject,handles); 
          
            %Manual selection step for ref ROI
            if ~isempty(find(strcmp(param.CF.refROI,'1'))) && ...
                fileIndex==find(strcmp(param.CF.refROI,'1')) && ...
                param.CF.refManSelect
                
            % Insert param if the variables exist
                if ~isfield(param.CF,'columnPerView')
                    columnPerView=param.CF.columnPerView;
                else
                    columnPerView=30;
                end
                if ~isfield(param.CF,'rowPerView')
                    rowPerView=param.CF.rowPerView;
                else
                    rowPerView=8;
                end
                if ~isfield(param.CF,'defaultManualSelect')
                    default=param.CF.defaultManualSelect;
                else
                    default='off';%default 'on' or 'off'
                end
                
                globalGood=GoodManualSelect(param.CF.fileProcPath,file, param.CF.subtraction, sepWidth, columnPerView, rowPerView, default,hObject,handles);%'on' or 'off' to chose mode of selection 
                good=globalGood;
            end
            
            % Curve fitting algorithum
            % Fit a selected domain (only for axial fit)
            %   1 at the end turns off plotting, remove variable to include plotting
            results = scblotfitsckTDtranBound(roiSub,roiSubSetData, background,good,sepLength,sepWidth,roiData,1,hObject,handles);

            %output variables
            output={roi,roiSubSet, good, devPerBlock,sepLength,sepWidth,param.CF.fileProcPath,file,roiData,roiSubSetData,background,param.CF.subtraction};
            
            %Output Array
            lineOut(sprintf('%d/%d good/total ROIs',sum(good),length(good)),hObject, handles);
            
            %save output
            save(fullfile(param.CF.fileProcPath,[file,'_RoiData.mat']),'output', 'results');
        
        %save output for all files with extractROI (no curve fitting
        %   performed)
        elseif strcmp(param.CF.extractROI{fileIndex},'1');
            save(fullfile(param.CF.fileProcPath,[file,'_RoiData.mat']),'output');
        end
        %iterate startSlide
        slideStart=slideStart+1;
    end
    close;
    web('https://www.youtube.com/watch?v=dQw4w9WgXcQ');
end

%------------------------------------------------------------------------
%-----------------------------EXTRACT ROI--------------------------------
%------------------------------------------------------------------------

function [roi, roiSub, roiSubAvg, roiData,roiDataSub,roiDataSubAvg,background] = roiTransformTD3(...
    folder,image_file,numBlocks,devPerBlock,devPerRow,x0,y0,horzDevSpace,vertDevSpace,...
    sepLength,sepWidth,spacerWidth,spacerHeight,horzBlockSpace,vertBlockSpace,flankSubtract,overlayArray)
    global param
    
    % DEFINE THE FOLLOWING PARAMETERS
    % location of top left corner of ROI for top left device on slide
    %   x0, y0
    % number of blocks
    %   numBlocks = 1;
    % number of devices per block
    %   devPerBlock = cols*rows;
    % number of devices per row of block
    %   devPerRow = cols;
    % horizontal spacing between devices
    %   horzDevSpace = wellSpacing_px;
    % vertical spacing between devices
    %   vertDevSpace = laneSpacing_px;
    % length of separation lane
    %   sepLength = laneSpacing_px - spacerHeight;
    % width of separation lane
    %   sepWidth = wellSpacing_px - spacerWidth;
    % width (number of columns) of spacer between devices in row
    %   spacerWidth = 0;????
    % length (number of rows) of spacer between rows in the ROI
    %   spacerHeight = 0;

    
    % horizontal spacing between blocks (first well to first well)
    %   horzBlockSpace = 0;
    % vertical spacing between blocks (first well to first well)
    %   vertBlockSpace = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % load image
    image = imread(fullfile(folder,image_file));
    
    % exports overlayed array's into separate folder
    %   1 if yes 0 if no
    if param.CF.overlayArray 
        mkdir(folder,'overlay');
        roiOver = image(y0:(y0+devPerBlock/devPerRow*vertDevSpace),x0:(x0+devPerRow*horzDevSpace));
        imwrite(roiOver,fullfile(folder,'overlay',['roiOverlay',image_file]),'tiff');
    end

    % initialize output image
    roi = [];
    roiSub = [];
    roiSubAvg = [];

    % rearrange each block of wells into a row
    for i = 1:numBlocks

        % location of top left corner of first ROI in each block
        x1 = x0 + horzBlockSpace*(ceil(i/2)-1);
        y1 = y0 + vertBlockSpace*(mod(i-1, 2));

        % initialize placeholder matrix for concatenating devices into a row
        roi_block = [];
        roiSub_block = [];
        roiSubAvg_block = [];

        % loop through each device of the block
        for j = 1:devPerBlock

            % determine row and column number of device
            numRow = floor((j-1)/devPerRow)+1;
            numCol = mod(j-1, devPerRow)+1;

            % subset device from whole image
            xSpan = round(x1+(numCol-1)*horzDevSpace):round((x1+(numCol-1)*horzDevSpace)+(sepWidth-1));
            ySpan = (y1+(numRow-1)*vertDevSpace):(y1+(numRow-1)*vertDevSpace+(sepLength-1));

            %keyboard;
            if max(ySpan)>size(image,1)
             %   keyboard;
            end
            device = image(ySpan, xSpan);

            % sum rows for separation and spacer ROIs,subtract mean spacer intensity from separation intensity
            %and append to blockData
            deviceData = sum(device, 2);  
            deviceDataT = sum(device, 1);

            %Axial background subtraction for each pixel in ROI
                %find spacer ROIs to the left and right of the separation lane
                leftspacer=(xSpan(1):(xSpan(1) + (flankSubtract-1)));
                rightspacer=(xSpan(length(xSpan))+(1-flankSubtract)):xSpan(length(xSpan));
                deviceleftspacer=image(ySpan,leftspacer);
                devicerightspacer=image(ySpan,rightspacer);

                leftspacerint=sum(deviceleftspacer,2);
                rightspacerint=sum(devicerightspacer,2);
                avgAxialbackground=(leftspacerint+rightspacerint)/2;

                for k=1:size(device,2)
                    deviceSub(:,k)=device(:,k)-uint16((1/flankSubtract)*(smooth(avgAxialbackground,2)));
                end      

                %sum of background subtraction for axial & transverse
                deviceDatabsub = sum(deviceSub,2);
                deviceDatabsubT = sum(deviceSub,1);

            %Avg Background subtraction
                avgBackground = mean(avgAxialbackground);
                for k=1:size(device,2)
                    deviceSubAvg(:,k)=device(:,k)-uint16((1/flankSubtract)*(avgBackground));
                end

                %sum of background subtraction for axial & transverse
                deviceDatabsubAvg = deviceDatabsub-avgBackground;
                deviceDatabsubAvgT = deviceDatabsubT-avgBackground;


            % add to placeholder with spacer
            roi_block = [roi_block device];
            roiSub_block = [roiSub_block deviceSub];
            roiSubAvg_block = [roiSubAvg_block deviceSubAvg];

            roiDataSub{1,1}(:,j) = deviceDatabsub;
            roiDataSub{1,2}(:,j) = deviceDatabsubT;

            roiDataSubAvg{1,1}(:,j) = deviceDatabsubAvg;
            roiDataSubAvg{1,2}(:,j) = deviceDatabsubAvgT;

            roiData{1,1}(:,j)=deviceData;
            roiData{1,2}(:,j)=deviceDataT;
            background{1}(:,j)=avgAxialbackground;

        end

        % append row placeholder to final matrix with spacer
        roi = [roi; roi_block];
        roiSub = [roiSub; roiSub_block];
        roiSubAvg = [roiSubAvg; roiSubAvg_block];

    end
end

%------------------------------------------------------------------------
%---------------------------Curve fitting loop-----------------------------
%------------------------------------------------------------------------

function [results] = scblotfitsckTDtranBound(roiSub, roiDataSub, background, good,sepLength,sepWidth,rawroiData,...
    offPlot,hObject,handles)
    global h1 param h3 fig
    
    %scblotfits is a function which outputs cell arrays of gaussian fit
    %properties and calculated area under the curve measurements. results are
    %the output parameters based on an initial gaussian fit of the data, while
    %results2 uses iterative gaussian fitting and a region 4-5sigma away from the
    %peak center to determine the background (for improved background
    %subtraction and higher r-squared value gaussian fits). The AUC is determined
    %by summing the pixel intensities of the peak plus or minus 4 sigma from the peak center
    %as determined from the gaussian fit. It is recommended
    %you analyze the data in results2, not results.

    %Important Output:
    %results2 [cell array]: results 2 is a cell array containing 8 fields
    %results2{1}: matrix of the gaussian fit parameter, a (scaling parameter)
    %for each intensity profile
    %results2{2}: matrix of the gaussian fit parameter, b (peak center parameter)
    %for each intensity profile
    %results2{3}: matrix of the gaussian fit parameter, c (sigma parameter)
    %for each intensity profile
    %results2{4}: matrix of the integrated peak intensity (or AUC) over the range of -4
    %sigma to +4sigma from the peak center for each intensity profile
    %results2{5}: matrix of the calculated SNR values of each intensity profile
    %results2{6}: matrix of the R-squared values of the gaussian fit for each
    %intensity profile
    %results2{7}: matrix of the SSE of the gaussian fit for each intensity
    %profile
    %results2{8}: matrix of the mean background intensity value in the region
    %4-5sigma away from the peak center for each intensity profile

    %Initialize matrices/cell arrays for data storage
    xval=(5*(1:sepLength)');%scales to microns, as 1 px is 5 microns
    yval=(5*(1:sepWidth)');%scales to microns, as 1 px is 5 microns

    median=[];
    results={};
    results2={};
    goodfirstfit=zeros(size(good));
    mean4sigma=[];
    sigma4forfits=[];
    sigma4l=[];
    sigma5l=[];
    meansigma4=[];
    meansigma5=[];
    sample=[];
    maxpeak=[];
    blockfitData=[];
    blockfitData2=[];

    % set block index
    i = 1;
    
    % subset data Axial
    rawblock_data=rawroiData{i,1};
    block_data = roiDataSub{i,1};
    back = background{i,1};
    block_good = good; %selects all devices in the block
    ind_good = find(block_good);%finds the index of the 'good' fits
    %xvalues=repmat(xval,1,length(ind_good));%creates 
    block_data_good=block_data(:,ind_good);

    rawblock_data_good=rawblock_data(:,ind_good);
    background_good=back(:,ind_good);
     
    % plot "good" data
    delete(h1); h1=subplot(2,4,[2:4 6:8]);
    plot(h1,xval,block_data(:, ind_good), 'k-'); hold on;
    maxI=max(max(block_data_good));
    minI=min(max(block_data_good));
    axis([0 max(xval) minI maxI]);
   
    % initial "guesses" for curvefitting paramaters
    title(['Select Mean location ', int2str(i)]);
    [xg(1),yg(1)] = ginputCustom(1,'y','Select Mean Location');
    plot(h1,[xg(1) xg(1)],[minI maxI], 'color', 'g','linewidth',5);
    
    %Function that verifies that the selections are on the correct position
    %relative to one another
    function response = checkGuess(x1,x2,message)
        response=0;
        if x1<x2
            uiwait(msgbox(message));
        else
            response=1;
        end
    end
    
    title(['Select Right STD location ', int2str(i)]);
    done=0;
    while done==0
        [xg(2),yg(2)] = ginputCustom(1,'y','Select Right STD');
        done = checkGuess(xg(2),xg(1),'Right STD must be positioned to the right of the mean');
    end
    plot(h1,[xg(2) xg(2)],[minI maxI], 'color', 'b','linewidth',5);
    
    title(['Select Left STD location ', int2str(i)]);
    done=0;
    while done==0
        [xg(3),yg(3)] = ginputCustom(1,'y','Select Left STD');
        done = checkGuess(xg(1),xg(3),'Left STD must be positioned to the left of the mean');
    end
    plot(h1,[xg(3) xg(3)],[minI maxI], 'color', 'b','linewidth',5);
    
    title(['Select Right Bound location ', int2str(i)]);
    done=0;
    while done==0
        [xg(4),yg(4)] = ginputCustom(1,'y','Select Right Bound');
        done = checkGuess(xg(4),xg(2),'Right bound must be positioned to the right of the right STD');
    end
    plot(h1,[xg(4) xg(4)],[minI maxI],'color', 'r','linewidth',5);

    title(['Select Left Bound location ', int2str(i)]);
    done=0;
    while done==0
        [xg(5),yg(5)] = ginputCustom(1,'y','Select Left Bound');
        done = checkGuess(xg(3),xg(5),'Left bound must be positioned to the left of the left STD');
    end
    plot(h1,[xg(5) xg(5)],[minI maxI],'color', 'r','linewidth',5);    

    if xg(5)<=0
        xg(5)=xval(1);
    end
    if xg(4)>max(xval)
        xg(4)=max(xval);
    end
    
    a_guess = yg(1);% a is the y-scaling
    b_guess = xg(1); % b is the median
    c_guess = abs(xg(2)-xg(3)); % c is the sd
    bound_high = round(xg(4)/5);
    bound_low = round(xg(5)/5);
    bound = bound_low:bound_high;
    
    % define range for peak location
    minPeakLoc = 0;
    maxPeakLoc = max(xval);
    % fitting parameters
    
    ft = fittype('gauss1');
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    opts.Display = 'Off';
    opts.Lower = [-Inf;minPeakLoc;0];
    opts.Upper = [Inf;maxPeakLoc;Inf];
    opts.StartPoint = [a_guess;b_guess;c_guess];
    
    %TRANSVERSE FITTING
    % subset data Transverse with bounds selected

    for k=1:size(roiSub,2)/sepWidth
       x1=1+(k-1)*sepWidth;
       x2=(k)*sepWidth;
       block_dataT(:,k)=sum(roiSub(bound_low:bound_high,x1:x2),1);
    end
    block_data_goodT=block_dataT(:,ind_good);
   
    % plot "good" TRANSVERSE data
    delete(h1); h1=subplot(2,4,[2:4 6:8]);
    plot(h1,yval,block_dataT(:, ind_good), 'k-'); hold on;
    maxI=max(max(block_data_goodT));
    minI=min(min(block_data_goodT));
    axis([0 max(yval) minI maxI]);
    
    title(['Transverse Block ', int2str(i)]);
        
    % initial "guesses"
    title(['Select Mean location ', int2str(i)]);
    [xg(1),yg(1)] = ginputCustom(1,'y','Select Mean Location');
    plot(h1,[xg(1) xg(1)],[minI maxI], 'color','g','linewidth',5);
    
    title(['Select Right STD location ', int2str(i)]);
    done=0;
    while done==0
        [xg(2),yg(2)] = ginputCustom(1,'y','Select Right STD');
        done = checkGuess(xg(2),xg(1),'Right STD must be positioned to the right of the mean');
    end
    plot(h1,[xg(2) xg(2)],[minI maxI], 'color','b','linewidth',5);
    
    title(['Select Left STD location ', int2str(i)]);
    done=0;
    while done==0
        [xg(3),yg(3)] = ginputCustom(1,'y','Select Left STD');
        done = checkGuess(xg(1),xg(3),'Left STD must be positioned to the left of the mean');
    end
    plot(h1,[xg(3) xg(3)],[minI maxI], 'color','b','linewidth',5);
    
    %remove instructions
    delete(h3);
    
    a_guessT = yg(1);% a is the y-scaling
    b_guessT = xg(1); % b is the median
    c_guessT = abs(xg(2)-xg(3)); % c is the sd

    % define range for peak location
    minPeakLoc = 0;
    maxPeakLoc = max(yval);
    % fitting parameters
    ftT = fittype('gauss1');
    optsT = fitoptions('Method', 'NonlinearLeastSquares');
    optsT.Display = 'Off';
    optsT.Lower = [-Inf;minPeakLoc;0];
    optsT.Upper = [Inf;maxPeakLoc;Inf];
    optsT.StartPoint = [a_guessT;b_guessT;c_guessT];
    
    %if counter is on turn on display
    if param.CF.counter
        set(handles.counter,'visible','on');
        set(handles.counterTitle,'visible','on');
    end
    
    for j=1:length(ind_good);
    
        %Counter display
        statusBoxTimer(j,length(ind_good),hObject,handles);
        
    %AXIAL FIT    
    %background offset for fit
        [~,mi]=min(block_data_good(:,j));
        %offset_data(bound,j) =block_data_good(bound,j)-block_data_good(mi);
        
        xvalB=(5*(bound)');%scales to microns, as 1 px is 5 microns
        
        % perform initial fit
        [f, gof] = fit(xvalB, block_data_good(bound,j), ft, opts);
        fv=coeffvalues(f);

    % Background region defined as region 4-5 times the standard
    % deviation of the gaussian away from the peak center
        sigma4=(fv(3)/2)*4;
        sigma5=(fv(3)/2)*5;
        lbound=fv(2)-sigma4;
        rbound=fv(2)+sigma4;
        [l,li]=min(abs(xval-lbound));
        [r,ri]=min(abs(xval-rbound));
        xrange=xval(li:ri);
     
        fluorInt=sum(block_data_good(li:ri,j));
        fluorIntBound=sum(block_data_good(bound_low:bound_high,j));

        snr=(max(smooth(block_data_good(li:ri,j))))/std(background_good(:,j));
        rsq=gof.rsquare;
        sse=gof.sse;

    %Transverse FIT    
    %background offset for fit
        [~,mi]=min(block_data_goodT(:,j));
        offset_dataT(:,j) =block_data_goodT(:,j)-block_data_goodT(mi);
        
    % perform initial fit
        [fT, gofT] = fit(yval, offset_dataT(:,j), ftT, optsT);
        
        fvT=coeffvalues(fT);
    % Background region defined as region 4-5 times the standard
    % deviation of the gaussian away from the peak center
        sigma4T=(fvT(3)/2)*4;
        sigma5T=(fvT(3)/2)*5;
        lboundT=fvT(2)-sigma4T;
        rboundT=fvT(2)+sigma4T;
        [lT,liT]=min(abs(yval-lboundT));
        [rT,riT]=min(abs(yval-rboundT));
        yrange=yval(liT:riT);
     
        fluorIntT=sum(block_data_goodT(liT:riT,j));
 
        snrT=(max(smooth(block_data_goodT(liT:riT,j))))/std(background_good(:,j));
        rsqT=gofT.rsquare;
        sseT=gofT.sse;
                
        %Since some initial gaussian fits will not be very good due to poor
        %background subtraction, the best fits are used to determine the
        %median 4-5sigma range for the data set. 
        if exist('offPlot')
        elseif rsq>0.9 && rsqT>0.9
            
            goodfirstfit(j)=1;
            sigma4forfits=[sigma4forfits,sigma4];
            sigma4l=[sigma4l,(fv(2)+(fv(3)/2)*4)];
            sigma5l=[sigma5l,(fv(2)+(fv(3)/2)*5)];
            median=[median,fv(2)];
            maxpeak=[maxpeak,(max(smooth(block_data_good(li:ri,j))))];
            
            h2=subplot(2,4,[2:4]); hold on;
            plot(h2,xrange,block_data_good(li:ri,j));
            plot(h2,xvalB,offset_data(bound,j) ,'b')
            plot(h2,f,'r')
            title('Good Fits');
            legend('');
            
            h4 = subplot(2,4,[6:8]);
            plot(h4,yrange,block_data_goodT(liT:riT,j));
            plot(h4,yval,offset_dataT(:,j) ,'b')
            plot(h4,fT,'r')
            title('Trans');
            legend('');       
        end
        
        blockfitData(:,j) = [fv(1);fv(2);fv(3);fluorInt;snr;rsq;sse;bound_high;bound_low;fluorIntBound];
        blockfitDataT(:,j) = [fvT(1);fvT(2);fvT(3);fluorIntT;snrT;rsqT;sseT;];
        %plot(xvalues(:,j),offset_data(:,j) ,'b')
        %hold on
        %plot(f,'r')

    end

results{i,1}=blockfitData;
results{i,2}=blockfitDataT;

%turn off counter visibility
set(handles.counter,'visible','off');    
set(handles.counterTitle,'visible','off');    

end

%Status box timer
function statusBoxTimer (current, total ,hObject,handles)
    global param        
    if param.CF.counter
        newLine = ['Current/total ', num2str(current),'/',num2str(total)];
        set(handles.counter,'string',newLine);    
    end
end

%------------------------------------------------------------------------
%----------------MANUAL REMOVAL --------------
%------------------------------------------------------------------------

function good = manualBlockRemove(file, roi, roiData, devPerBlock,vertDevSpace,horzDevSpace,good, hObject, handles)

%scBlotRemove2 is an upgraded version of scBlotRemove, which allows users
%to manually inspect devices that are selected with a drag and drop
%selection box.

% good is an optional parameter
% good = ones(devPerBlock,1);

% loop through blocks

    % subset data
    block_data = roiData{1,1};
    block_roi = roi;
    block_good = good;
    
    % run BlockRemove on block
    block_good_new = BlockRemove(file,block_data, block_roi, block_good, horzDevSpace,hObject, handles);
    
    % append output back to good
    good = block_good_new;
end

function block_good = BlockRemove(file,blockData, blockROI, block_good, horzDevSpace, hObject, handles)
    global h1 btnValue
    
    % index for which lane ROIs are "good"
    ind_good = find(block_good);

    % set data in 'non-good' lane ROIs as 0
    blockData(:, setdiff((1:size(blockData,2)),ind_good)) = 0;

    % plot good lane ROIs data
    delete(h1); h1= subplot(2,4,[2:4 6:8]);
    hold on;
    plot(blockData(:, ind_good), 'k-');
    title(['Filename: ',file]); 
    xlabel('Location [Pixels]'); 
    
    %optional axis scalining, not applied here because seeing negative data
    %can be helpful in the manual removal 
    %maxI=max(max(blockData(:,ind_good)));
    %minI=min(min(blockData(:,ind_good)));
    %axis([0 length(blockData(:,1)) minI maxI]);
    
    % Initiate output
    devsForReview = [];
    
    %Ask if manual removal of data is desired
    done = uibuttons('Manually remove poor lanes?','yes','no',hObject, handles);
    done=done-1;
    %If manual removal of lanes was slected
    if done==0
        % Loop for seleciton
        while done == 0
            devsForReview = [devsForReview, locator(blockData)];
            done = uibuttons('Done with selection?','yes','no',hObject, handles);
        end    

        % reset plot
        hold off;

        % remove nonunique values
        devsForReview = unique(devsForReview);

        % modify block_good to remove thrown out devices
        if ~isempty(devsForReview)
            block_good(devsForReview) = 0;
        end

        % plot final data
        plot(blockData(:,find(block_good)), 'k-');
        title(['Filename: ',file]); 
        xlabel('Location [Pixels]'); 
        hold off;
        %keyboard;
        % ask if done, if not, recursively run function again
        allDone = uibuttons('Done with block"','yes','no',hObject, handles);
        if allDone == 0 
            block_good = BlockRemove(file,blockData, blockROI, block_good, horzDevSpace,hObject, handles);
        end
    end
    
    % function to select devices for review
    function devNums = locator(blockData)
        %operate external funciton selectdata to select
        [pointslist,xselect,yselect] = selectdata('selectionmode','rect');
        i=cellfun('isempty',xselect);
        ind=find(i==0);

        % initialize array of checked indices
        indsClicked = [];

        % loop through each clicked point
        for k = 1:size(ind,1)
            xclick=xselect{ind(k)};
            yclick=yselect{ind(k)};
            for v=1:length(xclick)
            % determine which device is closes to the click point
                [minValue, minInd] = min(abs(blockData(xclick(v),:)-yclick(v)));

            % add index to indsClicked
            indsClicked = [indsClicked minInd];
            end
        end

        % replot selected devices as red
        plot(blockData(:,indsClicked), 'r-');

        % output is the cliked devices
        devNums = indsClicked;
    end
end

%------------------------------------------------------------------------
%----------------------STATUS UPDATING FUNCTION--------------------------
%------------------------------------------------------------------------

%Displayes status updates in command line and in GUI statusBox
function lineOut(newLine,hObject, handles)
    global resetAll
    %output to command line
    disp(newLine);
    
    %saves memory of 'log' in this funciton only
    persistent log;
    if resetAll==1
        log=[];
        resetAll=0;
    end
    
    %appends new line
    log=[log,{newLine}];
    %updates statusBox log
    set(handles.statusBox,'string',log);
    set(handles.statusBox,'value',length(log));
    %returns the loglist
    logout=log;
end

%statusbox update at the begining of a new file
function statusBoxUpdate(fileIndex, hObject, handles)
    global param
    lineOut(['Current File: ',param.CF.name{fileIndex}],hObject, handles);
    if strcmp(param.CF.extractROI{fileIndex},'1');
        er ='Yes';
    else
        er= 'No';
    end
    if strcmp(param.CF.curveFit{fileIndex},'1');
        cf ='Yes';
    else
        cf= 'No';
    end
    if strcmp(param.CF.refROI{fileIndex},'1');
        rr ='Yes';
    else
        rr= 'No';
    end
    if ~isempty(find(strcmp(param.CF.refROI,'1')))
        re ='Yes';
    else
        re= 'No';
    end
 
    lineOut(sprintf('Extract ROI:%s   Curve Fit:%s', er,cf),hObject, handles);
    lineOut(sprintf('Reference ROI:%s   Reference Extracted:%s',rr,re),hObject, handles);
    lineOut(sprintf('Background Subtraciton:%s    Lanes to process:%s',param.CF.subtraction,num2str(param.CF.lanes2process)),hObject, handles);
end

%reorders the files, such that a refROI file is first
function fileOrder()
    global param
    param.CF.runOrder = 1:length(param.CF.name);    
    if ~isempty(find(strcmp(param.CF.refROI,'1')))
        refIndex = find(strcmp(param.CF.refROI,'1'));
        if refIndex>1 && refIndex < length(param.CF.runOrder)
            param.CF.runOrder = [refIndex, param.CF.runOrder(1:refIndex-1),...
                param.CF.runOrder(refIndex+1:length(param.CF.runOrder))];
        elseif refIndex==1
            param.CF.runOrder = [refIndex, param.CF.runOrder(refIndex+1:length(param.CF.runOrder))];
        elseif refIndex==length(param.CF.runOrder)
            param.CF.runOrder = [refIndex, param.CF.runOrder(1:refIndex-1)];
        end
    end  
end
%------------------------------------------------------------------------
%----------------ON IMAGE SELECTIONS -(assorted functions) --------------
%------------------------------------------------------------------------

function [good] = GoodManualSelect(folder,name,subtraction, sepWidth,columnOnScreen,rowOnScreen,mode,hObject,handles)
%Allows for manual selection of 'good' rois directly from the image
%click to the right of the image to move to the next section of the slide
    global param h1

    if strcmpi(subtraction,'axial')
       sub = 'roiSUB';
    elseif strcmpi(subtraction,'average')
       sub = 'roiSUBAVG';
    else
       sub = 'roi';
    end

    image = imread(fullfile(folder,[sub,name,'.tif']));

    index=(size(image,2)/sepWidth);
    if strcmpi(mode,'on')
        good=ones(1,index);
        mode='Keep Mode';
    else
        good=zeros(1,index);
        mode='Exclude Mode';
    end

    done2='n';
    section=1;
    while strcmpi(done2,'n')
        %turn off reset buttons
        set(handles.startSlide,'visible','off');
        set(handles.start,'visible','off');
        
        %turn on and set buttons
        set(handles.pbtn3,'visible','on');
        set(handles.pbtn4,'visible','on');
        set(handles.pbtn5,'visible','on');
        set(handles.section,'visible','on');
        set(handles.laneCount,'visible','on');
        
        set(handles.pbtn3,'string','Invert');
        set(handles.pbtn4,'string','Previous');
        set(handles.pbtn5,'string','Next');
        
        %sets the base index for a new row
        index0=(section-1)*columnOnScreen*rowOnScreen+1;
        x0=(index0-1)*sepWidth+1;
        if x0>=size(image,2)
           done2='y';
        else
            index02=index0+columnOnScreen*rowOnScreen-1;
            if(index02>size(good,2))
                index02=size(good,2);
            end        
            x02=index02*sepWidth;
            y0=1;
            y02=size(image,1);
            
            rows=floor((index02-index0)/columnOnScreen);
            if rows<((index02-index0)/columnOnScreen)
                rows=rows+1;
            end
            if rows>rowOnScreen
                rows=rowOnScreen;
            end
            %clears previous full image
            clear imagecmp;
            %creates new composite image image
            for i=1:rows
                %indexes for each column of imagecmp
                indexc1=index0+(i-1)*columnOnScreen;
                indexc2=index0+i*columnOnScreen-1;
                %checks for boundary overrun
                if indexc2 > size(good,2)
                    indexc2 = size(good,2);
                end
                %sets image values for the corresponding indexes
                xc1=(indexc1-1)*sepWidth+1;
                xc2=(indexc2)*sepWidth;
                
                %clears previous row image
                clear imagerow;
                %iterates throgh all columns to create row image
                for l=indexc1:indexc2
                   x1dev=1+(l-1)*sepWidth;
                   x2dev=x1dev+sepWidth-1;
                   
                   %extracts each roi
                   imagedev=image(y0:y02,x1dev:x2dev);
                   %performs a contrast optimization on each individual roi
                   imagedev=imadjust(imagedev,stretchlim(imagedev),[0.1 1]);
                   
                   %after the first step,sets base point of the image row
                   %equal to its existind domain
                   if exist('imagerow','var')
                       pos=size(imagerow,2);
                   else
                       pos=0;
                   end
                   %appends new column image to the row
                   imagerow(1:(y02-y0+1),(pos+1):(pos+sepWidth))=imagedev;
                end
                %calculates preexisting dimensions of imagecmp
                c1=1;
                c2=xc2-xc1+1;
                r1=1+size(image,1)*(i-1);
                r2=size(image,1)*(i);
                %appends imagerow on imagecmp
                imagecmp(r1:r2,c1:c2)=imagerow;
            end
            
            %verifies imagecmp exists
            %(it should always under normal operation)
            if exist('imagecmp','var')
                delete(h1); h1=subplot(2,4,[2:4 6:8]);
                imshow(imagecmp,'parent',h1); hold on;
                title(strcat(name,' - Section',num2str(section)));

                %Plot rectangles
                for i=1:rows
                    for k=1:(size(imagecmp,2)/sepWidth)
                        indexCurrent=index0+(i-1)*columnOnScreen+(k-1);
                        r1=1+size(image,1)*(i-1);
                        c1=1+sepWidth*(k-1);
                        r=rectangle('Position',[c1+2,r1-2,sepWidth-4,y02-4]);
                        if indexCurrent<=length(good)
                            if good(indexCurrent)==0
                               set(r,'edgecolor','r');
                            else
                               set(r,'edgecolor','g');
                            end
                        end
                    end  
                end
                
                
                %Position indexing of buttons off of the axes
                %assigns common uni 'pixels
                set(gcf,'units','pixels');
                set(handles.pbtn3,'units','pixels');
                set(handles.pbtn4,'units','pixels');
                set(handles.pbtn5,'units','pixels');
                
                %extracts pixel locaiton for each
                invertBtn=get(handles.pbtn3,'position');
                prevBtn=get(handles.pbtn4,'position');
                nextBtn=get(handles.pbtn5,'position');
                
                %user interface loop
                done='n';
                set(handles.section,'string',sprintf('Section %d',section));
                while done=='n'
                    set(handles.laneCount,'string',sprintf('%d/%d ROIs included',sum(good),length(good)));
                    
                    %ginputCustom is modified version of MATLAB's ginput, it enables
                    % pix = [x,y] pixel locations of the seleciton on the
                    % gcf
                    
                    [xg,yg,pix] = ginputCustom(1,'','Select ROIs');
                    
                    %if the seleciton is within field
                    %invert the 'good' value of the selected
                    if 0<xg && xg<size(imagecmp,2) && yg>0 && yg<size(imagecmp,1)
                        colSel = floor(xg/sepWidth)+1;  
                        rowSel = floor(yg/size(image,1))+1;
                        indexSel=index0+(rowSel-1)*columnOnScreen+(colSel-1);
                        r1=1+size(image,1)*(rowSel-1);
                        c1=1+sepWidth*(colSel-1);
                        r=rectangle('Position',[c1+2,r1-2,sepWidth-4,y02-4]);
                        if indexSel<=length(good)
                            if good(indexSel)==0
                                good(indexSel)=1;
                                set(r,'edgecolor','g');
                            else
                                good(indexSel)=0;
                                set(r,'edgecolor','r');
                            end
                        end
                        
                        
                    %invert
                    elseif pix(1)>invertBtn(1) && pix(1) < invertBtn(1)+invertBtn(3)...
                            && pix(2)>invertBtn(2) && pix(2) < invertBtn(2)+invertBtn(4)
                        for i=1:rows
                            for k=1:(size(imagecmp,2)/sepWidth)
                                indexCurrent=index0+(i-1)*columnOnScreen+(k-1);
                                r1=1+size(image,1)*(i-1);
                                c1=1+sepWidth*(k-1);
                                r=rectangle('Position',[c1+2,r1-2,sepWidth-4,y02-4]);
                                if indexCurrent<=length(good)
                                    good(indexCurrent)=abs(good(indexCurrent)-1);
                                    if good(indexCurrent)==0
                                       set(r,'edgecolor','r');
                                    else
                                       set(r,'edgecolor','g');
                                    end
                                end
                            end
                        end
                        
                    elseif pix(1)>prevBtn(1) && pix(1) < prevBtn(1)+prevBtn(3)...
                            && pix(2)>prevBtn(2) && pix(2) < prevBtn(2)+prevBtn(4)...
                            && section > 1
                        section = section-2;
                        done='y';
                    
                    elseif pix(1)>nextBtn(1) && pix(1) < nextBtn(1)+nextBtn(3)...
                            && pix(2)>nextBtn(2) && pix(2) < nextBtn(2)+nextBtn(4)   
                        done='y';
                    end
                    
                end
            end
            section=section+1;
        end
    end
    
    %reset figure units
    set(gcf,'units','normalized');
    set(handles.pbtn3,'units','normalized');
    set(handles.pbtn4,'units','normalized');
    set(handles.pbtn5,'units','normalized');
                
    %turn on reset buttons
    set(handles.startSlide,'visible','on');
    set(handles.start,'visible','on');
    
    %turn off reset buttons
    set(handles.pbtn3,'visible','off');
    set(handles.pbtn4,'visible','off');
    set(handles.pbtn5,'visible','off');
    set(handles.laneCount,'visible','off');
    set(handles.section,'visible','off');
end
     




%------------------------------------------------------------------------
%------------------CHECK FUNCITON??-----------------
%------------------------------------------------------------------------


%------------------------------------------------------------------------
%-------------------------UIBUTTONS GUI OPERATION------------------------
%------------------------------------------------------------------------

% Displays 2 buttons and a question uibuttons on figure,
% exports selection in 'select' (1 for btn1, 0 for btn2)
% reports response with lineOut
function select = uibuttons(question, btnTitle1, btnTitle2, hObject, handles)  
    global btnValue
    % displays question
    set(handles.question,'string',question);
    set(handles.question,'visible','on');
    
    set(handles.pbtn1,'string',btnTitle1);
    set(handles.pbtn1,'visible','on');
    set(handles.pbtn1,'value',0);
    
    set(handles.pbtn2,'string',btnTitle2);
    set(handles.pbtn2,'visible','on');
    set(handles.pbtn2,'value',0);

    guidata(hObject, handles);
    
    uiwait();
    %btnValue is set in callback funcitons pbtn1 or pbtn2
    select = btnValue;
    
    %hide buttons
    set(handles.question,'visible','off');
    set(handles.pbtn1,'visible','off');
    set(handles.pbtn2,'visible','off');
    
    if select==1
        lineOut(strcat(question,': ',btnTitle1),hObject, handles);
    else
        select=0;
        lineOut(strcat(question,': ',btnTitle2),hObject, handles);
    end
end

% --- Executes on button press in pbtn1.
function pbtn1_Callback(hObject, eventdata, handles)
    global btnValue
    btnValue=1;
    guidata(hObject,handles);
    uiresume();
end

% --- Executes on button press in pbtn2.
function pbtn2_Callback(hObject, eventdata, handles)
    global btnValue
    btnValue=0;
    guidata(hObject,handles);
    uiresume();
end

%------------------------------------------------------------------------
%GUI FUNCTIONS ASSOCIATED WITH RESTARTING ANALYSIS (FOR ALL OR FOR SLIDE)
%------------------------------------------------------------------------

% Exits array extraction 
% this will incur errors, and may result in inviable exported data
function quitFig_Callback(hObject, eventdata, handles)
    if strcmpi('Yes',questdlg('Do you want quit curve fitting? (will result in an incomplete output)','Quit Dialog'))
        close;
    end
end

% Restarts array execution from first slide
function start_Callback(hObject, eventdata, handles)
    global slideStart globalGood   
    if strcmpi('Yes',questdlg('Do you want restart curve fitting (from first slide)?','Restart Dialog'))
        %reset initial functions
        slideStart=1;
        if exist('globalGood','var')
            clear globalGood;
        end
        
        lineOut('-----------------------------------------------------------', hObject, handles);    
        lineOut('RESTART ARRAY EXECUTION FROM BEGINNING', hObject, handles);
        lineOut('-----------------------------------------------------------', hObject, handles);
        restartFunction(hObject,eventdata, handles)
        operational_curvefit(hObject, eventdata, handles,[]);
    end
end

% Restarts array extraction at current slide
function startSlide_Callback(hObject, eventdata, handles)
    global globalGood
    if strcmpi('Yes',questdlg('Do you want restart curve fitting (from current slide)?','Current Slide Restart Dialog'))
        lineOut('-------------RESTART SLIDE--------------', hObject, handles);
        restartFunction(hObject,eventdata, handles)
        operational_curvefit(hObject, eventdata, handles,globalGood);
    end
end

function restartFunction(hObject,eventdata, handles)
    % reset buttons
    %manual select associated UI controls
    set(handles.laneCount,'visible','off');
    set(handles.section,'visible','off');
    set(handles.pbtn3,'visible','off');
    set(handles.pbtn4,'visible','off');
    set(handles.pbtn5,'visible','off');
    
    %uibuttons
    set(handles.pbtn1,'visible','off');
    set(handles.pbtn2,'visible','off');
    set(handles.question,'visible','off');
    
    %curvefit counting
    set(handles.counter,'visible','off');
    set(handles.counterTitle,'visible','off');    

end

%------------------------------------------------------------------------
%---------------------OTHER GUI OBJECTS FROM UICONTROL-------------------
%------------------------------------------------------------------------

% --- Executes on slider movement.
function sld_Callback(hObject, eventdata, handles)
    global btnValue
    btnValue='sld';
    uiresume();
end

% --- Executes on button press in selectROI.
function selectROI_Callback(hObject, eventdata, handles)
    global btnValue
    btnValue='selectROI';
    uiresume();
end

% --- Executes on button press in accept.
function accept_Callback(hObject, eventdata, handles)
    global btnValue
    btnValue='accept';
    uiresume();
end

% --- Executes on slider movement.
function sldY_Callback(hObject, eventdata, handles)
    global btnValue
    btnValue='sldY';
    uiresume();
end

% --- Executes on slider movement.
function sldX_Callback(hObject, eventdata, handles)
    global btnValue
    btnValue='sldX';
    uiresume();
end

% --- Executes on button press in skip.
function skip_Callback(hObject, eventdata, handles)
    global btnValue
    btnValue='skip';
    uiresume();
end

% --- Executes on button press in pbtn3.
function pbtn3_Callback(hObject, eventdata, handles)
    global btnValue fig
    btnValue='pbtn3';
end    
    
% --- Executes on button press in pbtn4.
function pbtn4_Callback(hObject, eventdata, handles)
    global btnValue
    btnValue='pbtn4';
end

% --- Executes on button press in pbtn5.
function pbtn5_Callback(hObject, eventdata, handles)
    global btnValue fig
    btnValue='pbtn5';
end

%------------------------------------------------------------------------
%----------------------ASSORTED GUI SYSTEM FUNCTIONS---------------------
%------------------------------DO NOT MODIFY-----------------------------

% --- Executes during object creation, after setting all properties.
function sld_CreateFcn(hObject, eventdata, handles)
    % Hint: slider controls usually have a light gray background.
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
end

% --- Executes during object creation, after setting all properties.
function sldY_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
end

% --- Executes during object creation, after setting all properties.
function sldX_CreateFcn(hObject, eventdata, handles)
    if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor',[.9 .9 .9]);
    end
end

function statusBox_Callback(hObject, eventdata, handles)
end

function statusBox_CreateFcn(hObject, eventdata, handles)
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end
