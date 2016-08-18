%------------------------------------------------------------------------
%-------------------SYSTEM FUNCTIONS - DO NOT MODIFY---------------------
%---------------REQUIRED FOR GUIDE MATLAB GUI CONSTRUCTOR----------------

function varargout = ArrayExtracitonExecution(varargin)
% ARRAYEXTRACITONEXECUTION MATLAB code for ArrayExtracitonExecution.fig
%      ARRAYEXTRACITONEXECUTION, by itself, creates a new ARRAYEXTRACITONEXECUTION or raises the existing
%      singleton*.
%
%      H = ARRAYEXTRACITONEXECUTION returns the handle to a new ARRAYEXTRACITONEXECUTION or the handle to
%      the existing singleton*.
%
%      ARRAYEXTRACITONEXECUTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ARRAYEXTRACITONEXECUTION.M with the given input arguments.
%
%      ARRAYEXTRACITONEXECUTION('Property','Value',...) creates a new ARRAYEXTRACITONEXECUTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ArrayExtracitonExecution_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ArrayExtracitonExecution_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ArrayExtracitonExecution

% Last Modified by GUIDE v2.5 08-Dec-2015 18:58:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ArrayExtracitonExecution_OpeningFcn, ...
                   'gui_OutputFcn',  @ArrayExtracitonExecution_OutputFcn, ...
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

function varagout = ArrayExtracitonExecution_OutputFcn(hObject, eventdata, handles)
    varagout='Array Extraction Completed';
end

%------------------------------------------------------------------------
%-------------------OPENING FUNCTION - (sets defaults)-------------------
%------------------------------------------------------------------------

% --- Executes just before ArrayExtracitonExecution is made visible.
function ArrayExtracitonExecution_OpeningFcn(hObject, eventdata, handles, varargin)
    handles.output = hObject;
    guidata(hObject, handles);
    
    global resetAll slideStart
    resetAll=1;
    
    %Sets the first slide to start processing at
    %used for restarting on a slide
    slideStart=1;
    
    initialize_gui(hObject, handles, false);
    
    array_orientation_extraction_GUI(hObject, eventdata, handles);
end

%------------------------------------------------------------------------
%----------------------MAIN ARRAY EXTRACITON FUNCTION--------------------
%------------------------------------------------------------------------

function array_orientation_extraction_GUI(hObject, eventdata, handles)
    clearvars -except hObject eventdata handles param
    %Assign figure variables
    global param fig slideStart h1
    %define fig and set as fullscreen
    fig = gcf;
    set(fig,'units','normalized');
    set(fig,'outerposition',[0.05 0.1 0.9 0.85]);
    
    %turn off specific warnings
    warning('off','images:initSize:adjustingMag');
    warning('off','MATLAB:MKDIR:DirectoryExists');
    
    %Load instructions INCOMPLETE
    instructImage=struct();
    instructImage.example = imread(fullfile('Images','arrayInstructionExample.jpg'));
    
    %Default x and y offsets, they are adjusted if the image is too small
    param.analysisStage='Array Orientation Execution';
    cropOffsetX=param.AO.wellSpacing_px;
    cropOffsetY=param.AO.laneSpacing_px*2;

    %Create directories
    if exist(fullfile(param.AO.base_folder,'prep'),'dir')~=7
        mkdir(param.AO.base_folder,'prep');
    end
    param.AO.filePrepPath = fullfile(param.AO.base_folder,'prep');

    %Load information from previous array
    if ~strcmp(param.AO.refArrayTemplate,'')
        [colRefArray,rowRefArray,refArrayName,imageRefArray] = loadProcessedArray(param.AO.refArrayTemplate);       
    end
    
    h1 = subplot(2,4,[2:4 6:8]); % main display
    h3 = subplot(2, 4, 1); % instruction images
    imshow(instructImage.example,'parent',h3);
      
    % Array Orientation Execution
    done_full='n';
    
    %Assign the start slide
    %used for restarting at a slide
    f=slideStart;
        
    %Excecution Loop
    while done_full=='n'
        slideStart=f;
        %check if all files are processed
        if f>length(param.AO.file)
            done_full='y';
        else
            % lineOut displays input string in command line & on figure listbox
            lineOut(strcat('-------------------file#',num2str(f),'--------------------'),hObject, handles);
            lineOut('Current File:',hObject, handles);
            lineOut(param.AO.file{f},hObject, handles);
            lineOut(['Output name: ', param.AO.name{f},' .tif/.mat'],hObject, handles);
                        
            %assign file names
            filename=param.AO.file{f};
            filename_out=param.AO.name{f};
            %values are exported later on
            path=param.AO.base_folder;
            name=filename_out;

            %Open image file & optmizes contrast
            orig=imread(fullfile(param.AO.base_folder,filename));
            J=imadjust(orig,stretchlim(orig),[0.1 1]);

            %Duplicate Transformations or Arrays to images
            %reset setting for array and image templates matching
            duplicateArray=0;
            duplicateTransform=0;

            % check for reference array template match
            if ~strcmp(param.AO.refArrayTemplate,'') && param.AO.arrayTemplate{f}==1
                duplicateArray=1;
                colPrevious = colRefArray;
                rowPrevious = rowRefArray;
                arrayTemplateName=refArrayName;
                imageP=imageRefArray;
                % display current file name in command line
                lineOut('Reference Array Template:',hObject, handles);
                lineOut(refArrayName,hObject, handles);
            % check for previous array template match
            elseif f>1 && param.AO.arrayTemplate{f}~=0
                %identify matching arrayTemplates
                arrayMatch=find([param.AO.arrayTemplate{:}]== ...
                    param.AO.arrayTemplate{f});
                % check if match has already been processed
                if arrayMatch(1)<f
                    duplicateArray=1;
                    prevArrayIndex = arrayMatch(1);
                    filepath = fullfile(param.AO.filePrepPath, ...
                        [param.AO.name{prevArrayIndex},'.mat']);
                    %loads previous settings from array match
                    [colPrevious,rowPrevious,prevArrayName,imageP] = ...
                        loadProcessedArray(filepath);
                    arrayTemplateName=prevArrayName;
                    lineOut('Array Template:',hObject, handles);
                    lineOut(prevArrayName,hObject, handles);                
                end
            end

            %chech for Transform Template from image match
            if f>1 && param.AO.transformTemplate{f}==param.AO.transformTemplate{f-1} && param.AO.transformTemplate{f}~=0
                duplicateTransform=1;%duplicates the entire transform
                lineOut('Transform Template:',hObject, handles)
                lineOut(param.AO.file{f-1},hObject, handles);
            end
          
          %Array Transformation
          if duplicateTransform==0
                % Rotate selection
                A=wellSelect(J,'Array Rotation: select well at top of image');
                B=wellSelect(J,'Array Rotation: select well at bottom of image on a common row as top well',A);

                %Rotate array to match selection
                rot=atand((A.yWell-B.yWell)/(A.xWell-B.xWell));
                if rot<0
                    rot=180+rot;
                end 
                delete(h1); h1=subplot(2,4,[2:4 6:8]);
                imshow(J,'parent',h1);
                set(gcf,'pointer','arrow');
                if uibuttons('Rotate','Right','Left',hObject, handles)
                    rot=rot+180;   
                end

                % Manual Rotation function
                rot = manualRotation(J,rot,hObject, handles);
                
                %Perform final rotation
                origR = imrotate(orig,rot,'bilinear');
                J=imrotate(J,rot,'bilinear');

                % Array definition
                if ~exist('colPrevious','var') || (exist('colPrevious','var') && duplicateArray==0)
                    %Select well locations
                    
                    Left =wellSelect(J,'Array Definition: select well at Left-most column of array');
                    Up=wellSelect(J,'Array Definition: select well at Top row of the array');
                    Down =wellSelect(J,'Array Definition: select well at Bottom row of array');
                    Right=wellSelect(J,'Array Definition: select well at Right-most column of array');
                    RightxWell=Right.xWell;
                    DownyWell=Down.yWell;
                    cols = 1+round((abs(Left.xWell-Right.xWell))/param.AO.wellSpacing_px);
                    rows = 1+round((abs(Up.yWell-Down.yWell))/param.AO.laneSpacing_px);
                else
                    %Adopts previous array conditions
                    cols = colPrevious;
                    rows = rowPrevious;
                    if exist('imageP','var')
                        %If the previous image (imageP) exists, shows previous
                        %image for help with alignment
                        fig2=figure;
                        a=axes();
                        imshow(imageP,'parent',a);
                        hold on;
                        title('Array Template Image');
                        xlabel(arrayTemplateName);
                        hold off;
                        
                        %set gcf to main figure
                        figure(fig);
                        Left = wellSelect(J,'Array Definition: select well at Left-most column of array');
                        
                        %close template image
                        delete(a);
                        %if exist('fig2','var')
                        %    close(fig2);
                        %end
                        
                        %if selection is on the right side of the image (past half),
                        %transfers the location over accordinglys
                        if Left.xWell>size(origR,2)/2
                            Left.xWell=Left.xWell-(colPrevious-1)*param.AO.wellSpacing_px;
                        end
                    else% when imageP does not exist
                        % Select Upper left well (with EP 'down')
                        Left =wellSelect(J,'Array Definition: select well at Left-most column of array');
                    end
                    figure(fig);
                    %Upper well selection (with EP driection 'down')
                    Up=wellSelect(J,'Array Definition: select well at Top row of the array');

                    %Assigns right and down wells automatically
                    RightxWell=Left.xWell+(cols-1)*param.AO.wellSpacing_px;
                    DownyWell=Up.yWell+(rows-1)*param.AO.laneSpacing_px;

                end
                               
                %report array
                lineOut([num2str(cols),' columns, ',num2str(rows),' rows'],hObject, handles);
                
                delete(h1);
                h1 = subplot(2,4,[2:4 6:8]); % main display
                imshow(J,'parent',h1);
    
                % Manual shift function REMOVED FROM CODE
                % ROI find funciton
                [y0_U, y0_L, cropOffsetY_adj, x0_L, x0_R, cropOffsetX_adj] =  findROI(...
                    Up, DownyWell, cropOffsetY, Left, RightxWell, cropOffsetX,...
                    J, cols, rows);                
                zoom('on');
                
                %Verifies if the array fits in the image, and lets you go back
                go=dimCheck(cropOffsetX_adj,param.AO.wellSpacing_px,cropOffsetY_adj,param.AO.minYCropOffset_px,hObject, eventdata,handles);
                %go=0 if the array does not pass dim Check
                if go==0; 
                elseif uibuttons('Accept alignment','Yes','No (restart slide)',hObject, handles)
                
                    % Crop and save image
                    Jcrop=J(y0_U:y0_L,x0_L:x0_R);
                    origcrop=origR(y0_U:y0_L,x0_L:x0_R);
                    imageP=Jcrop;
                    imwrite(origcrop,fullfile(param.AO.filePrepPath,[filename_out,'.tif']),'tiff');
                    slideParam=param;
                    % Save Settings
                    save(fullfile(param.AO.filePrepPath,[filename_out,'.mat']),'slideParam','cropOffsetX_adj','cropOffsetY_adj','cols','rows','filename','name','path');
                    %Save array dimensions as 'previous' dimensions
                    colPrevious=cols;
                    rowPrevious=rows;             
                else
                    %If the array is not accepted in the preivous prompt, repeats
                    %closeh1 =subplot(2,4,[2:4 6:8]);,subplot(1,4,4)); %Close figures
                    lineOut('-------------RESTART SLIDE--------------', hObject, handles);
                    restartFunction(hObject,eventdata, handles);
                end
          else
              %Verifies if the array fits in the image, and lets you go back
              go=dimCheck(cropOffsetX_adj,param.AO.wellSpacing_px,cropOffsetY_adj,param.AO.minYCropOffset_px,hObject, eventdata,handles);
                if go==0; %go=0 if the array does not fit
                else % Duplicate treatment to different flourescent channel
                    origR = imrotate(orig,rot,'bilinear');
                    origcrop=origR(y0_U:y0_L,x0_L:x0_R);
                    imageP=Jcrop;
                    imwrite(origcrop,fullfile(param.AO.filePrepPath,[filename_out,'.tif']),'tiff');
                    slideParam=param;
                    save(fullfile(param.AO.filePrepPath,[filename_out,'.mat']),'slideParam','cropOffsetX_adj','cropOffsetY_adj','cols','rows','filename','name','path');
                end
          end
          f=f+1;  
        end
    end
    close;
end  

%------------------------------------------------------------------------
%--------------------------------FIND ROI--------------------------------
%------------------------------------------------------------------------

function [y0_U, y0_L, cropOffsetY_adj, x0_L, x0_R, cropOffsetX_adj] = ...
            findROI(Up, DownyWell, cropOffsetY, Left, RightxWell, cropOffsetX,...
                    imageIn, cols, rows)
        global param
        
        x0_L=Left.xWell-cropOffsetX;
        x0_R=RightxWell+cropOffsetX;
        cropOffsetX_adj=cropOffsetX;
        if x0_L<1 || x0_R>size(imageIn,2)
            if x0_L<(size(imageIn,2)-x0_R)
                cropOffsetX_adj=Left.xWell-1;
            else
                cropOffsetX_adj=(size(imageIn,2)-RightxWell);
            end
        x0_L=Left.xWell-cropOffsetX_adj;
        x0_R=RightxWell+cropOffsetX_adj;
        end

        y0_U=Up.yWell-cropOffsetY;
        y0_L=DownyWell+param.AO.laneSpacing_px+cropOffsetY;
        cropOffsetY_adj=cropOffsetY;
        if y0_U<1 || y0_L>size(imageIn,1)
            if y0_U<(size(imageIn,1)-y0_L)
                cropOffsetY_adj=Up.yWell-1;
            else
                cropOffsetY_adj=(size(imageIn,1)-DownyWell-param.AO.laneSpacing_px);
            end
        y0_U=Up.yWell-cropOffsetY_adj;
        y0_L=DownyWell+param.AO.laneSpacing_px+cropOffsetY_adj;
        end
        
        %plot array on current axes
        hold on;
        r=rectangle('Position',[x0_L,y0_U,x0_R-x0_L,y0_L-y0_U]);
        set(r,'edgecolor','blue');
        for i=1:cols
            x=Left.xWell+(i-1)*param.AO.wellSpacing_px;
            xbound=Left.xWell+(i-1)*param.AO.wellSpacing_px+param.AO.wellSpacing_px/2;
            plot([x; x],[y0_U; y0_L],'r');
        end
        for i=1:(rows+1)
            y1=Up.yWell+(i-1)*param.AO.laneSpacing_px;
            plot([x0_L; x0_R],[y1; y1],'black');
        end
end

%------------------------------------------------------------------------
%---------------------------MANUAL ROTATION -----------------------------
%------------------------------------------------------------------------

%opens a manual rotation user interface
function rot = manualRotation(J,rot, hObject, handles)
    global h1 h2
    %Rotates image
    JR = imrotate(J,rot,'bilinear');             
    %Displays rotated image
    manRotateDisp(JR);
    %User prompt to perform manual curve fitting
    if uibuttons('Manual Rotation','Yes','No',hObject, handles)
        zoom('off');
        rot=manualRotationSlide(J,rot, hObject, handles);
        lineOut(['  rotated: ',num2str(rot),' degrees'],hObject, handles);
    end
    
    %Displays image with horizontal lines
    function manRotateDisp(JR)
            delete(h1); h1 =  subplot(2,4,[2:4 6:8]);
            imshow(JR,'Parent',h1);
            hold on;
            title('Array Rotation Verification');
            xlabel('Verify that the array is horizontal');
            ylabel('<------ <------ Protein separations <------ <------');
            zoom('on');
            % Draw periodic horizontal line
            for i=1:floor(size(JR,1)/50);
                plot(h1,[1; size(JR,2)],[i*50; i*50],'r');
            end
    end

    function rot=manualRotationSlide(J,rot_in, hObject, handles)
        global  imr im dim btnValue
        rot=rot_in;
        im = J;
        btnValue='';
        
        delete(h1);
        h2=subplot(2,4,[6 7 8]);
        h1=subplot(2,4,[2:4 6:8]);
        imr = imrotate(im,rot_in,'bilinear');  
        imshow(imr,'parent',h1); hold on;
        
        %Select zoomed retion button       
        set(handles.selectROI,'string','Select Zoomed ROI');
        set(handles.selectROI,'backgroundcolor',[0 1 0]);
        set(handles.accept,'string','Cancel');
       
        set(handles.selectROI,'visible','on');
        set(handles.accept,'visible','on');
        
        %uiresume is called in ptbn1 or ptbnt2 selections
        uiwait();
        %Select ROI on image (or cancel out of funciton)
        if strcmp(btnValue,'selectROI')
            selectROIfunc;
        else strcmp(btnValue,'accept')        
            set(handles.selectROI,'visible','off');
            set(handles.accept,'visible','off');
            return;
        end
        
        %make buttons visible for rotation
        set(handles.selectROI,'backgroundcolor',[.94 .94 .94]);
        set(handles.sldTxt,'string','0 deg');
        set(handles.sld,'value',0.5);
        set(handles.accept,'string','Accept');
        set(handles.sldTxt,'visible','on');
        set(handles.sld,'visible','on');
        
        done=0;
        while done==0
            uiwait();
            if strcmp(btnValue,'selectROI')
                selectROIfunc;
            elseif strcmp(btnValue,'sld')
                sldRotate(handles);
            elseif strcmp(btnValue,'accept')
                done=1;
            end
        end
        
        set(handles.selectROI,'visible','off');
        set(handles.accept,'visible','off');
        set(handles.sldTxt,'visible','off');
        set(handles.sld,'visible','off');
            
        function selectROIfunc()
            % plot 'fullscreen' array image
            delete(h1); 
            h1= subplot(2,4,[2:4 6:8]);
            imshow(imr,'parent',h1); hold on;

            %select ROI zoom
            dim=rectSelect(imr,h1);

            %Plot zoomed selection
            [h1,h2,crit] = plotBoxes(imr,dim,h1,h2);
            plotLines(h2,crit);
        end

        function sldRotate(handles)
            %extract value from slider
            rotShift = get(handles.sld,'value')*10-5;
            %displays the degrees rotated on figure
            set(handles.sldTxt,'string',strcat(num2str(rotShift),' deg'));
            %rotates figure accordingly
            imr=imrotate(im,rot_in+rotShift);
            rot = rot_in+rotShift;
            [h1,h2,crit] = plotBoxes(imr,dim,h1,h2);
            plotLines(h2,crit);
        end
        
        function plotLines(hin,crit)
            axes(hin);
            title('Left Zoomed Box/Right Zoomed Box');
            % plot lines on boxes
            plot([crit crit],[0 crit],'blue');
            lines=20;
            for i=1:lines
                plot([0 crit*2],[i/lines*crit i/lines*crit],...
                    'color','r','linestyle','--','parent',h2);
            end
        end
    end
end

%plots left and right boxes from selected ROI
%returns the axis handles, this was requied for seamless operation with
%other nested functions (e.g. sldRotate)
function [h1,h2,crit] = plotBoxes(imr,dim,h1,h2)
    delete(h1); h1= subplot(2,4,[2 3 4]);
    imshow(imr,'parent',h1(1)); hold on;

    %find the smallest dimension of the box
    crit = min([dim(3),dim(4)]);
    %plot boxes on main image
    rectangle('parent',h1(1),'position',[dim(1) dim(2) crit crit],'edgecolor','r');
    rectangle('parent',h1(1),'position',[dim(1)+dim(3)-crit dim(2)+dim(4)-crit crit crit],'edgecolor','r');

    %Extract left and right images
    boxLx=(dim(1):(dim(1)+crit));
    boxLy=(dim(2):(dim(2)+crit));
    imboxL=imr(boxLy,boxLx);
    boxRx=((dim(1)+dim(3)-crit):(dim(1)+dim(3)));
    boxRy=((dim(2)+dim(4)-crit):(dim(2)+dim(4)));
    imboxR=imr(boxRy,boxRx);

    %contrast each sub image, then combine
    imboxL=imadjust(imboxL,stretchlim(imboxL),[0.1 1]);
    imboxR=imadjust(imboxR,stretchlim(imboxR),[0.1 1]);
    imbox=[imboxL, imboxR];

    %plot box sizes
    delete(h2);
    h2 = subplot(2,4,[6 7 8]);
    imshow(imbox,'parent',h2); hold on;
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
    %returns the loglist
    logout=log;
end

%------------------------------------------------------------------------
%----------------ON IMAGE SELECTIONS -(assorted functions) --------------
%------------------------------------------------------------------------
%function for selecting well locations on array
function [roiData] = wellSelect(image,stepTitle,previous)
    global fig h1
    
    % Creates roiData struct to store wellSelect outcome
    roiData=struct();
    roiData.title=stepTitle;
    roiData.x=1:size(image,2);
    roiData.y=1:size(image,1);
    
    % if previous roiData is passed, it will mark the the xWell location
    % with a red line (helps with ballpark alignment)
    if exist('previous','var')
        roiData.previousX=previous.xWell;
    end
    
    %%% ZOOM IN ON REGION OF INTEREST %%%
    %zoom1
    roiData.instructions = ... 
        'Zoom 1/2: Left click and drag cursor to zoom';
    roiData=zoomSelect(image,roiData);
    
    %zoom2
    roiData.instructions = ... 
        'Zoom 2/2: Left click and drag cursor to zoom';
    roiData=zoomSelect(roiData.image,roiData);
        
    %%% SELECT WELL CENTER &&&&
    % invert image to help with selection with crosshairs
    % (crosshairs are black and cannot be changed to my knowledge)
    
    %if invert image is desired
    %invertImage=max(max(roiData.image))-roiData.image;
    %contrastImage=imadjust(invertImage,stretchlim(invertImage),[0.1 1]);
    
    %maximize contrast
    contrastImage=imadjust(roiData.image,stretchlim(roiData.image),[0.1 1]);
    
    % plot & label zoomed image
    figure(fig);
    delete(h1); h1 =  subplot(2,4,[2:4 6:8]);
    imshow(contrastImage,'Parent',h1);
    set(gcf, 'units','normalized','outerposition',[0.05 0.1 0.9 0.85]);
    set(gcf, 'pointer', 'crosshair');
    title(roiData.title);
    xlabel('Select well center with a left click')
    
    %ginput to select well center
    [xT,yT]=ginputCustom(1);
    % ADD QUERY HERE??
    roiData.xWell=roiData.x(round(xT));
    roiData.yWell=roiData.y(round(yT));
end
%sub funciton of wellSelect for zooming in a well (2 zoom steps_
function [roiData] = zoomSelect(image,roiData)
    global h1 fig
    % optimizes image contrast
    J=imadjust(image,stretchlim(image),[0.1 1]);
    
    % shows image
    figure(fig);
    delete(h1); h1 =  subplot(2,4,[2:4 6:8]);
    imshow(J,'Parent',h1); hold on;
    
    % if previous roiData was passed to wellSelect, draws vertical red line
    % corresponding to previous.xWell location
    if isfield(roiData,'previousX')
        plot([roiData.previousX roiData.previousX], [0 size(image,1)],'r');
    end
    
    % instructions for user
    title(roiData.title);
    xlabel(roiData.instructions);
    
    %rect select function, verifys boundaries of the image
    rect = rectSelect(image,h1);    
    
    % crops image
    cordy=rect(2):(rect(2)+rect(4));
    cordx=rect(1):(rect(1)+rect(3));
    image=image(cordy,cordx);
    
    % shifts previousX such that it can be displayed again in new regiou
    if isfield(roiData,'previousX')
        roiData.previousX=roiData.previousX-rect(1);
    end
    
    % saves changes
    roiData.x=roiData.x(cordx);
    roiData.y=roiData.y(cordy);
    roiData.rect=rect;
    roiData.image=image;
end

%sub function of zoomSelect & others
%verifies that the rectangle selected is valid and exports error messages
%when it is outside of the image boundaries (or too small, 10px)
function [rect] = rectSelect(image,hin)
    valid=0;
    while valid==0
        % User input for zoom rectangle
        rect=getrectCustom(hin);
        rect=round(rect);
        
        %error when section is out of bounds on image
        if rect(1)<0 || rect(1)+rect(3)>size(image,2)...
                || rect(2)<0 || rect(2)+rect(4)>size(image,1)
            msg=errordlg('Invalid selection (out of bounds) - try again');
            uiwait(msg);
            
        %error when dimension of selection is less than 10 pixels in size
        elseif rect(3)<10 || rect(4)<10
            msg=errordlg('Invalid selection (ROI dimension < 10 pixels) - try again');
            uiwait(msg);
        else
            valid=1;
        end
        
    end
    
    % Code below modifies the selected region so that it is within the domain of image,
    % this section of code is no longer needed with the errordlg's above...
        % I left it here because it can be useful if one would like the
        % operation of the code to automatically set the rectangle to being
        % valid rather than producing error messages. (simply remove all
        % code above remove comment marks below:)
    %{
    rect=getrect();
    rect=round(rect);
    if rect(1)<0
        rect(3) = rect(3)+rect(1);
        rect(1)=1;
    end
    if rect(1)+rect(3)>size(image,2)
        rect(3) = size(image,2)-rect(1);
    end

    if rect(2)<0
        rect(4) = rect(4)+rect(2);
        rect(2)=1;
    end
    if rect(2)+rect(4)>size(image,1)
        rect(4) = size(image,1)-rect(2);
    end
    %}
end

%------------------------------------------------------------------------
%------------------ARRAY DIMENSION CHECK (PRIOR TO CROP)-----------------
%------------------------------------------------------------------------
%checs array dimensions prior to crop, if an error will proced it skips
%displays an erro
function [go] = dimCheck(xOffset,well2well,yOffset,minimumOffset,hObject, eventdata, handles)
    go=0;
    if xOffset<=well2well/2
        lineOut(strcat('CropOffset =',num2str(xOffset),' minimum=',num2str(well2well)),hObject, handles);
        statement='Your arrays x dimension is too large for this image size!';
        lineOut(statement,hObject, handles);
        lineOut('Please use a larger image or a smaller array',hObject, handles);
        dimArrayQuestion(statement,hObject, eventdata, handles);
    elseif yOffset<=minimumOffset/2
        lineOut(strcat('CropOffset =',num2str(yOffset),' minimum=',num2str(minimumOffset)),hObject, handles);
        statement='Your arrays x dimension is too large for this image size!';
        lineOut(statement ,hObject, handles);
        lineOut(statement,hObject, handles);
        lineOut('Please use a larger image or a smaller array',hObject, handles);
        dimArrayQuestion(statement,hObject, eventdata, handles);
    else
        go=1;    
    end
end

%asks user what they want to do if slide dimensions do not fit
function dimArrayQuestion(statement,hObject, eventdata ,handles)
    global startSlide    
    redo = 'Redo slide';
    restart = 'Restart all slides';
    con = 'Continue (slide will not be saved)';
    answer = questdlg(statement,'Out of bounds arrayd dialog',... %Statement, dialog title 
          redo,restart,con,... %buttons
          con); %Default
    if strcmp(answer,restart)
        startSlide=1;
        lineOut('-----------------------------------------------------------', hObject, handles);    
        lineOut('RESTART ARRAY EXECUTION FROM BEGINNING', hObject, handles);
        lineOut('-----------------------------------------------------------', hObject, handles);
        restartFunction(hObject,eventdata, handles)
    elseif strcmp(answer,redo)
        lineOut('-----------------------------------------------------------', hObject, handles);    
        lineOut('RESTART ARRAY EXECUTION FROM CURRENT SLIDE', hObject, handles);
        lineOut('-----------------------------------------------------------', hObject, handles);
        restartFunction(hObject,eventdata, handles)
    else
        lineOut('SLIDE NOT EXPORTED', hObject, handles);    
    end 
end

%------------------------------------------------------------------------
%----------------------LOAD PREVIOUS ARRAY DETAILS-----------------------
%------------------------------------------------------------------------

%loads previous array, load in the .mat filepath
function [colPrevious,rowPrevious,prevArrayName,imagePrevArray] = loadProcessedArray(filepath)
    load(filepath,'cols','rows','filename');
    prevArrayName=filename;
    colPrevious=cols;
    rowPrevious=rows;

    %check if refarray image can be found
    [refpath,refname,~]=fileparts(filepath);
    if exist(fullfile(refpath,strcat(refname,'.tif')),'file') == 2
        imagePrevArray=imread(fullfile(refpath,strcat(refname,'.tif')));
        imagePrevArray=imadjust(imagePrevArray,stretchlim(imagePrevArray),[0.1 1]);
    end
end

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
    if strcmp('Yes',questdlg('Do you want quit Array Extraction? (will result in an incomplete output)','Quit Dialog'))
        close;
    end
end

% Restarts array execution from first slide
function start_Callback(hObject, eventdata, handles)
    global slideStart    
    if strcmp('Yes',questdlg('Do you want restart Array Extraction (from first slide)?','Restart Dialog'))
        slideStart=1;
        lineOut('-----------------------------------------------------------', hObject, handles);    
        lineOut('RESTART ARRAY EXECUTION FROM BEGINNING', hObject, handles);
        lineOut('-----------------------------------------------------------', hObject, handles);
        restartFunction(hObject,eventdata, handles)
    end
end

% Restarts array extraction at current slide
function startSlide_Callback(hObject, eventdata, handles)
    if strcmp('Yes',questdlg('Do you want restart Array Extraction (from current slide)?','Current Slide Restart Dialog'))
        lineOut('-------------RESTART SLIDE--------------', hObject, handles);
        restartFunction(hObject,eventdata, handles)
    end
end

function restartFunction(hObject,eventdata, handles)
    % reset buttons
    set(handles.sld,'visible','off');
    set(handles.sldTxt,'visible','off');
    set(handles.selectROI,'visible','off');
    set(handles.accept,'visible','off');
    set(handles.pbtn1,'visible','off');
    set(handles.pbtn2,'visible','off');
    set(handles.question,'visible','off');
    guidata(hObject, handles);
    %Run main funciton
    array_orientation_extraction_GUI(hObject, eventdata, handles);
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
