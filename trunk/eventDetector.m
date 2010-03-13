function varargout = eventDetector(varargin)
% EVENTDETECTOR M-file for eventDetector.fig
%      EVENTDETECTOR, by itself, creates a new EVENTDETECTOR or raises the existing
%      singleton*.
%
%      H = EVENTDETECTOR returns the handle to a new EVENTDETECTOR or the handle to
%      the existing singleton*.
%
%      EVENTDETECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EVENTDETECTOR.M with the given input arguments.
%
%      EVENTDETECTOR('Property','Value',...) creates a new EVENTDETECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs
%      are
%      applied to the GUI before eventDetector_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to eventDetector_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help eventDetector

% Last Modified by GUIDE v2.5 28-Dec-2009 15:28:24

%% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @eventDetector_OpeningFcn, ...
    'gui_OutputFcn',  @eventDetector_OutputFcn, ...
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


%% --- Executes just before eventDetector is made visible.
function eventDetector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to eventDetector (see VARARGIN)

% Choose default command line output for eventDetector
handles.output = hObject;

%add figure toolbar to window
%set(hObject,'toolbar','figure');
set(handles.figMain,'KeyPressFcn',@keyboardShortcuts);
set(handles.lstSweeps,'KeyPressFcn',@keyboardShortcuts);

%init values
%set dimensions of plot windows
handles.plotWidthMin = 0;
handles.plotWidthMax = 1;
handles.plotVHeightMin = 0;
handles.plotVHeightMax = 1;

handles.plotCHeightMin = 0;
handles.plotCHeightMax = 1;

handles.plotEventWidthMin = 0;
handles.plotEventWidthMax = 1;
handles.plotEventHeightMin = 0;
handles.plotEventHeightMax = 1;

handles.plotEventWidthMinDefault  = 0.02;
handles.plotEventWidthMaxDefault  = 20000;
% handles.plotEventHeightMinDefault = 0;
% handles.plotEventHeightMaxDefault = 80;
handles.plotEventHeightMinDefault = 3;
handles.plotEventHeightMaxDefault = 40;

%set dimensions of plot windows
handles.plotWidthMinDefault = 0;
handles.plotWidthMaxDefault = 3;
handles.plotVHeightMinDefault = -250;
handles.plotVHeightMaxDefault = 250;

handles.plotCHeightMinDefault = -100;
handles.plotCHeightMaxDefault = 100;

handles.onFigure = 0;

%init file and pathname as null arrays
handles.pathname = [];
handles.filename = [];
handles.currentFileIndex = 0;

handles.currentSweep = 1;
handles.lineBaseline = 0;
handles.firstLoad = 1;
handles.firstMeanLoad = 1;
handles.firstDiffLoad = 1;
handles.detectedEvents_ms = [];
handles.detectedEvents = [];
handles.detectedEventsVoltages = [];
handles.axes5 = 0;
handles.analyzeEvents = 0;
handles.numOpenFiles = 0;
handles.numDetectedEvents = zeros(handles.numOpenFiles,1);
handles.numEvents = zeros(handles.numOpenFiles,1);
handles.baseline = 0;
handles.fishingBaseline = 0;
handles.defaultPathname = '';
handles.filename = '';
handles.buttonUp = [0, 0];
handles.buttonDown = [0, 0];
handles.numSignals = 0;
handles.sweepsAnalyzed = [];
handles.rmsNoise = 0;
handles.fishingFile = 0;
handles.rampingFile = 0;
handles.primerFile = 0;
handles.analyzeBounds = 0;
handles.analyzeUpperBound = 1;
handles.analyzeLowerBound = 0;
handles.limitAnalysisBetweenCursors = 0;
handles.minEventTime = 0.2;
handles.probingVoltage = -20;
handles.holdingVoltage = 50;
handles.ejectVoltage = -120;
handles.cutoffEventSum = 0;
handles.runningSumCutoffEvent = 0;
handles.keepCutoffEvents = 0;
handles.currentDNAnumber = 0;
handles.oldDNAnumber = 0;
handles.currentMoleculeColor = [0,0,1];
handles.rmsThreshold = 7;
handles.filtered = 0;
handles.cutoffEventFoundRamp = [0,0];
handles.cutoffEventVoltageTimingPrimer = [0,0,0];

%zero out open files listbox
set( handles.lstOpenFiles, 'String', {} );

%initialize cutoff event numbers
handles.cutoffEventSweep = 0;
handles.cutoffEventStartTime = 0;
handles.beginHighOut = 0;
handles.vChangeStart = 0;

%calculate progress bar location above main window
tempUnits = get(handles.figMain, 'Units');
set(handles.figMain, 'Units', 'pixels');
currentMainWindowPosition = get(handles.figMain,'Position');
screenSize = get(0, 'ScreenSize');
%normalize to be less than 1
handles.progressBarPos = [currentMainWindowPosition(1)/screenSize(3),(currentMainWindowPosition(2)+currentMainWindowPosition(4)+35)/screenSize(4)];
set(handles.figMain, 'Units', tempUnits);

%zero out lstSweeps
set( handles.lstSweeps, 'String', {} );
set( handles.lstSweepsExcluded, 'String', {} );
set( handles.edtThreshold, 'String', '7' );
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes eventDetector wait for user response (see UIRESUME)
% uiwait(handles.figMain);


%% --- Outputs from this function are returned to the command line.
function varargout = eventDetector_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%% --- Executes on selection change in lstSweeps.
function lstSweeps_Callback(hObject, eventdata, handles)
% hObject    handle to lstSweeps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns lstSweeps contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstSweeps

%get selected sweep number for lstSweeps
listSweep = get(handles.lstSweeps, 'String');
indexSweep = get(handles.lstSweeps, 'Value');
selectedSweep = listSweep{indexSweep};

selectedSweepInfo = sscanf(selectedSweep, 'Sweep %d.%d');
handles.currentSweep = selectedSweepInfo(2);
handles.currentFileIndex = selectedSweepInfo(1);


%set dimensions of plot windows
handles.plotWidthMin = handles.plotWidthMinDefault;
handles.plotWidthMax = handles.plotWidthMaxDefault;
handles.plotVHeightMin = handles.plotVHeightMinDefault;
handles.plotVHeightMax = handles.plotVHeightMaxDefault;

handles.plotCHeightMin = handles.plotCHeightMinDefault;
handles.plotCHeightMax = handles.plotCHeightMaxDefault;

%update current sweep and plot
loadCurrentSweep( hObject, eventdata, handles );

%get changes to handles made in loadCurrentSweep
handles = guidata( hObject );

% Update handles structure
guidata(hObject, handles);


%% --- Executes during object creation, after setting all properties.
function lstSweeps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstSweeps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Executes on selection change in lstSweepsExcluded.
function lstSweepsExcluded_Callback(hObject, eventdata, handles)
% hObject    handle to lstSweepsExcluded (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns lstSweepsExcluded contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstSweepsExcluded


%% --- Executes during object creation, after setting all properties.
function lstSweepsExcluded_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstSweepsExcluded (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Executes on button press in btnExclude.
function btnExclude_Callback(hObject, eventdata, handles)
% hObject    handle to btnExclude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listSweep = get(handles.lstSweeps, 'String');
listSweepExcluded = get(handles.lstSweepsExcluded, 'String');
indexSweep = get(handles.lstSweeps, 'Value');

if( indexSweep <= length( listSweep ) )
    selectedSweep = listSweep{indexSweep};

    listSweepExcluded{length(listSweepExcluded)+1} = selectedSweep;
    listSweep(indexSweep) = [];

    %rebuild and sort cell array
    sweepsNumbers = zeros( length(listSweepExcluded), 2 );
    for i = 1:length(listSweepExcluded)        
        sweepsNumbers(i,:) = sscanf(listSweepExcluded{i}, 'Sweep %d.%d');
    end

    %sort!
    sweepsNumbers = sortrows( sweepsNumbers );

    %rebuild cell array
    for i = 1:length( listSweepExcluded );
        listSweepExcluded{i} = ['Sweep ',num2str(sweepsNumbers(i,1)), '.', num2str(sweepsNumbers(i,2))];
    end

    %update list
    set(handles.lstSweeps, 'String', listSweep );
    set(handles.lstSweepsExcluded, 'String', listSweepExcluded);
    set(handles.lstSweepsExcluded, 'Enable', 'on');

    %if last element is selected, select new last
    if( indexSweep >= length( listSweep ) && indexSweep > 1 )
        indexSweep = indexSweep - 1;
        set(handles.lstSweeps, 'Value', indexSweep );
    end

    handles.numEvents(handles.currentFileIndex) = handles.numEvents(handles.currentFileIndex) - 1; %length(listSweep);
    %disp( handles.numEvents );
    
    if( isequal( handles.sweepsAnalyzed{ handles.currentFileIndex }(handles.currentSweep) , 1 ) )
        sweepEvents = find( ((handles.detectedEvents(:,5) == handles.currentSweep) & (handles.detectedEvents(:,10) == handles.currentFileIndex)) );
        %load pre-calculated event info into timingInfo
        timingInfo = handles.detectedEvents( sweepEvents, 2:3 );
        [m,n] = size(timingInfo);
        handles.numDetectedEvents(handles.currentFileIndex) = handles.numDetectedEvents(handles.currentFileIndex) - m;
        set( handles.txtNumEvents, 'String', ['Detected Events: ', num2str(sum(handles.numDetectedEvents))] );
    end
        
    lstSweeps_Callback(hObject, eventdata, handles)
    handles = guidata( hObject );
end

% Update handles structure
guidata(hObject, handles);


%% --- Executes on button press in btnInclude.
function btnInclude_Callback(hObject, eventdata, handles)
% hObject    handle to btnInclude (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

listSweep = get(handles.lstSweeps, 'String');
listSweepExcluded = get(handles.lstSweepsExcluded, 'String');
indexSweep = get(handles.lstSweepsExcluded, 'Value');

%make sure list is non-empty
if( indexSweep <= length( listSweepExcluded ) )
    selectedSweep = listSweepExcluded{indexSweep};

    listSweep{length(listSweep)+1} = selectedSweep;
    listSweepExcluded(indexSweep) = [];

    %rebuild and sort cell array
    sweepsNumbers = zeros( length(listSweep), 2 );
    for i = 1:length(listSweep)
        sweepsNumbers(i,:) = sscanf( listSweep{i}, 'Sweep %d.%d');
    end

    %sort!
    sweepsNumbers = sortrows( sweepsNumbers );

    %rebuild cell array
    for i = 1:length( listSweep );
        listSweep{i} = ['Sweep ', num2str(sweepsNumbers(i,1)), '.', num2str(sweepsNumbers(i,2))];
    end

    set(handles.lstSweeps, 'String', listSweep);
    set(handles.lstSweepsExcluded, 'String', sort(listSweepExcluded) );
    set(handles.lstSweepsExcluded, 'Enable', 'on');

    %if last element is selected, select new last
    if( indexSweep >= length( listSweepExcluded ) && indexSweep > 1 )
        indexSweep = indexSweep - 1;
        set(handles.lstSweepsExcluded, 'Value', indexSweep );
    end

    handles.numEvents(handles.currentFileIndex) = handles.numEvents(handles.currentFileIndex) + 1; %length(listSweep);
    selectedSweepInfo = sscanf(selectedSweep, 'Sweep %d.%d');
    selectedSweepNum = selectedSweepInfo(2);
    selectedSweepFile = selectedSweepInfo(1);
    
    if( isequal( handles.sweepsAnalyzed{handles.currentFileIndex}( selectedSweepNum ) , 1 ) )
        sweepEvents = find( ((handles.detectedEvents(:,5) == selectedSweepNum) & (handles.detectedEvents(:,10) == selectedSweepFile)) );
        %load pre-calculated event info into timingInfo
        timingInfo = handles.detectedEvents( sweepEvents, 2:3 );
        [m,n] = size(timingInfo);
        handles.numDetectedEvents(handles.currentFileIndex) = handles.numDetectedEvents(handles.currentFileIndex) + m;
        set( handles.txtNumEvents, 'String', ['Detected Events: ', num2str(sum(handles.numDetectedEvents))] );
    end
end

% Update handles structure
guidata(hObject, handles);


%% --- Executes on changed of edtABFDescription
function edtABFDescription_Callback(hObject, eventdata, handles)
% hObject    handle to edtABFDescription (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtABFDescription as text
%        str2double(get(hObject,'String')) returns contents of edtABFDescription as a double


%% --- Executes during object creation, after setting all properties.
function edtABFDescription_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtABFDescription (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Executes on button press in btnExportEvent.
function btnExportEvent_Callback(hObject, eventdata, handles)
% hObject    handle to btnExportEvent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%disableButtons(handles);

%update baseline from edit box
handles.axes5 = figure(1);
handles.baseline = str2double( get( handles.edtBaseline, 'String') );

%init figure
set(handles.axes5, 'Position', [200 100 1000 450], 'Name', 'Detected Events Export');

%poreEventData = handles;
poreEventData.detectedEvents = handles.detectedEvents;
poreEventData.detectedEvents_ms = handles.detectedEvents_ms;
% poreEventData.listSweeps = get( handles.lstSweeps, 'String');
% poreEventData.listSweepsExcluded = get( handles.lstSweepsExcluded, 'String');


%build excluded sweep numbers
excludedSweeps = get(handles.lstSweepsExcluded, 'String');
excludedSweepsNumbers = zeros( length( excludedSweeps ), 1 );
if(  ~isempty(excludedSweeps) )


    for i = 1:length( excludedSweeps )
        excludedSweepInfo = sscanf(excludedSweeps{i}, 'Sweep %d.%d');
        excludedSweepsNumbers(i) = excludedSweepInfo(2);
        excludedFileNumbers(i) = excludedSweepInfo(1);
        poreEventData.sweepsAnalyzed(excludedSweepsNumbers(i)) = 0;
    end

    indsToDel = [];

    for sweepNum = 1:size(excludedSweepsNumbers,1)

        [inds, vals] = find( (handles.detectedEvents(:,1) == excludedSweepsNumbers(sweepNum)) & (handles.detectedEvents(:,10) == excludedFileNumbers(sweepNum)));
        if( ~isempty( inds ) )
            indsToDel = [indsToDel; inds];
        end
        
    end

    %trim off un-used space due to pre-allocating
    if( ~isempty(indsToDel) )
        poreEventData.detectedEvents(indsToDel,:) = [];
        poreEventData.detectedEvents_ms(indsToDel,:) = [];
    end

end

%plot events
[m,n] = size( poreEventData.detectedEvents_ms );
if( m >= 1) %if any events have been detected
    if( n >= 9 ) %dna molecule information included
        numberDNAmolecules = max( poreEventData.detectedEvents_ms(:,9) );
        for i = 0:numberDNAmolecules
            if( i == 0 )
                moleculeColor = [0,0,1];
            else
                moleculeColor = [rand, rand, rand];
            end
            
            try
		if( handles.fishingFile == 1 )
 		    nonCutoffEvents = find( (poreEventData.detectedEvents_ms(:,8) == 0) & (poreEventData.detectedEvents_ms(:,9) == i) & (poreEventData.detectedEvents_ms(:,10) == 3));
%             nonCutoffEvents = find( (poreEventData.detectedEvents_ms(:,8) == 0) & (poreEventData.detectedEvents_ms(:,9) == i) );
 		    nonCutoffEventsInitial = find( (poreEventData.detectedEvents_ms(:,8) == 0) & (poreEventData.detectedEvents_ms(:,9) == i) & (poreEventData.detectedEvents_ms(:,10) ~= 3));

		else
% 		    nonCutoffEvents = find( (poreEventData.detectedEvents_ms(:,8) == 0) & (poreEventData.detectedEvents_ms(:,9) == i) & (poreEventData.detectedEvents_ms(:,10) ~= 2));
% 		    nonCutoffEventsInitial = find( (poreEventData.detectedEvents_ms(:,8) == 0) & (poreEventData.detectedEvents_ms(:,9) == i) & (poreEventData.detectedEvents_ms(:,10) == 2));
             
%             nonCutoffEvents = find( (poreEventData.detectedEvents_ms(:,8) == 0) & (poreEventData.detectedEvents_ms(:,9) == i));
            nonCutoffEvents = find( (poreEventData.detectedEvents_ms(:,8) == 0) & (poreEventData.detectedEvents_ms(:,9) == i) );
		    nonCutoffEventsInitial = [];
        end		
                semilogx( poreEventData.detectedEvents_ms(nonCutoffEventsInitial ,4), poreEventData.detectedEvents_ms(nonCutoffEventsInitial,5), 's', 'Color', moleculeColor );
                semilogx( poreEventData.detectedEvents_ms(nonCutoffEvents ,4), poreEventData.detectedEvents_ms(nonCutoffEvents,5), '.', 'Color', moleculeColor );

            	%nonCutoffEventsInitial = find( (poreEventData.detectedEvents_ms(:,8) == 0) & (poreEventData.detectedEvents_ms(:,9) == i) & (poreEventData.detectedEvents_ms(:,10) == 2));
%             	semilogx( poreEventData.detectedEvents_ms(nonCutoffEvents ,4), poreEventData.detectedEvents_ms(nonCutoffEvents,5), '.', 'Color', moleculeColor );
            	%semilogx( poreEventData.detectedEvents_ms(nonCutoffEventsInitial ,4), poreEventData.detectedEvents_ms(nonCutoffEventsInitial,5), 's', 'Color', moleculeColor );
            catch
                semilogx( poreEventData.detectedEvents_ms(: ,4), poreEventData.detectedEvents_ms(:,5), '.', 'Color', moleculeColor );
            end
            
            if( i == 1 )
                hold on;
            end
        end
    else
        try
            nonCutoffEvents = find( poreEventData.detectedEvents_ms(:,8) == 0 );
            semilogx( poreEventData.detectedEvents_ms(nonCutoffEvents ,4), poreEventData.detectedEvents_ms(nonCutoffEvents,5), 'b.' );
        catch
            semilogx( poreEventData.detectedEvents_ms(: ,4), poreEventData.detectedEvents_ms(:,5), 'b.' );
        end

        hold on;
    end
    
    try
        cutoffEvents = find( poreEventData.detectedEvents_ms(:,8) == 1);
        semilogx( poreEventData.detectedEvents_ms(cutoffEvents ,4), poreEventData.detectedEvents_ms(cutoffEvents,5), 'r.' );
    catch
    end

    hold off;
end

% grid on;
%get number of fishing events
windowFilename = makeWindowTitleString(hObject, eventdata, handles);

if( handles.fishingFile == 1 )
    %nonCutoffEvents = find( (poreEventData.detectedEvents_ms(:,8) == 0) & (poreEventData.detectedEvents_ms(:,10) == 3));
    nonCutoffEvents = find( (poreEventData.detectedEvents_ms(:,8) == 0) & (poreEventData.detectedEvents_ms(:,9) ~= 0) );
    numFishEvents = length( poreEventData.detectedEvents_ms(nonCutoffEvents ,4) )
    title( { [windowFilename, ' - ', get(handles.edtABFDescription, 'String')]; [ 'baseline:  ', num2str(handles.baseline), ' pA', ',  ', num2str( numFishEvents ), ' total fishing events, analyzed  ', date]} );
else	
    title( { [windowFilename, ' - ', get(handles.edtABFDescription, 'String')]; [ 'baseline:  ', num2str(handles.baseline), ' pA', ',  ', num2str( sum(handles.numDetectedEvents) ), ' total events, analyzed  ', date]} );
end

xlabel('Mean Dwell Time (ms)');
ylabel('Amplitude (pA)');
%axis([0.02 20000 10 40]);

axis([ handles.plotEventWidthMin, handles.plotEventWidthMax, handles.plotEventHeightMin, handles.plotEventHeightMax ]);

%save figure as eps
maximize(gcf);
xlabel('Mean Dwell Time (10^x ms)');
set(gca, 'XScale', 'log')


handles.dateAnalyzed = date;
handles.dataDescription = get(handles.edtABFDescription, 'String');

descriptionLength = length( handles.dataDescription );
if( descriptionLength >= 50 )
    filenameDescriptLength = 50;
elseif( descriptionLength < 50 )
    filenameDescriptLength = descriptionLength;
end

windowFilename = makeWindowTitleString(hObject, eventdata, handles);
        
saveFilename = [windowFilename, ' - analyzed ', handles.dateAnalyzed, '_', handles.dataDescription(1:filenameDescriptLength)];
exportfig(gcf, [saveFilename, '.eps'], 'bounds', 'tight', 'FontSize', 1.8, 'Color', 'rgb');
eps2pdf([saveFilename, '.eps'], 'C:\gs\bin\gswin32.exe')
%delete eps file and leave pdf
delete([saveFilename, '.eps']);
close(gcf);

enableButtons(handles);

% Update handles structure
guidata(hObject, handles);


%% --- Executes on button press in btnExportTrace.
function btnExportTrace_Callback(hObject, eventdata, handles)
% hObject    handle to btnExportTrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%update baseline from edit box
handles.baseline = str2double( get( handles.edtBaseline, 'String') );

%init figure
handles.axesTracePlot = figure(2);

%fprintf('[%f:%f:%f] %f total points\n', handles.plotIndexBegin, handles.plotDelta, handles.plotIndexEnd, length( handles.traceData(handles.plotIndexBegin:handles.plotDelta:handles.plotIndexEnd) ) );

%update baseline from edit box
handles.baseline = str2double( get( handles.edtBaseline, 'String') );


%plot voltage signal if present
if( handles.numSignals > 1 )

    hVoltage = subplot(handles.numSignals,1,2);
    plot( (handles.timeVector(handles.plotIndexBegin:handles.plotDelta:handles.plotIndexEnd)-handles.plotWidthMin)*1000, handles.traceData(handles.plotIndexBegin:handles.plotDelta:handles.plotIndexEnd,2), 'r-', 'HitTest', 'off' );

    grid on;
    ylabel('mV');
    xlabel('Time (ms)');

    axis([(handles.plotWidthMin-handles.plotWidthMin)*1000, (handles.plotWidthMax-handles.plotWidthMin)*1000, handles.plotVHeightMin, handles.plotVHeightMax]);

    %hold off
end

%plot first sweep
%make current plot active
%plot current signal
hCurrent = subplot(handles.numSignals,1,1);
plot( (handles.timeVector(handles.plotIndexBegin:handles.plotDelta:handles.plotIndexEnd)-handles.plotWidthMin).*1000, handles.traceData(handles.plotIndexBegin:handles.plotDelta:handles.plotIndexEnd,1), 'r-', 'HitTest', 'off' );
ylabel('pA');

if( handles.numSignals == 1 )    
    xlabel('Time (ms)');
end

axis([(handles.plotWidthMin-handles.plotWidthMin)*1000, (handles.plotWidthMax-handles.plotWidthMin)*1000, handles.plotCHeightMin, handles.plotCHeightMax]);
grid on;
hold on;

%     %plot calculated baseline
%     if( handles.baseline > handles.plotCHeightMin && handles.baseline < handles.plotCHeightMax)
%         handles.lineBaseline = line( [handles.plotWidthMin handles.plotWidthMax], [handles.baseline handles.baseline], 'LineStyle', '--', 'Color', 'k' );
%     end

%%%%%%%%%%%%%%%%%%%%%%%%
%loop over data and apply exponential mean filter
%calculate the exponentially weighted moving average
if( isequal(get(handles.chkShowMean, 'Value'), 1) )
    subplot(handles.numSignals,1,1), plot( (handles.timeVector(handles.plotIndexBegin:handles.plotDelta:handles.plotIndexEnd)-handles.plotWidthMin).*1000, handles.meanArray(handles.plotIndexBegin:handles.plotDelta:handles.plotIndexEnd,1), 'b', 'HitTest', 'off' );
end
%end exp filter
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate finite difference
if( isequal(get( handles.chkShowDiff, 'Value'), 1) )
    subplot(handles.numSignals,1,1),plot( (handles.timeVector(handles.plotIndexBegin:handles.plotDelta:handles.plotIndexEnd)-handles.plotWidthMin).*1000, handles.diffArray(handles.plotIndexBegin:handles.plotDelta:handles.plotIndexEnd,1), 'g-', 'HitTest', 'off' );
end
%end finite difference
%%%%%%%%%%%%%%%%%%%%%%%%%%

%resize figures (values are percentages of figure window) 
%and adjust axes label (e.g. remove it) 
if( handles.numSignals > 1 )
     set(hCurrent, 'Position',[0.1300    0.3406    0.7750    0.5844], 'XTickLabel','');
     set(hVoltage, 'Position',[0.1300    0.1100    0.7750    0.1706]); 
end

hold off;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot detected events if previously calculated
% if( isequal(handles.sweepsAnalyzed( handles.currentSweep ), 1) )
%     %use single sweep analyze callback
%     btnAnalyze_Callback(hObject, eventdata, handles)
%     %get changes to handles made in btnAnalyze_Callback
%     handles = guidata( hObject );
% else
%     set( handles.txtCurrentPlot, 'String', 'Current Signal Trace (Not Analyzed)' );
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%save figure as eps
maximize(gcf);

handles.dateAnalyzed = date;
handles.dataDescription = get(handles.edtABFDescription, 'String');

descriptionLength = length( handles.dataDescription );
if( descriptionLength >= 50 )
    filenameDescriptLength = 50;
elseif( descriptionLength < 50 )
    filenameDescriptLength = descriptionLength;
end

windowFilename = makeWindowTitleString(hObject, eventdata, handles);
saveFilename = [windowFilename, ' - analyzed ', handles.dateAnalyzed, '_', num2str(handles.plotWidthMin), 'ms-', num2str(handles.plotWidthMax), 'ms Sweep ', num2str(handles.currentSweep), ' ', handles.dataDescription];
exportfig(gcf, [saveFilename, '.eps'], 'bounds', 'tight', 'FontSize', 1.8, 'Color', 'rgb');
eps2pdf([saveFilename, '.eps'], 'C:\gs\bin\gswin32.exe')
%delete eps file and leave pdf
delete([saveFilename, '.eps']);
close(handles.axesTracePlot);

%save trace data to mat file
timingInfo = handles.timeVector(handles.plotIndexBegin:handles.plotIndexEnd)';
currentData = handles.traceData(handles.plotIndexBegin:handles.plotIndexEnd,1);
if( handles.numSignals >= 2 )
    voltageData = handles.traceData(handles.plotIndexBegin:handles.plotIndexEnd,2);
end
save([saveFilename '.mat'], 'timingInfo', 'currentData', 'voltageData');
clear timingInfo;
clear currentData;
clear voltageData;


enableButtons(handles);

% Update handles structure
guidata(hObject, handles);


%% --- Executes on button press in btnAnalyze.
function btnAnalyze_Callback(hObject, eventdata, handles)
% hObject    handle to btnAnalyze (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%update baseline from edit box
handles.baseline = str2double( get( handles.edtBaseline, 'String') );

%set flag to show either events or current trace depending on what is
%visible
if( strcmp( get( handles.axes6,  'Visible'), 'on' ) )
    showEventPlot = 1;
    axes( handles.axes6 );
else
    showEventPlot = 0;
end

%extract events out of current sweep only if they haven't been analyzed
%before
if( isequal(handles.sweepsAnalyzed{handles.currentFileIndex}(handles.currentSweep), 0) )
    handles.minEventTime = str2double( get( handles.edtMinEventTime, 'String' ) );
    handles.holdingVoltage = str2double( get( handles.edtHoldingVoltage, 'String' ) );
    handles.probingVoltage = str2double( get( handles.edtProbingVoltage, 'String' ) );
    handles.ejectVoltage = str2double( get( handles.edtEjectVoltage, 'String' ) );
    
    if isnan( handles.minEventTime )
        handles.minEventTime = 0;
    end
    
    if( handles.fishingFile == 0 ) handles.fishingBaseline = -250; end
    if( handles.rampingFile == 1 ) handles.fishingBaseline = -500; end
    if( handles.primerFile == 1 ) handles.fishingBaseline = -750; end
    
    %check if in multisweep analyze mode
    if( isequal(handles.analyzeEvents, 0) )
       handles.cutoffEventSweep = 0;
       handles.cutoffEventStartTime = 0;
    end
    
    %bound event check to between cursors
    if( get(handles.chkAnalyzeBounds, 'Value') == 1 )        
        %send only amount of trace that is selected
        analyzeStart = round(handles.analyzeLowerBound / handles.samplePeriod);
        analyzeEnd = round(handles.analyzeUpperBound / handles.samplePeriod);
        
        [timingInfo, cutoffEventSweep, cutoffEventStartTime, currentDNAnumber, beginHighOut, cutoffVChangeStart,cutoffEventFoundRamp, cutoffEventFoundPrimer, eventVoltageTimingPrimer] ...
            = parseEvents( handles.traceData(analyzeStart:analyzeEnd,:), handles.samplePeriod, handles.baseline, handles.rmsNoise, handles.minEventTime, handles.fishingBaseline, handles.currentSweep, 0, 0, 0, 0, 0, 0, handles.rmsThreshold, handles.holdingVoltage, handles.probingVoltage, handles.ejectVoltage );
        
        %add on cutoff event flag and DNA number for events and
        %classification
         if( ~isempty(timingInfo) )
            timingInfo = [timingInfo, zeros(length(timingInfo(:,1)),4)];
	    
         end
         
         %put in if cutoffEventSweep and cutoffEventStartTime are non-zero
         %then an event has been cutoff, flag accordingly
         if( handles.keepCutoffEvents == 1 && cutoffEventSweep ~= 0 && cutoffEventStartTime ~= 0 )
             timingInfo = [timingInfo; cutoffEventStartTime, analyzeEnd-(analyzeStart-1), cutoffEventSweep, cutoffEventSweep, 1, cutoffVChangeStart, 0, 0, 0, 0, 0, 0, 0];
         end

         if( ~isempty(timingInfo) )
            %re-bias to starting time equal to zero;
            timingInfo(:,1:2) = timingInfo(:,1:2) + (analyzeStart-1);
        end
    else
%         fprintf( '%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n', handles.samplePeriod, handles.baseline, handles.rmsNoise, handles.minEventTime, handles.fishingBaseline, handles.currentSweep, handles.cutoffEventSweep, handles.cutoffEventStartTime, handles.currentDNAnumber, handles.beginHighOut, handles.vChangeStart, handles.rmsThreshold);
        [timingInfo, cutoffEventSweep, cutoffEventStartTime, currentDNAnumber, beginHighOut, cutoffVChangeStart,cutoffEventFoundRamp, cutoffEventFoundPrimer, eventVoltageTimingPrimer] ...
            = parseEvents( handles.traceData, handles.samplePeriod, handles.baseline, handles.rmsNoise, handles.minEventTime, handles.fishingBaseline, handles.currentSweep, handles.cutoffEventSweep, handles.cutoffEventStartTime, handles.currentDNAnumber, handles.beginHighOut, handles.vChangeStart, handles.cutoffEventFoundRamp, handles.cutoffEventVoltageTimingPrimer, handles.rmsThreshold, handles.holdingVoltage, handles.probingVoltage, handles.ejectVoltage );
        
        handles.currentDNAnumber = currentDNAnumber;
        handles.beginHighOut = beginHighOut;
        handles.vChangeStart = cutoffVChangeStart;

        %add on cutoff event flag and DNA number for events
        if( ~isempty(timingInfo) )
            timingInfo = [timingInfo(:,1:4), zeros(length(timingInfo(:,1)),1), timingInfo(:,5:8)];
        end

        %if sweep info is the same, continued sweep event
        if( (isequal( cutoffEventSweep, 0 ) && isequal( cutoffEventStartTime, 0 )) )
            %previous cutoff event possible ended
            %save for amp calculation
            handles.cutoffEventSum = handles.runningSumCutoffEvent;            
            handles.runningSumCutoffEvent = 0;

        elseif( (isequal( cutoffEventSweep, handles.cutoffEventSweep ) && isequal( cutoffEventStartTime, handles.cutoffEventStartTime )) )
            %add sweep sum to running total
            handles.runningSumCutoffEvent = handles.runningSumCutoffEvent + sum( handles.traceData(:,1) );            
            
        else %cutoff event info changes (new cutoff event)
            %previous cutoff event possible ended
            %save for amp calculation
            handles.cutoffEventSum = handles.runningSumCutoffEvent;
            
            %start new running total
            handles.runningSumCutoffEvent = sum( handles.traceData( cutoffEventStartTime-cutoffVChangeStart:end, 1 ) );
            
        end

        handles.cutoffEventSweep = cutoffEventSweep;
        handles.cutoffEventStartTime = cutoffEventStartTime;        
        handles.cutoffEventFoundRamp = cutoffEventFoundRamp;
        handles.cutoffEventVoltageTimingPrimer = cutoffEventFoundPrimer;

    end
else
    %just load the previously saved results
    %get indices of events from current sweep    
    try
        %find what found events are in currentSweep
%         temp = handles.detectedEvents;
%         disp(temp)
       
        sweepEvents = find( ((handles.detectedEvents(:,1) == handles.currentSweep) & (handles.detectedEvents(:,10) == handles.currentFileIndex)) | ((handles.detectedEvents(:,4) == handles.currentSweep) & (handles.detectedEvents(:,10) == handles.currentFileIndex)) | ((handles.detectedEvents(:,5) == handles.currentSweep) & (handles.detectedEvents(:,10) == handles.currentFileIndex)) );
        %load pre-calculated event info into timingInfo
        timingInfo = [handles.detectedEvents( sweepEvents, 2:3 ), handles.detectedEvents( sweepEvents, 4:11 )];
	
        
        %go through
    catch
        lasterror
        timingInfo = [];
    end
end
%set current signal plot active
%axes(handles.axesCurrent);

%loop over detected events for current sweep and annotate figure
[m,n] = size(timingInfo);

%get number of events that end in the current sweep
%use that to associate event with sweep
if( m > 0 )
    numEventsInSweep = length( find(timingInfo(:,4) == handles.currentSweep) );
else
    numEventsInSweep = 0;
end

%update detected event count
set( handles.txtCurrentPlot, 'String', ['Current Signal Trace (', num2str(numEventsInSweep),' detected events)'] );

for i = 1:m

    eventStart = timingInfo(i,1)*handles.samplePeriod;
    eventEnd = timingInfo(i,2)*handles.samplePeriod;
    
    %set start/end of event for current sweep
    eventStartSamples = timingInfo(i,1);
    eventEndSamples = timingInfo(i,2);
  
  
    %check to see if event data has start and stop sweeps (old mat files
    %won't
    if( n >= 4 )
        
%         disp( timingInfo(i,:) )
        
        %calculate dwell time for events ending in this sweep
        if( timingInfo(i,3) ~= timingInfo(i,4) )
            %add in sweep time if event crosses one whole sweep before
            %ending
            if( (timingInfo(i,4)-timingInfo(i,3)) > 1 ) %event continues through whole sweep
                eventDwellSamples = ( (timingInfo(i,4)-timingInfo(i,3)-1).*handles.numSamples + (eventEndSamples) + (handles.numSamples-eventStartSamples+1) );                
                                
            else
                %don't add in extra sweep lengths
                eventDwellSamples = ( (eventEndSamples) + (handles.numSamples-eventStartSamples+1) );
                
            end
            
            %calculate amplitude for cutoff events
            avgAmplitude = (handles.cutoffEventSum + sum( handles.traceData( 1:eventEndSamples,1 ) )) ./ eventDwellSamples;
            
            %calculate dwell time
            eventDwellTime = eventDwellSamples*handles.samplePeriod;
         
%              fprintf('eventDwellTime: %f, sum(traceData): %d, end sweep: %d\n', eventDwellTime, sum( handles.traceData(1:eventEndSamples,1 ) ), timingInfo(i,4) );
        else
            %events are in same sweep
%            avgAmplitude = mean( handles.traceData(eventStartSamples:eventEndSamples, 1) );
	    %iqrAmplitudeIndices = find( (handles.traceData(eventStartSamples:eventEndSamples,1) >= prctile(handles.traceData(eventStartSamples:eventEndSamples,1),25)) & (handles.traceData(eventStartSamples:eventEndSamples,1) <= prctile(handles.traceData(eventStartSamples:eventEndSamples,1),75)) );
            %avgAmplitude = mean( handles.traceData( iqrAmplitudeIndices, 1 ) );
           
       	     avgAmplitude = mean( handles.traceData(eventStartSamples+timingInfo(i,8):eventEndSamples, 1) );

            eventDwellTime = (eventEnd - eventStart);
        end

        %adjust start and end indices for plotting purposes
        if( timingInfo(i,3) < handles.currentSweep )
            eventStartSamples = 1;                  %event starts in previous sweep
        elseif( timingInfo(i,4) > handles.currentSweep )
            eventEndSamples = handles.numSamples;   %event starts in next sweep
        end
    
        %calculate current event length
%         eventDwellTime = ( ((timingInfo(i,4)-timingInfo(i,3)).*handles.numSamples) + (eventEndSamples-eventStartSamples) )*handles.samplePeriod;
    else
%        avgAmplitude = mean( handles.traceData(eventStartSamples:eventEndSamples, 1) );
	    %iqrAmplitudeIndices = find( (handles.traceData(eventStartSamples:eventEndSamples,1) >= prctile(handles.traceData(eventStartSamples:eventEndSamples,1),25)) & (handles.traceData(eventStartSamples:eventEndSamples,1) <= prctile(handles.traceData(eventStartSamples:eventEndSamples,1),75)) );
            %avgAmplitude = mean( handles.traceData( iqrAmplitudeIndices, 1 ) );

	%offset amplitue calculation by the size of the 
	%capacitive transient, timingInfo(i,8)
              
        avgAmplitude = mean( handles.traceData(eventStartSamples+timingInfo(i,8):eventEndSamples, 1) );
        eventDwellTime = (eventEnd - eventStart);
    end

%     if( ~isempty(timingInfo) )
%         fprintf('event [%d,%d], cutoff event: %d\n', timingInfo(i,1), timingInfo(i,2), timingInfo(i,5) )
%     end

    %color detected events
    hold on;  
    if( (isequal(handles.analyzeEvents, 0) || isequal(get(handles.chkReview, 'Value'), 1)) && (showEventPlot == 0) ) %only show events if doing single sweep analysis or review events is on
        if( isequal(get(handles.chkShowMean, 'Value'), 1) )
            %plot mean of detected event
            plot(  handles.timeVector(eventStartSamples:handles.plotDelta:eventEndSamples), handles.meanArray(eventStartSamples:handles.plotDelta:eventEndSamples), 'm-', 'HitTest', 'off');
        else
            if( timingInfo(i,7) == 2  || timingInfo(i,7) == 0  || timingInfo(i,7) == 4 || timingInfo(i,7) == 10 )
                plot(  handles.timeVector(eventStartSamples:handles.plotDelta:eventEndSamples), handles.traceData(eventStartSamples:handles.plotDelta:eventEndSamples, 1), 'b-', 'HitTest', 'off');               
            elseif (timingInfo(i,7) == 11)
                %Plot eject primer events as green
                plot(  handles.timeVector(eventStartSamples:handles.plotDelta:eventEndSamples), handles.traceData(eventStartSamples:handles.plotDelta:eventEndSamples, 1), 'g-', 'HitTest', 'off');               
            else
                %plot raw signal of detected event
                plot(  handles.timeVector(eventStartSamples:handles.plotDelta:eventEndSamples), handles.traceData(eventStartSamples:handles.plotDelta:eventEndSamples, 1), 'm-', 'HitTest', 'off');
            end            
        end
    end
   
    %make matrix of detected event statistics   
    if( isequal(handles.sweepsAnalyzed{handles.currentFileIndex}(handles.currentSweep), 0) ) %if events haven't been detected yet (first time through)
        %plot fitted curve
        %avgAmplitude = plotFittedEventCurve( hObject, eventdata, handles, timingInfo(i,1), timingInfo(i,2) );
        
        if( showEventPlot == 0 )
            line( [handles.timeVector(eventStartSamples), handles.timeVector(eventEndSamples)], [avgAmplitude, avgAmplitude], 'LineStyle','--', 'Color', 'g','HitTest', 'off');
%             line( [handles.timeVector(eventStartSamples+timingInfo(i,8)), handles.timeVector(eventStartSamples+timingInfo(i,8))], [avgAmplitude-5, avgAmplitude+5], 'LineStyle','--', 'Color', 'c','HitTest', 'off');
        else
            %plot event on event plot
            hold on;
            if( timingInfo(i,5) )
                %plot cutoff event red
                semilogx( eventDwellTime*1000, avgAmplitude, 'r.' );
            else
                %get appropriate color based on DNA molecule number
%                 if( timingInfo(i,6) ~= handles.oldDNAnumber  && handles.fishingFile == 1 )
%                     %new DNA molecule, change color
%                     handles.lastMoleculeColor = handles.currentMoleculeColor;
%                     handles.currentMoleculeColor = [ rand, rand, rand ];
%                 else
%                     %color regular events blue
%                     handles.currentMoleculeColor = [ 0, 0, 1 ];
%                 end
                %color events by min amplitude
                if(timingInfo(i,9) < -20)
                    handles.currentMoleculeColor = [0,1,0];
                else
                    handles.currentMoleculeColor = [1,1,0];
                end
         
                %save "old" DNA molecule number
                handles.oldDNAnumber = timingInfo(i,6);
 
                if( (timingInfo(i,7) == 3 )  && (handles.fishingFile == 1) )
                    semilogx( eventDwellTime*1000, avgAmplitude, '.', 'Color', handles.currentMoleculeColor );
                    
%                     if( timingInfo(i,6) > 0 )
%                         timingInfo(i,6) = timingInfo(i,6) - 1;
%                     end
                else
                    semilogx( eventDwellTime*1000, avgAmplitude, 's', 'Color', handles.currentMoleculeColor );

                end
            end
        end
       
        handles.detectedEvents = [handles.detectedEvents; handles.currentSweep, timingInfo(i,1), timingInfo(i,2), timingInfo(i,3), timingInfo(i,4), timingInfo(i,5), timingInfo(i,6), timingInfo(i,7), timingInfo(i,8), handles.currentFileIndex, timingInfo(i,9)];
        handles.detectedEvents_ms = [handles.detectedEvents_ms; handles.currentSweep, eventStart*1000, eventEnd*1000, eventDwellTime*1000, avgAmplitude, timingInfo(i,3), timingInfo(i,4), timingInfo(i,5), timingInfo(i,6), timingInfo(i,7), timingInfo(i,8), handles.currentFileIndex];
        handles.detectedEventsVoltages = [handles.detectedEventsVoltages; eventVoltageTimingPrimer(i,:)];
    else %events have been found, just bplot saved amplitude
        
        if( showEventPlot == 0 )
            %plot saved avg amplitude for each event
            line( [handles.timeVector(eventStartSamples), handles.timeVector(eventEndSamples)], [handles.detectedEvents_ms( sweepEvents(i), 5), handles.detectedEvents_ms( sweepEvents(i), 5)], 'LineStyle','--', 'Color', 'g','HitTest', 'off');
%             line( [handles.timeVector(eventStartSamples+timingInfo(i,8)), handles.timeVector(eventStartSamples+timingInfo(i,8))], [avgAmplitude-5, avgAmplitude+5], 'LineStyle','--', 'Color', 'c','HitTest', 'off');
        else
            %plot event on event plot
            if( handles.detectedEvents(i,5) == handles.currentSweep && handles.detectedEvents(i,11) == handles.currentFileIndex)
                if( handles.detectedEvents_ms(i,7) )
                    %plot cutoff event red
                    semilogx( eventDwellTime*1000, handles.detectedEvents_ms( sweepEvents(i), 5), 'r.' );
                else
                    %get appropriate color based on DNA molecule number
                    if( timingInfo(i,6) ~= handles.oldDNAnumber )
                        %new DNA molecule, change color
                        handles.currentMoleculeColor = [ rand, rand, rand ];
                    end
                    %save "old" DNA molecule number
                    handles.oldDNAnumber = timingInfo(i,6);
                    
                    %plot cutoff event blue
                    semilogx( eventDwellTime*1000, handles.detectedEvents_ms( sweepEvents(i), 5), '.', 'Color', handles.currentMoleculeColor );
                end
            end
        end
    end

    hold off;
 
end %end for loop over detected events

if( isequal(handles.sweepsAnalyzed{handles.currentFileIndex}(handles.currentSweep), 0) )
    %update found event count
    handles.numDetectedEvents(handles.currentFileIndex) = handles.numDetectedEvents(handles.currentFileIndex) + numEventsInSweep;
    set( handles.txtNumEvents, 'String', ['Detected Events: ', num2str(sum(handles.numDetectedEvents))] );

    %mark sweep as analyzed
    handles.sweepsAnalyzed{handles.currentFileIndex}(handles.currentSweep) = 1;
end

%set the current axes as active again
%axes( handles.axesCurrent );

% Update handles structure
guidata(hObject, handles);


%% --- Executes on button press in btnExportData.
function btnExportData_Callback(hObject, eventdata, handles)
% hObject    handle to btnExportData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% disableButtons(handles);

handles.dateAnalyzed = date;
handles.dataDescription = get(handles.edtABFDescription, 'String');

descriptionLength = length( handles.dataDescription );
if( descriptionLength >= 100 )
    filenameDescriptLength = 100;
elseif( descriptionLength < 100 )
    filenameDescriptLength = descriptionLength;
end

windowFilename = makeWindowTitleString(hObject, eventdata, handles);

saveFilename = [windowFilename, ' - analyzed ', handles.dateAnalyzed, '_', handles.dataDescription(1:filenameDescriptLength)];

%copy data and remove trace info
poreEventData = handles;
poreEventData.traceData = [];
poreEventData.timeVector = [];
poreEventData.detectedEvents = handles.detectedEvents;
poreEventData.detectedEvents_ms = handles.detectedEvents_ms;
poreEventData.detectedEventsVoltages = handles.detectedEventsVoltages;
poreEventData.listSweeps = get( handles.lstSweeps, 'String');
poreEventData.listSweepsExcluded = get( handles.lstSweepsExcluded, 'String');

%build excluded sweep numbers
excludedSweeps = get(handles.lstSweepsExcluded, 'String');
excludedSweepsNumbers = zeros( length( excludedSweeps ), 1 );
if(  ~isempty(excludedSweeps) )
    
    for j = 1:handles.numOpenFiles
        for i = 1:length( excludedSweeps )
            excludedSweepInfo = sscanf(excludedSweeps{i}, 'Sweep %d.%d');
            excludedSweepsNumbers(i) = excludedSweepInfo(2);
            excludedSweepsFile(i) = excludedSweepInfo(1);
            poreEventData.sweepsAnalyzed{j}(excludedSweepsNumbers(i)) = 0;
        end
    end

    indsToDel = [];

    for sweepNum = 1:size(excludedSweepsNumbers,1)

        [inds, vals] = find( (handles.detectedEvents(:,1) == excludedSweepsNumbers(sweepNum)) & (handles.detectedEvents(:,10) == excludedSweepsFile(sweepNum)) );
        if( ~isempty( inds ) )
            indsToDel = [indsToDel; inds];
        end
        
    end

    %trim off un-used space due to pre-allocating
    if( ~isempty(indsToDel) )
        poreEventData.detectedEvents(indsToDel,:) = [];
        poreEventData.detectedEvents_ms(indsToDel,:) = [];
        poreEventData.detectedEventsVoltages(indsToDel,:) = [];
    end

end

try
%will need to update to save more data, ok for now
save( [saveFilename, '.mat'], 'poreEventData')
catch
    msgbox( 'Error saving file. Ensure file desription does not contain slashes.', 'Error', 'error' );
end
clear poreEventData

enableButtons(handles);

% Update handles structure
guidata(hObject, handles);



%% --- Executes on button press in btnOpenABF.
function btnOpenABF_Callback(hObject, eventdata, handles)
% hObject    handle to btnOpenABF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axesCurrent);
cla;
axes(handles.axesVoltage);
cla;

%show load ABF file window
set(handles.pnlOpenABFFiles, 'Visible', 'on');
disableButtons(handles);

% Update handles structure
guidata(hObject, handles);


%% --- Executes on button press in chkShowMean.
function chkShowMean_Callback(hObject, eventdata, handles)
% hObject    handle to chkShowMean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkShowMean
%update current sweep and plot
if( handles.firstLoad == 0) loadCurrentSweep( hObject, eventdata, handles ); end
%get changes to handles made in loadCurrentSweep
handles = guidata( hObject );

%% --- Executes on button press in chkShowDiff.

function chkShowDiff_Callback(hObject, eventdata, handles)
% hObject    handle to chkShowDiff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkShowDiff
%update current sweep and plot
if( handles.firstLoad == 0) loadCurrentSweep( hObject, eventdata, handles ); end
%get changes to handles made in loadCurrentSweep
handles = guidata( hObject );


%% --- Executes on change in edtAlpha
function edtAlpha_Callback(hObject, eventdata, handles)
% hObject    handle to edtAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtAlpha as text
%        str2double(get(hObject,'String')) returns contents of edtAlpha as a double


%% --- Executes during object creation, after setting all properties.
function edtAlpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: editcontrols usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%% --- Executes on button press in btnShowEventPlot.
function btnShowEventPlot_Callback(hObject, eventdata, handles)
% hObject    handle to btnShowEventPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%update baseline from edit box
handles.baseline = str2double( get( handles.edtBaseline, 'String') );

%clear figure
cla(handles.axesCurrent);
cla(handles.axesVoltage);

%init figure
%make event plot visible and hide other plots
set(handles.axes6, 'Visible', 'on');
set(handles.axesCurrent, 'Visible', 'off');
set(handles.axesVoltage, 'Visible', 'off');
set(handles.txtCurrentPlot, 'Visible', 'off');
set(handles.txtVoltagePlot, 'Visible', 'off');

axes(handles.axes6);
%set(handles.axes5, 'Position', [200 100 1000 450], 'Name', 'Detected Events');

%poreEventData = handles;
poreEventData.detectedEvents = handles.detectedEvents;
poreEventData.detectedEvents_ms = handles.detectedEvents_ms;
% poreEventData.listSweeps = get( handles.lstSweeps, 'String');
% poreEventData.listSweepsExcluded = get( handles.lstSweepsExcluded, 'String');


%build excluded sweep numbers
excludedSweeps = get(handles.lstSweepsExcluded, 'String');
excludedSweepsNumbers = zeros( length( excludedSweeps ), 1 );
if(  ~isempty(excludedSweeps) )

    for j = 1:handles.numOpenFiles
        for i = 1:length( excludedSweeps )
            excludedSweepInfo = sscanf(excludedSweeps{i}, 'Sweep %d.%d');
            excludedSweepsNumbers(i) = excludedSweepInfo(2);
            excludedSweepsFile(i) = excludedSweepInfo(1);
            poreEventData.sweepsAnalyzed{j}(excludedSweepsNumbers(i)) = 0;
        end
    end

    indsToDel = [];

    %set x axis to log
    set( handles.axes6, 'XScale', 'log' );
    for sweepNum = 1:size(excludedSweepsNumbers,1)
        
        %don't search if no events
        if( min(size( handles.detectedEvents )) > 0 )
            [inds, vals] = find( (handles.detectedEvents(:,1) == excludedSweepsNumbers(sweepNum)) & (handles.detectedEvents(:,10) == excludedSweepsFile(sweepNum)));

            if( ~isempty( inds ) )
                indsToDel = [indsToDel; inds];
            end
        end
                
    end

    %trim off un-used space due to pre-allocating
    if( ~isempty(indsToDel) )
        poreEventData.detectedEvents(indsToDel,:) = [];
        poreEventData.detectedEvents_ms(indsToDel,:) = [];
    end

end

%plot events
set( handles.axes6, 'XScale', 'log' );
[m,n] = size( poreEventData.detectedEvents_ms );
if( m >= 1) %if any events have been detected
    if( n >= 9 ) %dna molecule information included
        numberDNAmolecules = max( poreEventData.detectedEvents_ms(:,9) );
        for i = 0:numberDNAmolecules
            
            if( i == 0 )
                moleculeColor = [0,0,1];
            else
                moleculeColor = [rand, rand, rand];
            end

            try
		if( handles.fishingFile == 1 )
 		    nonCutoffEvents = find( (poreEventData.detectedEvents_ms(:,8) == 0) & (poreEventData.detectedEvents_ms(:,9) == i) & (poreEventData.detectedEvents_ms(:,10) == 3));
%             nonCutoffEvents = find( (poreEventData.detectedEvents_ms(:,8) == 0) & (poreEventData.detectedEvents_ms(:,9) == i) );
 		    nonCutoffEventsInitial = find( (poreEventData.detectedEvents_ms(:,8) == 0) & (poreEventData.detectedEvents_ms(:,9) == i) & (poreEventData.detectedEvents_ms(:,10) ~= 3));
		else
% 		    nonCutoffEvents = find( (poreEventData.detectedEvents_ms(:,8) == 0) & (poreEventData.detectedEvents_ms(:,9) == i) & (poreEventData.detectedEvents_ms(:,10) ~= 2));
% 		    nonCutoffEventsInitial = find( (poreEventData.detectedEvents_ms(:,8) == 0) & (poreEventData.detectedEvents_ms(:,9) == i) & (poreEventData.detectedEvents_ms(:,10) == 2));
             
%             nonCutoffEvents = find( (poreEventData.detectedEvents_ms(:,8) == 0) & (poreEventData.detectedEvents_ms(:,9) == i));
            nonCutoffEvents = find( (poreEventData.detectedEvents_ms(:,8) == 0) & (poreEventData.detectedEvents_ms(:,9) == i) );
		    nonCutoffEventsInitial = [];
	   	end

                semilogx( poreEventData.detectedEvents_ms(nonCutoffEventsInitial ,4), poreEventData.detectedEvents_ms(nonCutoffEventsInitial,5), 's', 'Color', moleculeColor );
                semilogx( poreEventData.detectedEvents_ms(nonCutoffEvents ,4), poreEventData.detectedEvents_ms(nonCutoffEvents,5), '.', 'Color', moleculeColor );
%             	semilogx( poreEventData.detectedEvents_ms(nonCutoffEvents ,4), poreEventData.detectedEvents_ms(nonCutoffEvents,5), '.', 'Color', moleculeColor );
            catch
                semilogx( poreEventData.detectedEvents_ms(: ,4), poreEventData.detectedEvents_ms(:,5), '.', 'Color', moleculeColor );
            end

            hold on;
        end
    else
        try
            nonCutoffEvents = find( poreEventData.detectedEvents_ms(:,8) == 0 );
            semilogx( poreEventData.detectedEvents_ms(nonCutoffEvents ,4), poreEventData.detectedEvents_ms(nonCutoffEvents,5), 'b.' );
        catch
            semilogx( poreEventData.detectedEvents_ms(: ,4), poreEventData.detectedEvents_ms(:,5), 'b.' );
        end
        hold on;
    
    end
    
    try
        cutoffEvents = find( poreEventData.detectedEvents_ms(:,8) == 1);
        semilogx( poreEventData.detectedEvents_ms(cutoffEvents ,4), poreEventData.detectedEvents_ms(cutoffEvents,5), 'r.' );
    catch
    end

    hold off;
end

grid on;
%title( { [handles.filename, ' - ', get(handles.edtABFDescription, 'String')]; [ 'baseline:  ', num2str(handles.baseline), ' pA', ',  ', num2str(handles.numDetectedEvents), ' total events, analyzed  ', date]} );
xlabel('Mean Dwell Time (ms)');
ylabel('Amplitude (pA)');
axis([0.02 20000 0 200]);
% axis([ handles.plotEventWidthMin, handles.plotEventWidthMax, handles.plotEventHeightMin, handles.plotEventHeightMax ]);


%clear poreEventData;
axes(handles.axesCurrent);

% Update handles structure
guidata(hObject, handles);



%% --- Executes on button press in btnAnalyzeList.
function btnAnalyzeList_Callback(hObject, eventdata, handles)
% hObject    handle to btnAnalyzeList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%disableButtons(handles);

%save to return to after event detection
oldCurrentSweep = handles.currentSweep;
oldCurrentFileIndex = handles.currentFileIndex;

%initialize cutoff event numbers
handles.cutoffEventSweep = 0;
handles.cutoffEventStartTime = 0;
handles.currentDNAnumber = 0;
handles.currentMoleculeColor = [0,0,1];
handles.oldDNAnumber = 0;

%currently analyzing events flag
handles.analyzeEvents = 1;
handles.beginHighOut = 0;

%set flag if event plot is visible
if( strcmp( get( handles.axes6,  'Visible'), 'on' ) )
    showEventPlot = 1;
    axes( handles.axes6 );
else
    showEventPlot = 0;
end

%init progressbar
%progressbar(0);

% fprintf('numEvents: %d\n', handles.numEvents);

%Loop over sweeps in lstSweeps
stopBar = 0;

%calculate progress bar location above main window
tempUnits = get(handles.figMain, 'Units');
set(handles.figMain, 'Units', 'pixels');
currentMainWindowPosition = get(handles.figMain,'Position');
screenSize = get(0, 'ScreenSize');
%normalize to be less than 1
handles.progressBarPos = [currentMainWindowPosition(1)/screenSize(3),(currentMainWindowPosition(2)+currentMainWindowPosition(4)+35)/screenSize(4)];
set(handles.figMain, 'Units', tempUnits);

progressbar(0, handles.progressBarPos);

fprintf('number of open files: %d\n', handles.numOpenFiles);
for fileNum = 1:handles.numOpenFiles
    handles.currentFileIndex = fileNum;
    
    for sweepNum = 1:handles.numEvents(handles.currentFileIndex)
        if(stopBar)
            break;
        end
        
        %save previous sweep to look for continuity
        previousSweep = handles.currentSweep;

        %get selected sweep number for lstSweeps
        listSweep = get(handles.lstSweeps, 'String');
        selectedSweep = listSweep{sum(handles.numEvents(1:fileNum-1))+sweepNum};

        currentSweepInfo = sscanf(selectedSweep, 'Sweep %d.%d');
        handles.currentSweep = currentSweepInfo(2);
        handles.currentFileIndex = currentSweepInfo(1);
        
        %if sweeps aren't consecutive
        %zero out previous sweep info
        if( handles.currentSweep-1 ~= previousSweep )
            %fprintf('non-continuous sweep\n');
            handles.cutoffEventSweep = 0;
            handles.cutoffEventStartTime = 0;
            %increment DNA number to not associate across breakss
            handles.currentDNAnumber = handles.currentDNAnumber + 1;
            handles.beginHighOut = 0;
            handles.vChangeStart = 0;
        end

        %disp(sscanf(selectedSweep, 'Sweep %d'))
        %update current sweep and plot if on current plot
        if( showEventPlot == 0 )
            loadCurrentSweep( hObject, eventdata, handles );
            %get changes to handles made in loadCurrentSweep
            handles = guidata( hObject );
        else
            loadTraceData( hObject, eventdata, handles );
            %get changes to handles made in loadCurrentSweep
            handles = guidata( hObject );
        end

        %use single sweep analyze callback
        btnAnalyze_Callback(hObject, eventdata, handles)
        %get changes to handles made in btnAnalyze_Callback
        handles = guidata( hObject );

        %give some time for gui interaction
        pause(.1);
        %double(sweepNum)/double(handles.numEvents)
        %show progressbar
        %sum(handles.numEvents)
%         if( fileNum == 1 )
%             barPercentage = double(sweepNum)/double(sum(handles.numEvents));
%         else
%             barPercentage = double(sum(handles.numEvents(1:fileNum-1)+sweepNum))/double(sum(handles.numEvents));
%         end
%         stopBar = progressbar( double(sum(handles.numEvents(1:fileNum-1)+sweepNum))/double(sum(handles.numEvents)), handles.progressBarPos);
        barPercentage = double(sweepNum)/double(handles.numEvents(fileNum));
        stopBar = progressbar( barPercentage, handles.progressBarPos);

        if(stopBar)
            break;
        end
    end
end

beep;

%restore currentSweep
handles.currentSweep = oldCurrentSweep;
handles.currentFileIndex = oldCurrentFileIndex;

%done analyzing events
handles.analyzeEvents = 0;

%enable buttons again
enableButtons(handles);

if( showEventPlot == 0)
    %update current sweep and plot
    loadCurrentSweep( hObject, eventdata, handles );
    %get changes to handles made in loadCurrentSweep
    handles = guidata( hObject );
end

% for i = 1:handles.currentDNAnumber
%     figure(100+i);
%     hist(  handles.detectedEvents(find( (abs(handles.detectedEvents(:,11)) <= 30) & (handles.detectedEvents(:,7) == i) ),11), 120 );
% end

% Update handles structure
guidata(hObject, handles);



%% --- Executes on button press in chkReview.
function chkReview_Callback(hObject, eventdata, handles)
% hObject    handle to chkReview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkReview

%% --- Executes on change in edtFishingBaseline
function edtFishingBaseline_Callback(hObject, eventdata, handles)
% hObject    handle to edtFishingBaseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtFishingBaseline as text
%        str2double(get(hObject,'String')) returns contents of
%        edtFishingBaseline as a double
handles.fishingBaseline = str2double( get( handles.edtFishingBaseline, 'String') );

%plot fishing baseline
%fprintf('%d\n', ( handles.fishingBaseline > handles.plotCHeightMin && handles.fishingBaseline < handles.plotCHeightMax && handles.fishingFile == 1 ))
if( handles.fishingBaseline > handles.plotCHeightMin && handles.fishingBaseline < handles.plotCHeightMax && handles.fishingFile == 1 )
    try
	set(handles.lineFishingBaseline, 'YData', [handles.fishingBaseline handles.fishingBaseline]);
    catch
        handles.lineBaseline = line( [handles.plotWidthMin handles.plotWidthMax], [handles.baseline handles.baseline], 'LineStyle', '--', 'Color', 'k' );
    end
end


% Update handles structure
guidata(hObject, handles);


%% --- Executes during object creation, after setting all properties.
function edtFishingBaseline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtFishingBaseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%% --- Executes on button press in btnOpenMAT.
function btnOpenMAT_Callback(hObject, eventdata, handles)
% hObject    handle to btnOpenMAT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disableButtons(handles);

%set default search path for open dialog
if( ~exist(handles.defaultPathname) )
    handles.defaultPathname = [pwd '\'];
end
if( ~exist('handles.defaultMATPathname') )
    handles.defaultMATPathname = [pwd '\'];
end

[filename, pathname] = uigetfile('*.mat','Open Trace Data', [handles.defaultMATPathname], 'MultiSelect', 'on');

if isequal(filename,0) || isequal(pathname,0)
    %do nothing, user pressed cancel
    handles.defaultMATPathname = [pwd '\'];
else
    %update file info
    handles.defaultMATPathname = pathname;
    %handles.MATpathname = pathname;
    
    %check if multiple MAT files selected
    if( iscell( filename ) )
       numMATfiles = length(filename);
    %only one file selected
    else
       numMATfiles = 1;
    end

    for currentMATindex = 1:numMATfiles
        if( numMATfiles > 1 )
            load( [pathname filename{currentMATindex}] );
        else
            load( [pathname filename] );
        end
        
        %if multiple files are selected, need to update and combine them
        %into on multi-file compatible MAT file.
        %init cell arrays here:
        if( currentMATindex == 1 )           
            %reload variables
            handles.currentSweep = 1;
            handles.currentFileIndex = 1;
           
            handles.filename = poreEventData.filename;
            handles.pathname = poreEventData.pathname;

            %make filename/pathname cell
            if( ~iscell(handles.filename) )
                filename_temp = handles.filename;
                handles.filename = cell(1);
                handles.filename{1} = filename_temp;
            end

            if( ~iscell(handles.pathname) )
                pathname_temp = handles.pathname;
                handles.pathname = cell(1);
                handles.pathname{1} = pathname_temp;
            end
                
%                 handles.numEvents = poreEventData.numEvents
            
%             if( sum(size(handles.numEvents)) <= 2  )
                handles.numEvents = zeros(length(handles.filename),1);
%                 for k = 1:length(handles.filename)
%                     handles.numEvents(k) = length( poreEventData.listSweeps{k} );
%                 end
%             end
             
            handles.numOpenFiles = length(poreEventData.filename);
   
            handles.baseline = poreEventData.baseline;
            handles.analyzeEvents = poreEventData.analyzeEvents;
            handles.currentSweep = poreEventData.currentSweep;
            
            
            %add current file index to detectedEvents arrays
            %if mat file is old might not have all the columns. just
            %add the proper number of columns to make it the right
            %width and then put filenumber on the end
            %width is 11 for handles.detectedEvents
            %and 12 for handles.detectedEvents_ms
            [m,n] = size(poreEventData.detectedEvents);
            if(n<11)
                handles.detectedEvents = [poreEventData.detectedEvents, zeros(m,9-n), ones(m,1), zeros(m,1)];
            else
                handles.detectedEvents = poreEventData.detectedEvents;
            end
            [m,n] = size(poreEventData.detectedEvents_ms);
            if(n<12)
                handles.detectedEvents_ms = [poreEventData.detectedEvents_ms, zeros(m,11-n), ones(m,1)];
            else
                handles.detectedEvents_ms = poreEventData.detectedEvents_ms;
            end
            
            if( sum(size(handles.numDetectedEvents)) <= 2 )
                handles.numDetectedEvents = zeros(length(handles.filename),1);
                for k = 1:length(handles.filename)
                    handles.numDetectedEvents(k) = length( find( handles.detectedEvents(:,10) == k ) );
                end
            end

            handles.firstLoad = poreEventData.firstLoad;
            handles.firstMeanLoad = poreEventData.firstMeanLoad;
            handles.firstDiffLoad = poreEventData.firstDiffLoad;
            handles.numSamples = poreEventData.numSamples;
            handles.samplePeriod = poreEventData.samplePeriod;

            handles.sweepsAnalyzed = poreEventData.sweepsAnalyzed;
            
            if( ~iscell( handles.sweepsAnalyzed ) )
                %make sweepsAnalyzed cell
                a = handles.sweepsAnalyzed;
                handles.sweepsAnalyzed = cell(1);
                handles.sweepsAnalyzed{handles.currentFileIndex} = a;
            end
            
            handles.numSignals = poreEventData.numSignals;
            handles.rmsNoise = poreEventData.rmsNoise;
            handles.dataDescription = poreEventData.dataDescription;
            handles.minEventTime = poreEventData.minEventTime;
            handles.fishingBaseline = poreEventData.fishingBaseline;
            set( handles.edtFishingBaseline, 'String', num2str( handles.fishingBaseline ) );
            try
                handles.keepCutoffEvents = poreEventData.keepCutoffEvents;
                set( handles.chkKeepCutoffEvents, 'Value', handles.keepCutoffEvents );
            catch
                handles.keepCutoffEvents = 0;
            end

            try
                handles.filtered = poreEventData.filtered;
                if ( handles.filtered == 1 )
                    set(handles.chkFilterSignal, 'Value', 1)
                end
            end

            handles.fishingFile = poreEventData.fishingFile;
            set( handles.chkFishingData, 'Value', handles.fishingFile );

            try
                handles.rampingFile = poreEventData.rampingFile;
            catch
                handles.rampingFile = 0;
            end
            set( handles.chkRampData, 'Value', handles.rampingFile);
            
            try
                handles.primerFile = poreEventData.primerFile;
            catch
                handles.primerFile = 0;
            end
            
            set( handles.chkPrimerData, 'Value', handles.primerFile);
            
            %set dimensions of plot windows
            handles.plotWidthMin = poreEventData.plotWidthMin;
            handles.plotWidthMax = poreEventData.plotWidthMax;
            handles.plotVHeightMin = poreEventData.plotVHeightMin;
            handles.plotVHeightMax = poreEventData.plotVHeightMax;

            handles.plotCHeightMin = poreEventData.plotCHeightMin;
            handles.plotCHeightMax = poreEventData.plotCHeightMax;

            %set dimensions of plot windows
            handles.plotWidthMinDefault = poreEventData.plotWidthMinDefault;
            handles.plotWidthMaxDefault = poreEventData.plotWidthMaxDefault;
            handles.plotVHeightMinDefault = poreEventData.plotVHeightMinDefault;
            handles.plotVHeightMaxDefault = poreEventData.plotVHeightMaxDefault;

            handles.plotCHeightMinDefault = poreEventData.plotCHeightMinDefault;
            handles.plotCHeightMaxDefault = poreEventData.plotCHeightMaxDefault;


            %update text box with calculated baseline
            set( handles.edtBaseline, 'String', num2str( handles.baseline ) );

            %recreate timeVector
            handles.timeVector = (0:handles.numSamples-1).*handles.samplePeriod;
            
            %reload lstSweeps
            set( handles.lstSweeps, 'String', {} );
            sweepsCellList = get( handles.lstSweeps, 'String' );

            for k = 1:length( poreEventData.listSweeps );
                sweepsCellList{length(sweepsCellList)+1} = poreEventData.listSweeps{k};
            end

            %reload lstSweepsExcluded
            set( handles.lstSweepsExcluded, 'String', {} );
            sweepsCellListExcluded = get( handles.lstSweepsExcluded, 'String' );

            for k = 1:length( poreEventData.listSweepsExcluded );
                sweepsCellListExcluded{length(sweepsCellListExcluded)+1} = poreEventData.listSweepsExcluded{k};
            end
            
            %set selected item in lstSweeps to first item
            set(handles.lstSweepsExcluded, 'String', sweepsCellListExcluded, 'Value', 1, 'Enable', 'on');
            set(handles.lstSweeps, 'String', sweepsCellList, 'Value', 1);
            
            %read in lstSweeps
            listSweep = get(handles.lstSweeps, 'String');
            [m,n] = size(listSweep);
          
            %only convert list text formatting old version
            if( ~iscell( poreEventData.sweepsAnalyzed ) )
                %update lstSweeps
                handles.numEvents = m;
                for k = 1:m
                    selectedSweepInfo = sscanf(listSweep{k}, 'Sweep %d');
                    listSweep{k} = ['Sweep ' num2str(currentMATindex) '.' num2str(selectedSweepInfo)];
                end
            else
                %update lstSweeps
                for k = 1:m
                    selectedSweepInfo = sscanf(listSweep{k}, 'Sweep %d.%d');
                    %count number of sweeps from each file
                    handles.numEvents(selectedSweepInfo(1)) = handles.numEvents(selectedSweepInfo(1)) + 1;
                    listSweep{k} = ['Sweep ' num2str(selectedSweepInfo(1)) '.' num2str(selectedSweepInfo(2))];
                end
            end
            %update lstSweeps
            set( handles.lstSweeps, 'String', listSweep );

            
            %read in lstSweepsExcluded
            listSweep = get(handles.lstSweepsExcluded, 'String');
            [m,n] = size(listSweep);

            %only convert list text formatting old version
            if( ~iscell( poreEventData.sweepsAnalyzed ) )
                %update lstSweepsExcluded
                for k = 1:m
                    selectedSweepInfo = sscanf(listSweep{k}, 'Sweep %d');

                    listSweep{k} = ['Sweep ' num2str(currentMATindex) '.' num2str(selectedSweepInfo)];
                end
            else
                %update lstSweepsExcluded
                for k = 1:m
                    selectedSweepInfo = sscanf(listSweep{k}, 'Sweep %d.%d');
                    
                    listSweep{k} = ['Sweep ' num2str(selectedSweepInfo(1)) '.' num2str(selectedSweepInfo(2))];
                end
            end
            %update lstSweepsExcluded
            set( handles.lstSweepsExcluded, 'String', listSweep );
            
            %load current trace data and update data structures
            %update current sweep and plot
            handles.currentSweep = 1;
            handles.currentFile = 1;

            loadTraceData( hObject, eventdata, handles );
            
            %get changes to handles made in loadCurrentSweep
            handles = guidata( hObject );

            % Update handles structure
            guidata(hObject, handles);
            
        else %for MAT files after first one, append data to existing structures
            
            handles.currentFileIndex = length(handles.filename);
            if( ~iscell(poreEventData.filename) )
                %make filename/pathname cell
                handles.filename{handles.currentFileIndex+1} = poreEventData.filename;

                handles.pathname{handles.currentFileIndex+1} = poreEventData.pathname;
                
                handles.numEvents = [handles.numEvents; length( poreEventData.listSweeps )];
                handles.numDetectedEvents = [handles.numDetectedEvents; length( poreEventData.detectedEvents(:,1) )];
            else
               
                for k = 1:length(poreEventData.filename)
                    handles.filename{length(handles.filename)+1} = poreEventData.filename{k};
                    handles.filename{length(handles.pathname)+1} = poreEventData.pathname{k};
%                     handles.numEvents = [handles.numEvents; length( poreEventData.listSweeps{k} )];                    
                    handles.numDetectedEvents = [handles.numDetectedEvents; length( find( poreEventData.detectedEvents(:,10) == k ) )];
                end
            end
            handles.numOpenFiles = length(handles.filename);
            
            %add current file index to detectedEvents arrays
            %if mat file is old might not have all the columns. just
            %add the proper number of columns to make it the right
            %width and then put filenumber on the end
            %width is 11 for handles.detectedEvents
            %and 12 for handles.detectedEvents_ms
            [m,n] = size(poreEventData.detectedEvents);
            if(n<11)
                poreEventData.detectedEvents = [poreEventData.detectedEvents, zeros(m,9-n), currentMATindex.*ones(m,1), zeros(m,1)];
            end
            [m,n] = size(poreEventData.detectedEvents_ms);
            if(n<12)
                poreEventData.detectedEvents_ms = [poreEventData.detectedEvents_ms, zeros(m,11-n), currentMATindex.*ones(m,1)];
            end
            %concatenate
            handles.detectedEvents_ms = [handles.detectedEvents_ms; poreEventData.detectedEvents_ms];
            handles.detectedEvents = [handles.detectedEvents; poreEventData.detectedEvents];
            handles.numDetectedEvents = [handles.numDetectedEvents; poreEventData.numDetectedEvents];
            
            %check if numSamples is the same.
            %error out if not
            if( handles.numSamples ~= poreEventData.numSamples || handles.numSignals ~= poreEventData.numSignals || handles.samplePeriod ~= poreEventData.samplePeriod)
                fprintf('ABF file parameters do not match across MAT files. Hold on to your hats.\n');
            end
            
            handles.sweepsAnalyzed{currentMATindex} = poreEventData.sweepsAnalyzed;
            
            %reload lstSweeps
%             listSweep = cell(length( poreEventData.listSweeps ) );
%             for k = 1:length( poreEventData.listSweeps );
%                 listSweep{length(sweepsCellList)+1} = poreEventData.listSweeps{k};
%             end
           
            %update formating for lstSweeps
            listSweepHandles = get(handles.lstSweeps, 'String');
            [m,n] = size(poreEventData.listSweeps);
            listSweep = cell(m,n);

            handles.numOpenFiles = length(handles.filename);

            %only convert list text formatting old version
            if( ~iscell( poreEventData.sweepsAnalyzed ) )
                %update lstSweeps
                handles.numEvents(handles.currentFileIndex+1) = m;
                for k = 1:m
                    selectedSweepInfo = sscanf(poreEventData.listSweeps{k}, 'Sweep %d');

                    listSweep{k} = ['Sweep ' num2str(currentMATindex) '.' num2str(selectedSweepInfo)];
                end
            else
                %update lstSweeps
                for k = 1:m
                    selectedSweepInfo = sscanf(poreEventData.listSweeps{k}, 'Sweep %d.%d');

                    %count number of sweeps from each file
                    handles.numEvents(handles.currentFileIndex+selectedSweepInfo(1)) = handles.numEvents(selectedSweepInfo(1)) + 1;

                    listSweep{k} = ['Sweep ' num2str(handles.currentFileIndex+selectedSweepInfo(1)) '.' num2str(selectedSweepInfo(2))];
                end
            end
            %append new list to current list
            listSweep = [listSweepHandles; listSweep];
            %update lstSweeps
            set( handles.lstSweeps, 'String', listSweep );

            
            %read in lstSweepsExcluded
%             listSweep = cell( length( poreEventData.listSweepsExcluded ) );
% 
%             for k = 1:length( poreEventData.listSweepsExcluded );
%                 listSweep{length(sweepsCellListExcluded)+1} = poreEventData.listSweepsExcluded{k};
%             end

            listSweepHandles = get(handles.lstSweepsExcluded, 'String');
            [m,n] = size(poreEventData.listSweepsExcluded);
            listSweep = cell(m,n);

            %only convert list text formatting old version
            if( ~iscell( poreEventData.sweepsAnalyzed ) )
                %update lstSweepsExcluded
                for k = 1:m
                    selectedSweepInfo = sscanf(poreEventData.listSweepsExcluded{k}, 'Sweep %d');

                    listSweep{k} = ['Sweep ' num2str(currentMATindex) '.' num2str(selectedSweepInfo)];
                end
            else
                %update lstSweepsExcluded
                for k = 1:m
                    selectedSweepInfo = sscanf(poreEventData.listSweepsExcluded{k}, 'Sweep %d.%d');

                    listSweep{k} = ['Sweep ' num2str(currentMATindex) '.' num2str(selectedSweepInfo(2))];
                end
            end
            %append new list to current list
            listSweep = [listSweepHandles; listSweep];
            %update lstSweepsExcluded
            set( handles.lstSweepsExcluded, 'String', listSweep );
        
            handles.numOpenFilenames = length(handles.filename);
            
            % Update handles structure
            guidata(hObject, handles);

        end
        
        
        %update file info information in window
        windowFilename = makeWindowTitleString(hObject, eventdata, handles);

        set( handles.figMain, 'Name', ['Nanopore Event Detector - ' windowFilename] );
        set( handles.txtNumSweeps, 'String', ['Total Sweeps: ', num2str(sum(handles.numEvents))] );
        set( handles.txtSamplePeriod, 'String', ['Sample Period: ', num2str(handles.samplePeriod/(1e-6)), ' usec'] );
        set( handles.txtNumEvents, 'String', ['Detected Events: ', num2str(sum(handles.numDetectedEvents))] );
        set( handles.txtEstBaseline, 'String', ['Est. Baseline: ', num2str(handles.baseline) ,' pA'] );
        set( handles.txtRMSNoise, 'String', ['RMS Noise: ', num2str(handles.rmsNoise) ,' pA RMS'] );

        %update description box
        set( handles.edtABFDescription, 'String', handles.dataDescription );
        set( handles.edtMinEventTime, 'String', handles.minEventTime );

%         %reload lstSweeps
%         set( handles.lstSweeps, 'String', {} );
%         sweepsCellList = get( handles.lstSweeps, 'String' );
% 
%         for k = 1:length( poreEventData.listSweeps );
%             sweepsCellList{length(sweepsCellList)+1} = poreEventData.listSweeps{k};
%         end
% 
%         %reload lstSweepsExcluded
%         set( handles.lstSweepsExcluded, 'String', {} );
%         sweepsCellListExcluded = get( handles.lstSweepsExcluded, 'String' );
% 
%         for k = 1:length( poreEventData.listSweepsExcluded );
%             sweepsCellListExcluded{length(sweepsCellListExcluded)+1} = poreEventData.listSweepsExcluded{k};
%         end

        %reload lstOpenFiles
        set( handles.lstOpenFiles, 'String', {} );
        fileList = get( handles.lstOpenFiles, 'String' );

        try
            for k = length( handles.filename ):-1:1
                fileList{length(fileList)+1} = handles.filename{k};
            end
        catch
            fileList{length(fileList)+1} = handles.filename;
        end
        set(handles.lstOpenFiles, 'String', handles.filename, 'Value', 1);

        handles.numOpenFiles = length(handles.filename);

        %set selected item in lstSweeps to first item
        set(handles.lstSweepsExcluded, 'Value', 1, 'Enable', 'on');
        set(handles.lstSweeps, 'Value', 1);
size(handles.detectedEvents)
        %update current sweep and plot
        handles.currentSweep = 1;
        %     handles.currentFile = 1;
        loadCurrentSweep( hObject, eventdata, handles );
        %get changes to handles made in loadCurrentSweep
        handles = guidata( hObject );
    end
end

enableButtons(handles);

% Update handles structure
guidata(hObject, handles);



%% --- Executes on change in edtMinEventTime
function edtMinEventTime_Callback(hObject, eventdata, handles)
% hObject    handle to edtMinEventTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtMinEventTime as text
%        str2double(get(hObject,'String')) returns contents of edtMinEventTime as a double


%% --- Executes during object creation, after setting all properties.
function edtMinEventTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtMinEventTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%% --- Executes on button press in btnClearDetected.
function btnClearDetected_Callback(hObject, eventdata, handles)
% hObject    handle to btnClearDetected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.detectedEvents = [];
handles.detectedEvents_ms = [];
handles.detectedEventsVoltages = [];

%read current sweep from file and plot
%re-allocate numEvents and numEventsTot
tempCurrentFileIndex = handles.currentFileIndex;

%save current file
tempCurrentFileIndex = handles.currentFileIndex;
for j = 1:handles.numOpenFiles
     handles.currentFileIndex = j;
     loadTraceData( hObject, eventdata, handles );

     %get changes to handles made in loadCurrentSweep
     handles = guidata( hObject );

     handles.sweepsAnalyzed{handles.currentFileIndex} = zeros(handles.numEventsTot(handles.currentFileIndex), 1);
end
%restore current file
handles.currentFileIndex = tempCurrentFileIndex;

handles.currentFileIndex = tempCurrentFileIndex;

handles.numDetectedEvents = zeros(handles.numOpenFiles,1);
set( handles.txtNumEvents, 'String', 'Detected Events: 0' );
handles.currentDNAnumber = 0;

% Update handles structure
guidata(hObject, handles);

%update current sweep and plot depending on what is visible
if( strcmp( get( handles.axes6,  'Visible'), 'on' ) )    
    axes( handles.axes6 );
    cla(handles.axes6 );
else
    loadCurrentSweep( hObject, eventdata, handles );
end

%get changes to handles made in loadCurrentSweep
handles = guidata( hObject );

% Update handles structure
guidata(hObject, handles);

%% --- Update handles struct with handles.currentSweep from current file
function loadTraceData( hObject, eventdata, handles )

correctFile = 0;
counter = 0; %error out after 3 times
while( correctFile == 0 && counter < 3)
    
    %set sampleInt to 1 in case of breakage
    sampleInt = 1;
        
    %load sweep from disk
    try
        %check and see if filename is cell-array or not
        if( iscell(handles.filename) )
%             fprintf('current path: %s, file: %s, index: %d\n', handles.pathname{handles.currentFileIndex}, handles.filename{handles.currentFileIndex}, handles.currentFileIndex);
            [handles.traceData, sampleInt, handles.numEventsTot(handles.currentFileIndex)] = import_abf( [handles.pathname{handles.currentFileIndex} handles.filename{handles.currentFileIndex}], handles.currentSweep );        
        else
            [handles.traceData, sampleInt, handles.numEventsTot(handles.currentFileIndex)] = import_abf( [handles.pathname handles.filename], handles.currentSweep );        
            fprintf('should never go here\n');
        end
        
        %hack for weird recording gain issue with data
%         if( abs( rem(handles.traceData(1,2), 1) ) < 1 )
%             handles.traceData(:,2) = handles.traceData(:,2).*100;           
%         end
        
        %fprintf('%s: numSweeps: %d/%d\n', handles.filename{handles.currentFileIndex}, handles.currentSweep, handles.numEventsTot(handles.currentFileIndex));
        
        if ( get(handles.chkFilterSignal, 'Value') == 1 )
            handles.filtered = 1;
            cutoff_freq = 5000/((1/(sampleInt*1e-6))/2);
            [z,p,k] = butter(8,cutoff_freq,'low');
            [s,g] = zp2sos(z,p,k);
            Hb = dfilt.df2sos(s,g);           

%             cutoff_freq = (5000*sampleInt*1e-6*2);
%             Hf = fdesign.lowpass('N,F3db',8,cutoff_freq);
%             Hb = design(Hf,'butter');

            if (handles.currentSweep ~= 1)
                [tempTraceData, sampleInt, numEventsTot(handles.currentFileIndex)] = import_abf( [handles.pathname{handles.currentFileIndex} handles.filename{handles.currentFileIndex}], handles.currentSweep-1);
                tempTraceData = [tempTraceData((end-round(1000/sampleInt)):end,1); handles.traceData(:,1)];
                tempTraceData = filter(Hb, tempTraceData);
                handles.traceData(:,1) = tempTraceData((end-(length(handles.traceData(:,1))-1)):end);
            else
                tempTraceData = [handles.traceData(round(1000/sampleInt):-1:1,1); handles.traceData(:,1)];
                tempTraceData = filter(Hb, tempTraceData);
                handles.traceData(:,1) = tempTraceData((end-(length(handles.traceData(:,1))-1)):end);
            end

        else
            handles.filtered = 0;
        end

        correctFile = 1;

    catch        
        disp('file not found');
        
%       rethrow(lasterror)
        
        %increment counter and try to find file again
        counter = counter + 1;
    
        %set default search path for open dialog
         if( ~exist(handles.defaultPathname)  )
            handles.defaultPathname = [pwd '\'];
         end
         
         if( iscell(handles.filename) ) %in case of multi file
             [filename, pathname] = uigetfile('*.abf',['Please find ', handles.filename{handles.currentFileIndex}], [handles.defaultPathname]);             

             if isequal(filename,0) || isequal(pathname,0)
                 %do nothing, user pressed cancel
                 handles.defaultPathname = [pwd '\'];
                 counter = 3;
             else
                 %update file info
                 handles.defaultPathname = pathname;
                 %             handles.filename{handles.currentFileIndex} = filename;
                 temp_path = handles.pathname;
                 disp( temp_path )

                 handles.pathname{handles.currentFileIndex} = pathname;
             end

         else % in case of single file
             
             [filename, pathname] = uigetfile('*.abf',['Please find ', handles.filename], [handles.defaultPathname]);

             if isequal(filename,0) || isequal(pathname,0)
                 %do nothing, user pressed cancel
                 handles.defaultPathname = [pwd '\'];
                 counter = 3;
             else
                 %update file info
                 handles.defaultPathname = pathname;
%                  handles.filename= filename;
                 handles.pathname = pathname;
             end
         end        
         
    end     
end %while

handles.samplePeriod = sampleInt*1e-6; %convert sampleInt to usec

%update current sweep text box
set( handles.txtCurrentSweep, 'String', ['Current Sweep: ', num2str(handles.currentSweep)] );


% Update handles structure
guidata(hObject, handles);


%% --- Refreshes current/voltage plot
function loadCurrentSweep( hObject, eventdata, handles )

handles.plotDelta = 1;

%set title to show window is drawing
set( handles.txtCurrentPlot, 'String', 'Drawing...' );

%get current sweep data
loadTraceData( hObject, eventdata, handles );

%get changes to handles made in loadCurrentSweep
handles = guidata( hObject ); 

% keyboard
%fprintf('[%d,%d,%d,%d,%d]\n',
%handles.traceData(1,2),handles.traceData(2,2),handles.traceData(3,2),handles.traceData(4,2),handles.traceData(5,2) );


%create time vector for data from sample period

if( isequal(handles.firstLoad, 1) )
    
    handles.numEvents(handles.currentFileIndex) = handles.numEventsTot(handles.currentFileIndex);
    [handles.numSamples, handles.numSignals] = size( handles.traceData );
    
    %set dimensions of plot windows
    handles.plotWidthMin = 0;
    handles.plotWidthMax = length(handles.traceData(:,1))*handles.samplePeriod;
    handles.plotVHeightMin = -250;
    handles.plotVHeightMax = 250;
    
    if (get(handles.chkPrimerData, 'Value') == 1)
        handles.plotCHeightMin = -150;
        handles.plotCHeightMax = 150;
        handles.plotCHeightMinDefault = -150;
        handles.plotCHeightMaxDefault = 150;
    else
        handles.plotCHeightMin = -60;
        handles.plotCHeightMax = 250;
        handles.plotCHeightMinDefault = -60;
        handles.plotCHeightMaxDefault = 250;
    end
    %set dimensions of plot windows
    handles.plotWidthMinDefault = 0;
    handles.plotWidthMaxDefault = length(handles.traceData(:,1))*handles.samplePeriod;
    handles.plotVHeightMinDefault = -250;
    handles.plotVHeightMaxDefault = 250;

    %initialize time vector if first time loading sweep from file
    handles.timeVector = (0:handles.numSamples-1).*handles.samplePeriod;

    %%%%%%%%%%%%%%%%%%
    %estimate baseline
    %plot first order fit to data
    rmsNoiseWindowSize = 512;
    baseLineMax = 90;
    baseLineMin = 30;

    dataMedian = median( handles.traceData(find( (handles.traceData(:,1) < baseLineMax) & (handles.traceData(:,1) > baseLineMin)) ,1) );
    dataIQR = iqr(handles.traceData(:,1));
    detectionThreshold = dataMedian-(dataIQR/2);

    [I] = find( handles.traceData(:,1) > detectionThreshold );
    baseLineTimeVector = handles.timeVector(I)';
    baseLineData = handles.traceData(I,1);
    stdDev = std( baseLineData );
    handles.rmsNoise = sqrt(mean( (handles.traceData(1:rmsNoiseWindowSize,1)-mean(handles.traceData(1:rmsNoiseWindowSize,1))).^2 ));

    try
        [P,S] = polyfit(baseLineTimeVector, baseLineData, 1);
        handles.baseline = P(2)-(stdDev/1.5);
    catch        
    end
    
    if( isnan(handles.baseline) ) handles.baseline = 60; end

    %update text box with calculated baseline
    set( handles.edtBaseline, 'String', num2str( handles.baseline ) );
    set( handles.edtRMS, 'String', num2str(handles.rmsNoise) );
    set( handles.txtEstBaseline, 'String', ['Est. Baseline: ', num2str(handles.baseline) ,' pA'] );
    set( handles.txtRMSNoise, 'String', ['RMS Noise: ', num2str(handles.rmsNoise) ,' pA RMS'] );
    %end estimate baseline
    %%%%%%%%%%%%%%%%%%%%%%%%%
    handles.minEventData = 0.2;
    
    handles.firstLoad = 0;
end

%clear plot
cla(handles.axes6);

%show proper plots
%hide event plot and restore current plot
set(handles.axes6, 'Visible', 'off');
set(handles.axesCurrent, 'Visible', 'on');
set(handles.txtEventPlotTitle, 'Visible', 'off');
set(handles.axesVoltage, 'Visible', 'on');
set(handles.txtCurrentPlot, 'Visible', 'on');
set(handles.txtVoltagePlot, 'Visible', 'on');

set(handles.btnShowEventTrace, 'String', 'Show Event Plot', 'FontSize', 8.0);

%adjust plotDelta based on how much the signal is zoomed
handles.plotDelta = round(((handles.plotWidthMax - handles.plotWidthMin)/handles.samplePeriod)/70000);
if(handles.plotDelta <= 0) handles.plotDelta = 1; end

%use scale factor for changing x axis units when zooming
% if( (handles.plotWidthMax - handles.plotWidthMin) < 1 )
%     handles.plotTimeScaleing = 1000;
% else
%     handles.plotTimeScaling = 1;
% end

%get indices of displayed data
handles.plotIndexBegin = round(handles.plotWidthMin / handles.samplePeriod);
handles.plotIndexEnd = round(handles.plotWidthMax / handles.samplePeriod);
if( isequal(handles.plotIndexBegin, 0) ) handles.plotIndexBegin = 1; end
if( isequal(handles.plotIndexEnd, 0) ) handles.plotIndexEnd = 1; end

%fprintf('[%f:%f:%f] %f total points\n', handles.plotIndexBegin, handles.plotDelta, handles.plotIndexEnd, length( handles.traceData(handles.plotIndexBegin:handles.plotDelta:handles.plotIndexEnd) ) );

%update baseline from edit box
handles.baseline = str2double( get( handles.edtBaseline, 'String') );

if( get(handles.chkReview, 'Value') ~= 1 && isequal(handles.analyzeEvents, 1) )
    %skip plotting of events (speed up processing)
else

    %plot voltage signal if present
    if( handles.numSignals > 1 )
        %set voltage plot active
        axes(handles.axesVoltage);

        plot( handles.timeVector(handles.plotIndexBegin:handles.plotDelta:handles.plotIndexEnd), handles.traceData(handles.plotIndexBegin:handles.plotDelta:handles.plotIndexEnd,2), 'r-', 'HitTest', 'off' );

        grid on;
        ylabel('Amplitude (mV)');
        xlabel('Time (sec)');
        
        axis([handles.plotWidthMin, handles.plotWidthMax, handles.plotVHeightMin, handles.plotVHeightMax]);

        %hold off
    end

    %plot first sweep
    %make current plot active
    axes(handles.axesCurrent);
    %axis([handles.plotWidthMinDefault, handles.plotWidthMaxDefault, handles.plotCHeightMinDefault, handles.plotCHeightMaxDefault]);

    %plot current signal
    %             plot( handles.timeVector, handles.traceData(:,1), 'r-' );
    plot( handles.timeVector(handles.plotIndexBegin:handles.plotDelta:handles.plotIndexEnd), handles.traceData(handles.plotIndexBegin:handles.plotDelta:handles.plotIndexEnd,1), 'r-', 'HitTest', 'off' );
    ylabel('Amplitude (pA)');
    xlabel('Time (sec)');
    axis([handles.plotWidthMin, handles.plotWidthMax, handles.plotCHeightMin, handles.plotCHeightMax]);
    grid on;
    hold on;

    %plot calculated baseline
    if( handles.baseline > handles.plotCHeightMin && handles.baseline < handles.plotCHeightMax)
        handles.lineBaseline = line( [handles.plotWidthMin handles.plotWidthMax], [handles.baseline handles.baseline], 'LineStyle', '--', 'Color', 'k' );
    end
    
    %plot fishing baseline
    if( handles.fishingBaseline > handles.plotCHeightMin && handles.fishingBaseline < handles.plotCHeightMax && handles.fishingFile == 1 )
        handles.lineFishingBaseline = line( [handles.plotWidthMin handles.plotWidthMax], [handles.fishingBaseline handles.fishingBaseline], 'LineStyle', '--', 'Color', 'm' );
    end
    
    %plot analyze bounds
    if( handles.analyzeUpperBound < handles.plotWidthMax && handles.analyzeBounds == 1 )
        handles.lineAnalyzeBound2 = line( [handles.analyzeUpperBound handles.analyzeUpperBound], [handles.plotCHeightMin handles.plotCHeightMax ], 'LineStyle', '--', 'Color', 'c' );
    end
    if( handles.analyzeLowerBound > handles.plotWidthMin && handles.analyzeBounds == 1 )
        handles.lineAnalyzeBound1 = line( [handles.analyzeLowerBound handles.analyzeLowerBound], [handles.plotCHeightMin handles.plotCHeightMax ], 'LineStyle', '--', 'Color', 'c' );    
    end
    
    set( handles.txtCursorDiff, 'String', ['Cursors Diff: ', num2str( (handles.analyzeUpperBound-handles.analyzeLowerBound)*1000 ) ,' ms'] );


    %%%%%%%%%%%%%%%%%%%%%%%%
    %loop over data and apply exponential mean filter
    %calculate the exponentially weighted moving average
    if( isequal(get(handles.chkShowMean, 'Value'), 1) || isequal(get( handles.chkShowDiff, 'Value'), 1) )
        %axes(handles.axesCurrent);
        
        %check mean box is diff is checked
        if( isequal(get( handles.chkShowDiff, 'Value'), 1) )
            set( handles.chkShowMean, 'Value', 1 );
        end
        
        
        currentMean = handles.traceData(1,1);

        if( isequal(handles.firstMeanLoad, 1) )
            handles.meanArray = zeros(size(handles.traceData(:,1)));
            handles.firstMeanLoad = 0;
        end

        alpha = str2double( get( handles.edtAlpha, 'String' ) );

        for i = 1:handles.numSamples
            currentMean = (1-alpha).*handles.traceData(i,1) + alpha*currentMean;
            handles.meanArray(i) = currentMean;
        end

        plot( handles.timeVector(handles.plotIndexBegin:handles.plotDelta:handles.plotIndexEnd), handles.meanArray(handles.plotIndexBegin:handles.plotDelta:handles.plotIndexEnd,1), 'b', 'HitTest', 'off' );

	%rmsNoiseWindowSize = 512;
	%fprintf('RMS Noise of low pass filtered signal: %f\n', sqrt(mean( (handles.meanArray(1:rmsNoiseWindowSize)-mean(handles.meanArray(1:rmsNoiseWindowSize))).^2 )) );

    end
    %end exp filter
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %calculate finite difference
    if( isequal(get( handles.chkShowDiff, 'Value'), 1) )
        %axes(handles.axesCurrent);

        if( isequal(handles.firstDiffLoad, 1) )
            handles.diffArray = zeros(size(handles.traceData(:,1)));
            handles.firstDiffLoad = 0;
        end

        diffWindowSize = 40;

        for i = 1:(handles.numSamples-diffWindowSize)
            handles.diffArray(i+diffWindowSize) = sign(( handles.meanArray(i+diffWindowSize) - handles.meanArray(i) ))*( (handles.meanArray(i+diffWindowSize) - handles.meanArray(i))./2 ).^2;
        end

        plot( handles.timeVector(handles.plotIndexBegin:handles.plotDelta:handles.plotIndexEnd), handles.diffArray(handles.plotIndexBegin:handles.plotDelta:handles.plotIndexEnd,1), 'g-', 'HitTest', 'off' );
    end
    %end finite difference
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    hold off;

end %end review events

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot detected events if previously calculated

if( isequal(handles.sweepsAnalyzed{handles.currentFileIndex}(handles.currentSweep), 1) )
    %use single sweep analyze callback
    btnAnalyze_Callback(hObject, eventdata, handles)

    %get changes to handles made in btnAnalyze_Callback
    handles = guidata( hObject );
else
    set( handles.txtCurrentPlot, 'String', 'Current Signal Trace (Not Analyzed)' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set( handles.axesCurrent, 'ButtonDownFcn', @axesCurrent_ButtonDownFcn );

% Update handles structure
guidata(hObject, handles);


%% --- Grays out and disables window's controls
function disableButtons(handles)
%change the mouse cursor to an hourglass
set(handles.figMain,'Pointer','watch');
% uiwait(handles.figMain);

%disable all the buttons so they cannot be pressed
set(handles.edtThreshold,'Enable','off');
set(handles.edtRMS,'Enable','off');
set(handles.chkKeepCutoffEvents,'Enable','off');
set(handles.chkFilterSignal,'Enable','off');
set(handles.btnAnalyze,'Enable','off');
set(handles.btnAnalyzeList,'Enable','off');
set(handles.btnExclude,'Enable','off');
set(handles.btnExportData,'Enable','off');
set(handles.btnExportEvent,'Enable','off');
set(handles.btnExportTrace,'Enable','off');
set(handles.btnInclude,'Enable','off');
set(handles.btnOpenABF,'Enable','off');
set(handles.btnShowEventPlot,'Enable','off');
set(handles.chkReview,'Enable','off');
set(handles.chkShowMean,'Enable','off');
set(handles.chkShowDiff,'Enable','off');
set(handles.edtABFDescription,'Enable','off');
set(handles.edtAlpha,'Enable','off');
set(handles.lstSweeps,'Enable','off');
set(handles.lstSweepsExcluded,'Enable','off');
set(handles.btnOpenMAT,'Enable','off');
set(handles.edtFishingBaseline,'Enable','off');
set(handles.edtBaseline,'Enable','off');
set(handles.edtMinEventTime,'Enable','off');
set(handles.btnClearDetected,'Enable','off');
set(handles.btnResetPlot,'Enable','off');
set(handles.btnShowEventPlot,'Enable','off');
set(handles.chkFishingData,'Enable','off');
set(handles.chkShowCursors,'Enable','off');
set(handles.chkAnalyzeBounds,'Enable','off');
set(handles.btnShowEventTrace,'Enable','off');
set(handles.chkRampData, 'Enable', 'off');
set(handles.chkPrimerData, 'Enable', 'off');
%% --- Re-enables window's controls after a disableButtons() call
function enableButtons(handles)
%change the mouse cursor to an arrow
set(handles.figMain,'Pointer','arrow');
% uiresume(handles.figMain);

%enable all the buttons so they can be pressed
set(handles.edtThreshold,'Enable','on');
set(handles.edtRMS,'Enable','on');
set(handles.chkKeepCutoffEvents,'Enable','on');
set(handles.chkFilterSignal,'Enable','on');
set(handles.btnAnalyze,'Enable','on');
set(handles.btnAnalyzeList,'Enable','on');
set(handles.btnExclude,'Enable','on');
set(handles.btnExportData,'Enable','on');
set(handles.btnExportEvent,'Enable','on');
set(handles.btnExportTrace,'Enable','on');
set(handles.btnInclude,'Enable','on');
set(handles.btnOpenABF,'Enable','on');
set(handles.btnShowEventPlot,'Enable','on');
set(handles.chkReview,'Enable','on');
set(handles.chkShowMean,'Enable','on');
set(handles.chkShowDiff,'Enable','on');
set(handles.edtABFDescription,'Enable','on');
set(handles.edtAlpha,'Enable','on');
set(handles.lstSweeps,'Enable','on');
set(handles.btnOpenMAT,'Enable','on');
set(handles.edtFishingBaseline,'Enable','on');
set(handles.edtBaseline,'Enable','on');
set(handles.edtMinEventTime,'Enable','on');
set(handles.btnClearDetected,'Enable','on');
set(handles.lstSweepsExcluded,'Enable','on');
set(handles.btnResetPlot,'Enable','on');
set(handles.btnShowEventPlot,'Enable','on');
set(handles.chkFishingData,'Enable','on');
set(handles.chkShowCursors,'Enable','on');
set(handles.chkAnalyzeBounds,'Enable','on');
set(handles.btnShowEventTrace,'Enable','on');
set(handles.chkRampData, 'Enable', 'on');
set(handles.chkPrimerData, 'Enable', 'on');

%% --- Keyboard shortcuts for GUI functions
function keyboardShortcuts(src,evnt)
%keyPressFcn automatically takes in two inputs
%src is the object that was active when the keypress occurred
%evnt stores the data for the key pressed

%brings in the handles structure in to the function
handles = guidata(src);

k= evnt.Key; %k is the key that is pressed

%disp(k);

if strcmp(k,'e') %if enter was pressed
    pause(0.01) %allows time to update

    %define hObject as the object of the callback that we are going to use
    %in this case, we are mapping the enter key to the add_pushbutton
    %therefore, we define hObject as the add pushbutton
    %this is done mostly as an error precaution
    hObject = handles.btnExclude;

    %call the add pushbutton callback.
    %the middle argument is not used for this callback
    btnExclude_Callback(hObject, [], handles);
    
elseif strcmp(k,'period') %next sweep
    pause(0.01) %allows time to update

    %define hObject as the object of the callback that we are going to use
    %in this case, we are mapping the enter key to the add_pushbutton
    %therefore, we define hObject as the add pushbutton
    %this is done mostly as an error precaution
    hObject = handles.lstSweeps;
  
    %get selected sweep number for lstSweeps
    listSweep = get(handles.lstSweeps, 'String');
    indexSweep = get(handles.lstSweeps, 'Value');

    if( indexSweep < length( listSweep ) )
        set(handles.lstSweeps, 'Value', indexSweep + 1);
    end

    %call the add pushbutton callback.
    %the middle argument is not used for this callback
    lstSweeps_Callback(hObject, [], handles);
       
        
elseif strcmp(k,'comma') %prev sweep
    pause(0.01) %allows time to update
    
    %define hObject as the object of the callback that we are going to use
    %in this case, we are mapping the enter key to the add_pushbutton
    %therefore, we define hObject as the add pushbutton
    %this is done mostly as an error precaution
    hObject = handles.lstSweeps;
        
    %get selected sweep number for lstSweeps
    %listSweep = get(handles.lstSweeps, 'String');
    indexSweep = get(handles.lstSweeps, 'Value');
    
    if( indexSweep > 1)
        set(handles.lstSweeps, 'Value', indexSweep - 1);
    end
   
    %call the add pushbutton callback.
    %the middle argument is not used for this callback
    lstSweeps_Callback(hObject, [], handles);
end

%% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figMain_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%if( handles.onFigure == 1 )
curPos = get( handles.axesCurrent, 'CurrentPoint' );
curPosV = get( handles.axesVoltage, 'CurrentPoint' );
curPosEvent = get( handles.axes6, 'CurrentPoint' );


% positionStr = ['(' num2str( curPos(1,1) ) ', ' num2str( curPos(1,2) ) ')'];
% set(handles.text17, 'String', positionStr );

handles.buttonUp = [ curPos(1,1) curPos(1,2) ];
handles.buttonUpV = [curPosV(1,1) curPosV(1,2)];
handles.buttonUpEvent = [curPosEvent(1,1) curPosEvent(1,2)];


%if position didn't move
if( handles.buttonUp == handles.buttonDown )
    handles.onFigure = 0;
    
    try
        %remove cursor line
        delete(handles.resizeLine1);
        delete(handles.resizeLine2);
    catch
        %ignore if doesn't exist
    end
end

currentVisible = get( handles.axesCurrent, 'Visible');

%only have resize if current trace is visible
if( strcmp( currentVisible, 'on') )
    %if vertical current resize
    if( (handles.onFigure == 1 ) && (curPos(1,1) < handles.plotWidthMin) && (curPos(1,1) > (handles.plotWidthMin - (.07*abs(handles.plotWidthMax-handles.plotWidthMin))))&& (curPos(1,2) > handles.plotCHeightMin) && (curPos(1,2) < handles.plotCHeightMax) )
        set( handles.figMain, 'Pointer', 'top' );

        %move line to where mouse button is released
        set( handles.resizeLine2, 'YData', [curPos(1,2), curPos(1,2)] );

        if( handles.buttonDown(2) < handles.buttonUp(2) )
            handles.plotCHeightMin = handles.buttonDown(2);
            handles.plotCHeightMax = handles.buttonUp(2);
        else
            handles.plotCHeightMin = handles.buttonUp(2);
            handles.plotCHeightMax = handles.buttonDown(2);
        end

        %update current sweep and plot
        if( handles.firstLoad == 0) loadCurrentSweep( hObject, eventdata, handles ); end
        %get changes to handles made in loadCurrentSweep
        handles = guidata( hObject );

        %if horizontal current resize
    elseif( (handles.onFigure == 2) && (curPos(1,2) < handles.plotCHeightMin) && (curPos(1,2) > (handles.plotCHeightMin - (.1*abs(handles.plotCHeightMax-handles.plotCHeightMin)))) && (curPos(1,1) > handles.plotWidthMin) && (curPos(1,1) < handles.plotWidthMax) )
        set( handles.figMain, 'Pointer', 'left' );

        %move line to where mouse button is released
        set( handles.resizeLine2, 'XData', [curPos(1,1), curPos(1,1)] );

        if( handles.buttonDown(1) < handles.buttonUp(1) )
            handles.plotWidthMin = handles.buttonDown(1);
            handles.plotWidthMax = handles.buttonUp(1);
        else
            handles.plotWidthMin = handles.buttonUp(1);
            handles.plotWidthMax = handles.buttonDown(1);
        end

        %update current sweep and plot
        if( handles.firstLoad == 0) loadCurrentSweep( hObject, eventdata, handles ); end
        %get changes to handles made in loadCurrentSweep
        handles = guidata( hObject );

        %if vertical voltage resize
    elseif( (handles.onFigure == 3) && (curPosV(1,1) < handles.plotWidthMin) && (curPosV(1,1) > (handles.plotWidthMin - (.07*abs(handles.plotWidthMax-handles.plotWidthMin))))&& (curPosV(1,2) > handles.plotVHeightMin) && (curPosV(1,2) < handles.plotVHeightMax) )
        set( handles.figMain, 'Pointer', 'top' );

        %move line to where mouse button is released
        set( handles.resizeLine2, 'YData', [curPosV(1,2), curPosV(1,2)] );

        if( handles.buttonDownV(2) < handles.buttonUpV(2) )
            handles.plotVHeightMin = handles.buttonDownV(2);
            handles.plotVHeightMax = handles.buttonUpV(2);
        else
            handles.plotVHeightMin = handles.buttonUpV(2);
            handles.plotVHeightMax = handles.buttonDownV(2);
        end

        %update current sweep and plot
        if( handles.firstLoad == 0) loadCurrentSweep( hObject, eventdata, handles ); end
        %get changes to handles made in loadCurrentSweep
        handles = guidata( hObject );

        %if horizontal voltage resize
    elseif( (handles.onFigure == 4) && (curPosV(1,2) < handles.plotVHeightMin) && (curPosV(1,2) > (handles.plotVHeightMin - (.4*abs(handles.plotVHeightMax-handles.plotVHeightMin)))) && (curPosV(1,1) > handles.plotWidthMin) && (curPosV(1,1) < handles.plotWidthMax) )
        set( handles.figMain, 'Pointer', 'left' );

        %move line to where mouse button is released
        try
            set( handles.resizeLine2, 'XData', [curPosV(1,1), curPosV(1,1)] );
        catch
            handles.resizeLine2 = line( [handles.plotWidthMin handles.plotWidthMax], [curPos(1,2) curPos(1,2)], 'LineStyle', '--'  );
        end

        if( handles.buttonDownV(1) < handles.buttonUpV(1) )
            handles.plotWidthMin = handles.buttonDownV(1);
            handles.plotWidthMax = handles.buttonUpV(1);
        else
            handles.plotWidthMin = handles.buttonUpV(1);
            handles.plotWidthMax = handles.buttonDownV(1);
        end

        %update current sweep and plot
        if( handles.firstLoad == 0) loadCurrentSweep( hObject, eventdata, handles ); end
        %get changes to handles made in loadCurrentSweep
        handles = guidata( hObject );


        %if baseline change
    elseif( (handles.onFigure == 5) && (curPos(1,2) < handles.plotCHeightMax) && (curPos(1,2) > handles.plotCHeightMin) && (curPos(1,1) > handles.plotWidthMin) && (curPos(1,1) < handles.plotWidthMax) )
        set( handles.figMain, 'Pointer', 'top' );

        %make sure mouse didn't move out of range of figure
        if( (curPos(1,2) < handles.plotCHeightMax) && (curPos(1,2) > handles.plotCHeightMin) )
            %move line to where mouse button is released
            set( handles.lineBaseline, 'YData', [curPos(1,2), curPos(1,2)] );

            %update baseline
            handles.baseline = curPos(1,2);


            %update text boxes
            set( handles.edtBaseline, 'String', num2str( handles.baseline ) );
            set( handles.txtEstBaseline, 'String', ['Est. Baseline: ', num2str(handles.baseline) ,' pA'] );
            %     else
            %         %move line back to original spot
            %         set( handles.lineBaseline, 'YData', [handles.baseline handles.baseline] );
        end

        %if fishing baseline change
    elseif( (handles.onFigure == 6) && (curPos(1,2) < handles.plotCHeightMax) && (curPos(1,2) > handles.plotCHeightMin) && (curPos(1,1) > handles.plotWidthMin) && (curPos(1,1) < handles.plotWidthMax) )
        set( handles.figMain, 'Pointer', 'top' );

        %make sure mouse didn't move out of range of figure
        if( (curPos(1,2) < handles.plotCHeightMax) && (curPos(1,2) > handles.plotCHeightMin) )
            %move line to where mouse button is released
            set( handles.lineFishingBaseline, 'YData', [curPos(1,2), curPos(1,2)] );
            
            %update baseline
            handles.fishingBaseline = curPos(1,2);
            set( handles.edtFishingBaseline, 'String', num2str( handles.fishingBaseline ) );
        end

        %if analzye lower bound change
    elseif( (handles.onFigure == 7) && (curPos(1,1) < handles.plotWidthMax) && (curPos(1,1) > handles.plotWidthMin) && (curPos(1,2) > handles.plotCHeightMin) && (curPos(1,2) < handles.plotCHeightMax) )
        set( handles.figMain, 'Pointer', 'top' );

        %make sure mouse didn't move out of range of figure
        if( (curPos(1,1) < handles.plotWidthMax) && (curPos(1,1) > handles.plotWidthMin) )
            %move line to where mouse button is released
            set( handles.lineAnalyzeBound1, 'XData', [curPos(1,1), curPos(1,1)] );

            %update bound
            handles.analyzeLowerBound = curPos(1,1);

            %fix if lower bound > than upper bound
            if( handles.analyzeLowerBound > handles.analyzeUpperBound )
                tempBound = handles.analyzeUpperBound;
                handles.analyzeUpperBound = handles.analyzeLowerBound;
                handles.analyzeLowerBound = tempBound;

                %swap lines around
                set( handles.lineAnalyzeBound1, 'XData', [handles.analyzeLowerBound, handles.analyzeLowerBound] );
                set( handles.lineAnalyzeBound2, 'XData', [handles.analyzeUpperBound, handles.analyzeUpperBound] );

            end

            set( handles.txtCursorDiff, 'String', ['Cursors Diff: ', num2str( (handles.analyzeUpperBound-handles.analyzeLowerBound)*1000 ) ,' ms'] );
        end

        %if analyze upper bound change
    elseif( (handles.onFigure == 8) && (curPos(1,1) < handles.plotWidthMax) && (curPos(1,1) > handles.plotWidthMin) && (curPos(1,2) > handles.plotCHeightMin) && (curPos(1,2) < handles.plotCHeightMax) )
        set( handles.figMain, 'Pointer', 'top' );

        %make sure mouse didn't move out of range of figure
        if( (curPos(1,1) < handles.plotWidthMax) && (curPos(1,1) > handles.plotWidthMin) )
            %move line to where mouse button is released
            set( handles.lineAnalyzeBound2, 'XData', [curPos(1,1), curPos(1,1)] );

            %update bound
            handles.analyzeUpperBound = curPos(1,1);

            %fix if upper bound < than lower bound
            if( handles.analyzeLowerBound > handles.analyzeUpperBound )
                tempBound = handles.analyzeUpperBound;
                handles.analyzeUpperBound = handles.analyzeLowerBound;
                handles.analyzeLowerBound = tempBound;

                %swap lines around
                set( handles.lineAnalyzeBound1, 'XData', [handles.analyzeLowerBound, handles.analyzeLowerBound] );
                set( handles.lineAnalyzeBound2, 'XData', [handles.analyzeUpperBound, handles.analyzeUpperBound] );

            end
            set( handles.txtCursorDiff, 'String', ['Cursors Diff: ', num2str( (handles.analyzeUpperBound-handles.analyzeLowerBound)*1000 ) ,' ms'] );
        end

        %not a resize event
    else
        set( handles.figMain, 'Pointer', 'arrow' );
        try
            %move line back to original spot
            if( (handles.baseline > handles.plotCHeightMin) && (handles.baseline < handles.plotCHeightMax))
                set( handles.lineBaseline, 'YData', [handles.baseline handles.baseline] );
            end

            if( (handles.fishingFile) && (handles.fishingBaseline > handles.plotCHeightMin) && (handles.fishingBaseline < handles.plotCHeightMax))
                set( handles.lineFishingBaseline, 'YData', [handles.fishingBaseline handles.fishingBaseline] );
            end
	 catch
   	 end

         try
             %remove cursor line
             delete(handles.resizeLine1);
             delete(handles.resizeLine2);
         catch
             %ignore if doesn't exist
         end
	 
    end

%if a eventplot resize event
else

    %if vertical current resize
    if( (handles.onFigure == 9 ) && (curPosEvent(1,1) < handles.plotEventWidthMin) && (curPosEvent(1,1) > (handles.plotEventWidthMin - (.07*abs(handles.plotEventWidthMax-handles.plotEventWidthMin))))&& (curPosEvent(1,2) > handles.plotEventHeightMin) && (curPosEvent(1,2) < handles.plotEventHeightMax) )
        set( handles.figMain, 'Pointer', 'top' );

        %move line to where mouse button is released
        set( handles.resizeLine2, 'YData', [curPosEvent(1,2), curPosEvent(1,2)] );

        if( handles.buttonDownEvent(2) < handles.buttonUpEvent(2) )
            handles.plotEventHeightMin = handles.buttonDownEvent(2);
            handles.plotEventHeightMax = handles.buttonUpEvent(2);
        else
            handles.plotEventHeightMin = handles.buttonUpEvent(2);
            handles.plotEventHeightMax = handles.buttonDownEvent(2);
        end
	
        %update event plot
        if( handles.firstLoad == 0) 
	    btnShowEventPlot_Callback(hObject, eventdata, handles)
	end		
	
        %get changes to handles made in loadCurrentSweep
        handles = guidata( hObject );

        %if horizontal current resize
    elseif( (handles.onFigure == 10) && (curPosEvent(1,2) < handles.plotEventHeightMin) && (curPosEvent(1,2) > (handles.plotEventHeightMin - (.08*abs(handles.plotEventHeightMax-handles.plotEventHeightMin)))) && (curPosEvent(1,1) > handles.plotEventWidthMin) && (curPosEvent(1,1) < handles.plotEventWidthMax) )
        set( handles.figMain, 'Pointer', 'left' );

        %move line to where mouse button is released
        set( handles.resizeLine2, 'XData', [curPosEvent(1,1), curPosEvent(1,1)] );

        if( handles.buttonDownEvent(1) < handles.buttonUpEvent(1) )
            handles.plotEventWidthMin = handles.buttonDownEvent(1);
            handles.plotEventWidthMax = handles.buttonUpEvent(1);
        else
            handles.plotEventWidthMin = handles.buttonUpEvent(1);
            handles.plotEventWidthMax = handles.buttonDownEvent(1);
        end
	
        %update event plot
        if( handles.firstLoad == 0) 
	    btnShowEventPlot_Callback(hObject, eventdata, handles)
	end		

        %get changes to handles made in loadCurrentSweep
        handles = guidata( hObject );

    else
	 set( handles.figMain, 'Pointer', 'arrow' );

         try
             %remove cursor line
             delete(handles.resizeLine1);
             delete(handles.resizeLine2);
         catch
             %ignore if doesn't exist
         end
    end
end

handles.onFigure = 0;

% Update handles structure
guidata(hObject, handles);



%% --- Executes on mouse motion over figure - except title and menu.
function figMain_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 
    curPos = get( handles.axesCurrent, 'CurrentPoint' );
    curPosV = get( handles.axesVoltage, 'CurrentPoint' );
    curPosEvent = get( handles.axes6, 'CurrentPoint' );
    
     
    %positionStr = ['(' num2str( curPos(1,1) ) ', ' num2str( curPos(1,2) ) ')'];
    %set(handles.text17, 'String', positionStr );
  
    try
        currentVisible = get( handles.axesCurrent, 'Visible');

        %only have resize if current trace is visible
        if( strcmp( currentVisible, 'on') )
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %if vertical current resize
            if( (curPos(1,1) < handles.plotWidthMin) && (curPos(1,1) > (handles.plotWidthMin - (.07*abs(handles.plotWidthMax-handles.plotWidthMin))))&& (curPos(1,2) > handles.plotCHeightMin) && (curPos(1,2) < handles.plotCHeightMax) )
                set( handles.figMain, 'Pointer', 'top' );
                %move line to where mouse button is released
                if( handles.onFigure == 1) set( handles.resizeLine2, 'YData', [curPos(1,2), curPos(1,2)] ); end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %if horizontal current resize
            elseif( (curPos(1,2) < handles.plotCHeightMin) && (curPos(1,2) > (handles.plotCHeightMin - (.1*abs(handles.plotCHeightMax-handles.plotCHeightMin)))) && (curPos(1,1) > handles.plotWidthMin) && (curPos(1,1) < handles.plotWidthMax) )
                set( handles.figMain, 'Pointer', 'left' );
                %move line to where mouse button is released
                if( handles.onFigure == 2 ) set( handles.resizeLine2, 'XData', [curPos(1,1), curPos(1,1)] ); end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %if vertical voltage resize
            elseif( (curPosV(1,1) < handles.plotWidthMin) && (curPosV(1,1) > (handles.plotWidthMin - (.07*abs(handles.plotWidthMax-handles.plotWidthMin))))&& (curPosV(1,2) > handles.plotVHeightMin) && (curPosV(1,2) < handles.plotVHeightMax) )
                set( handles.figMain, 'Pointer', 'top' );
                %move line to where mouse button is released
                if( handles.onFigure == 3) set( handles.resizeLine2, 'YData', [curPosV(1,2), curPosV(1,2)] ); end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %if horizontal voltage resize
            elseif( (curPosV(1,2) < handles.plotVHeightMin) && (curPosV(1,2) > (handles.plotVHeightMin - (.4*abs(handles.plotVHeightMax-handles.plotVHeightMin)))) && (curPosV(1,1) > handles.plotWidthMin) && (curPosV(1,1) < handles.plotWidthMax) )
                set( handles.figMain, 'Pointer', 'left' );
                %move line to where mouse button is released
                if( handles.onFigure == 4 ) set( handles.resizeLine2, 'XData', [curPosV(1,1), curPosV(1,1)] ); end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %if baseline move
            elseif( (curPos(1,2) < (handles.baseline + (.02*abs(handles.plotCHeightMax-handles.plotCHeightMin)))) && (curPos(1,2) > (handles.baseline - (.02*abs(handles.plotCHeightMax-handles.plotCHeightMin)))) && (curPos(1,1) > handles.plotWidthMin) && (curPos(1,1) < handles.plotWidthMax) )
                set( handles.figMain, 'Pointer', 'top' );

                %move line to where mouse button is released
                if( handles.onFigure == 5 ) set( handles.lineBaseline, 'YData', [curPos(1,2), curPos(1,2)] ); end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %move line if mousebutton has been pressed
            elseif( handles.onFigure == 5)
                set( handles.lineBaseline, 'YData', [curPos(1,2), curPos(1,2)] );

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %if fishing baseline move
            elseif( (handles.fishingFile == 1) && (curPos(1,2) < (handles.fishingBaseline + (.02*abs(handles.plotCHeightMax-handles.plotCHeightMin)))) && (curPos(1,2) > (handles.fishingBaseline - (.02*abs(handles.plotCHeightMax-handles.plotCHeightMin)))) && (curPos(1,1) > handles.plotWidthMin) && (curPos(1,1) < handles.plotWidthMax) )
                set( handles.figMain, 'Pointer', 'top' );

                %move line to where mouse button is released
                if( handles.onFigure == 6 ) set( handles.lineFishingBaseline, 'YData', [curPos(1,2), curPos(1,2)] ); end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %move line if mousebutton has been pressed
            elseif( (handles.fishingFile == 1) && (handles.onFigure == 6) )
                set( handles.lineFishingBaseline, 'YData', [curPos(1,2), curPos(1,2)] );

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %if lower analyze bounds move
            elseif( (handles.analyzeBounds == 1) && (curPos(1,1) < (handles.analyzeLowerBound + (.02*abs(handles.plotWidthMax-handles.plotWidthMin)))) && (curPos(1,1) > (handles.analyzeLowerBound - (.02*abs(handles.plotWidthMax-handles.plotWidthMin)))) && (curPos(1,2) > handles.plotCHeightMin) && (curPos(1,2) < handles.plotCHeightMax) )
                set( handles.figMain, 'Pointer', 'left' );

                %move line to where mouse button is released
                if( handles.onFigure == 7 )
                    set( handles.lineAnalyzeBound1, 'XData', [curPos(1,1), curPos(1,1)] );
                    set( handles.txtCursorDiff, 'String', ['Cursors Diff: ', num2str( (handles.analyzeUpperBound-handles.analyzeLowerBound)*1000 ) ,' ms'] );
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %move line if mousebutton has been pressed
                %       elseif( (handles.analzyeBounds == 1) && (handles.onFigure == 7) && (curPos(1,2) > handles.plotCHeightMin) && (curPos(1,2) < handles.plotCHeightMax))
            elseif( handles.onFigure == 7 )
                set( handles.lineAnalyzeBound1, 'XData', [curPos(1,1), curPos(1,1)] );
                set( handles.txtCursorDiff, 'String', ['Cursors Diff: ', num2str( (handles.analyzeUpperBound-handles.analyzeLowerBound)*1000 ) ,' ms'] );

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %if upper analyze bounds move
            elseif( (handles.analyzeBounds == 1) && (curPos(1,1) < (handles.analyzeUpperBound + (.02*abs(handles.plotWidthMax-handles.plotWidthMin)))) && (curPos(1,1) > (handles.analyzeUpperBound - (.02*abs(handles.plotWidthMax-handles.plotWidthMin)))) && (curPos(1,2) > handles.plotCHeightMin) && (curPos(1,2) < handles.plotCHeightMax) )
                set( handles.figMain, 'Pointer', 'left' );

                %move line to where mouse button is released
                if( handles.onFigure == 8 )
                    set( handles.lineAnalyzeBound2, 'XData', [curPos(1,1), curPos(1,1)] );
                    set( handles.txtCursorDiff, 'String', ['Cursors Diff: ', num2str( (handles.analyzeUpperBound-handles.analyzeLowerBound)*1000 ) ,' ms'] );
                end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %move line if mousebutton has been pressed
                %         elseif( (handles.analzyeBounds == 1) && (handles.onFigure == 8) && (curPos(1,2) > handles.plotCHeightMin) && (curPos(1,2) < handles.plotCHeightMax))
            elseif( handles.onFigure == 8 )
                set( handles.lineAnalyzeBound2, 'XData', [curPos(1,1), curPos(1,1)] );
                set( handles.txtCursorDiff, 'String', ['Cursors Diff: ', num2str( (handles.analyzeUpperBound-handles.analyzeLowerBound)*1000 ) ,' ms'] );

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %if not in any hot zone
            else
                set( handles.figMain, 'Pointer', 'arrow' );
            end
    %event plot resize
    else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %if vertical current resize
            if( (curPosEvent(1,1) < handles.plotEventWidthMin) && (curPosEvent(1,1) > (handles.plotEventWidthMin - (.01*abs(handles.plotEventWidthMax-handles.plotEventWidthMin))))&& (curPosEvent(1,2) > handles.plotEventHeightMin) && (curPosEvent(1,2) < handles.plotEventHeightMax) )
                set( handles.figMain, 'Pointer', 'top' );
                %move line to where mouse button is released
                if( handles.onFigure == 9) set( handles.resizeLine2, 'YData', [curPosEvent(1,2), curPosEvent(1,2)] ); end

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %if horizontal current resize
            elseif( (curPosEvent(1,2) < handles.plotEventHeightMin) && (curPosEvent(1,2) > (handles.plotEventHeightMin - (.08*abs(handles.plotEventHeightMax-handles.plotEventHeightMin)))) && (curPosEvent(1,1) > handles.plotEventWidthMin) && (curPosEvent(1,1) < handles.plotEventWidthMax) )
                set( handles.figMain, 'Pointer', 'left' );
                %move line to where mouse button is released
                if( handles.onFigure == 10 ) set( handles.resizeLine2, 'XData', [curPosEvent(1,1), curPosEvent(1,1)] ); end
	    else
                set( handles.figMain, 'Pointer', 'arrow' );
			
	    end
	    
        end
    catch
    end

%% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figMain_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figMain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

curPos = get( handles.axesCurrent, 'CurrentPoint' );
curPosV = get( handles.axesVoltage, 'CurrentPoint' );
curPosEvent = get( handles.axes6, 'CurrentPoint' );

currentVisible = get( handles.axesCurrent, 'Visible');

%only have resize if current trace is visible
if( strcmp( currentVisible, 'on') )

    %if vertical current resize
    if( (curPos(1,1) < handles.plotWidthMin) && (curPos(1,1) > (handles.plotWidthMin - (.07*abs(handles.plotWidthMax-handles.plotWidthMin))))&& (curPos(1,2) > handles.plotCHeightMin) && (curPos(1,2) < handles.plotCHeightMax) )
        set( handles.figMain, 'Pointer', 'top' );
        handles.onFigure = 1;

        %draw line denoting resize region
        handles.resizeLine1 = line( [handles.plotWidthMin handles.plotWidthMax], [curPos(1,2) curPos(1,2)], 'LineStyle', '--' );
        handles.resizeLine2 = line( [handles.plotWidthMin handles.plotWidthMax], [curPos(1,2) curPos(1,2)], 'LineStyle', '--'  );

        %if horizontal current resize
    elseif( (curPos(1,2) < handles.plotCHeightMin) && (curPos(1,2) > (handles.plotCHeightMin - (.1*abs(handles.plotCHeightMax-handles.plotCHeightMin)))) && (curPos(1,1) > handles.plotWidthMin) && (curPos(1,1) < handles.plotWidthMax) )
        set( handles.figMain, 'Pointer', 'left' );
        handles.onFigure = 2;

        %draw line denoting resize region
        handles.resizeLine1 = line( [curPos(1,1) curPos(1,1)], [handles.plotCHeightMin handles.plotCHeightMax], 'LineStyle', '--'  );
        handles.resizeLine2 = line( [curPos(1,1) curPos(1,1)], [handles.plotCHeightMin handles.plotCHeightMax], 'LineStyle', '--'  );

        %if vertical voltage resize
    elseif( (curPosV(1,1) < handles.plotWidthMin) && (curPosV(1,1) > (handles.plotWidthMin - (.07*abs(handles.plotWidthMax-handles.plotWidthMin))))&& (curPosV(1,2) > handles.plotVHeightMin) && (curPosV(1,2) < handles.plotVHeightMax) )
        set( handles.figMain, 'Pointer', 'top' );
        handles.onFigure = 3;

        %set voltage plot active
        axes(handles.axesVoltage);

        %draw line denoting resize region
        handles.resizeLine1 = line( [handles.plotWidthMin handles.plotWidthMax], [curPosV(1,2) curPosV(1,2)], 'LineStyle', '--' );
        handles.resizeLine2 = line( [handles.plotWidthMin handles.plotWidthMax], [curPosV(1,2) curPosV(1,2)], 'LineStyle', '--'  );

        %set current plot active again
        axes(handles.axesCurrent);

        %if horizontal voltage resize
    elseif( (curPosV(1,2) < handles.plotVHeightMin) && (curPosV(1,2) > (handles.plotVHeightMin - (.4*abs(handles.plotVHeightMax-handles.plotVHeightMin)))) && (curPosV(1,1) > handles.plotWidthMin) && (curPosV(1,1) < handles.plotWidthMax) )
        set( handles.figMain, 'Pointer', 'left' );
        handles.onFigure = 4;

        %set voltage plot active
        axes(handles.axesVoltage);

        %draw line denoting resize region
        handles.resizeLine1 = line( [curPosV(1,1) curPosV(1,1)], [handles.plotVHeightMin handles.plotVHeightMax], 'LineStyle', '--'  );
        handles.resizeLine2 = line( [curPosV(1,1) curPosV(1,1)], [handles.plotVHeightMin handles.plotVHeightMax], 'LineStyle', '--'  );

        %set current plot active again
        axes(handles.axesCurrent);

        %if baseline
    elseif( (curPos(1,2) < (handles.baseline + (.02*abs(handles.plotCHeightMax-handles.plotCHeightMin)))) && (curPos(1,2) > (handles.baseline - (.02*abs(handles.plotCHeightMax-handles.plotCHeightMin)))) && (curPos(1,1) > handles.plotWidthMin) && (curPos(1,1) < handles.plotWidthMax) )
        set( handles.figMain, 'Pointer', 'top' );
        handles.onFigure = 5;

        %move line to where mouse button is released
        set( handles.lineBaseline, 'YData', [curPos(1,2), curPos(1,2)] );

        %if fishing baseline
    elseif( (handles.fishingFile == 1) && (curPos(1,2) < (handles.fishingBaseline + (.02*abs(handles.plotCHeightMax-handles.plotCHeightMin)))) && (curPos(1,2) > (handles.fishingBaseline - (.02*abs(handles.plotCHeightMax-handles.plotCHeightMin)))) && (curPos(1,1) > handles.plotWidthMin) && (curPos(1,1) < handles.plotWidthMax) )
        set( handles.figMain, 'Pointer', 'top' );
        handles.onFigure = 6;

        %move line to where mouse button is released
        set( handles.lineFishingBaseline, 'YData', [curPos(1,2), curPos(1,2)] );

        %if analzye lower bound
    elseif( (handles.analyzeBounds == 1) && (curPos(1,1) < (handles.analyzeLowerBound + (.02*abs(handles.plotWidthMax-handles.plotWidthMin)))) && (curPos(1,1) > (handles.analyzeLowerBound - (.02*abs(handles.plotWidthMax-handles.plotWidthMin)))) && (curPos(1,2) > handles.plotCHeightMin) && (curPos(1,2) < handles.plotCHeightMax) )
        set( handles.figMain, 'Pointer', 'top' );
        handles.onFigure = 7;

        %move line to where mouse button is released
        set( handles.lineAnalyzeBound1, 'XData', [curPos(1,1), curPos(1,1)] );

        %if analzye upper bound
    elseif( (handles.analyzeBounds == 1) && (curPos(1,1) < (handles.analyzeUpperBound + (.02*abs(handles.plotWidthMax-handles.plotWidthMin)))) && (curPos(1,1) > (handles.analyzeUpperBound - (.02*abs(handles.plotWidthMax-handles.plotWidthMin)))) && (curPos(1,2) > handles.plotCHeightMin) && (curPos(1,2) < handles.plotCHeightMax) )
        set( handles.figMain, 'Pointer', 'top' );
        handles.onFigure = 8;

        %move line to where mouse button is released
        set( handles.lineAnalyzeBound2, 'XData', [curPos(1,1), curPos(1,1)] );

    else
        set( handles.figMain, 'Pointer', 'arrow' );
        handles.onFigure = 0;
    end

%resize/zoom for event plot
else

    %if vertical event plot resize
    if( (curPosEvent(1,1) < handles.plotEventWidthMin) && (curPosEvent(1,1) > (handles.plotEventWidthMin - (.07*abs(handles.plotEventWidthMax-handles.plotEventWidthMin))))&& (curPosEvent(1,2) > handles.plotEventHeightMin) && (curPosEvent(1,2) < handles.plotEventHeightMax) )
        set( handles.figMain, 'Pointer', 'top' );
        handles.onFigure = 9;

        %draw line denoting resize region
        handles.resizeLine1 = line( [handles.plotEventWidthMin handles.plotEventWidthMax], [curPosEvent(1,2) curPosEvent(1,2)], 'LineStyle', '--' );
        handles.resizeLine2 = line( [handles.plotEventWidthMin handles.plotEventWidthMax], [curPosEvent(1,2) curPosEvent(1,2)], 'LineStyle', '--'  );

        %if horizontal event plot resize
    elseif( (curPosEvent(1,2) < handles.plotEventHeightMin) && (curPosEvent(1,2) > (handles.plotEventHeightMin - (.08*abs(handles.plotEventHeightMax-handles.plotEventHeightMin)))) && (curPosEvent(1,1) > handles.plotEventWidthMin) && (curPosEvent(1,1) < handles.plotEventWidthMax) )
        set( handles.figMain, 'Pointer', 'left' );
        handles.onFigure = 10;

        %draw line denoting resize region
        handles.resizeLine1 = line( [curPosEvent(1,1) curPosEvent(1,1)], [handles.plotEventHeightMin handles.plotEventHeightMax], 'LineStyle', '--'  );
        handles.resizeLine2 = line( [curPosEvent(1,1) curPosEvent(1,1)], [handles.plotEventHeightMin handles.plotEventHeightMax], 'LineStyle', '--'  );
    else
        set( handles.figMain, 'Pointer', 'arrow' );
        handles.onFigure = 0;
	
    end	

end

% positionStr = ['(' num2str( curPos(1,1) ) ', ' num2str( curPos(1,2) ) ')'];
% set(handles.text17, 'String', positionStr );

handles.buttonDown = [ curPos(1,1) curPos(1,2) ];
handles.buttonDownV = [ curPosV(1,1) curPosV(1,2) ];
handles.buttonDownEvent = [ curPosEvent(1,1) curPosEvent(1,2) ];

% Update handles structure
guidata(hObject, handles);



%% --- Executes on mouse press over axes background.
function axesCurrent_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesCurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




%% --- Executes on button press in btnResetPlot.
function btnResetPlot_Callback(hObject, eventdata, handles)
% hObject    handle to btnResetPlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%if event plot is visible
if( strcmp( get( handles.axes6,  'Visible'), 'on' ) )
    handles.plotEventWidthMin  = handles.plotEventWidthMinDefault; 
    handles.plotEventWidthMax  = handles.plotEventWidthMaxDefault; 
    handles.plotEventHeightMin = handles.plotEventHeightMinDefault; 
    handles.plotEventHeightMax = handles.plotEventHeightMaxDefault; 
    btnShowEventPlot_Callback(hObject, eventdata, handles)
%if current trace is visible
else
    %reload sweep
    lstSweeps_Callback(hObject, eventdata, handles)
end



%% --- Executes on mouse press over axes background.
function axesVoltage_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axesVoltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% --- Attempts to fit an exponential to the event to better 
% estimate the event amplitude
%fits curve to detected event rather than a line (yeah brett)
function adjustedAmp = plotFittedEventCurve(hObject, eventdata, handles, eventStart, eventEnd)

%3 parameter fit
% k0*( 1 + e^-a1*t - e^-a2*t )

%fit
time = handles.timeVector(eventStart:100:eventEnd);
y = handles.traceData(eventStart:100:eventEnd,1);

b0 = [15, 10, 0.5];
options = statset('MaxIter',1000);
beta_hat = nlinfit(time,y,@fun,b0,options)
y_hat = beta_hat(1) + exp( -beta_hat(2)*time ) + exp( -beta_hat(3)*time.^3 );

%return mean of adjusted signal
adjustedAmp = beta_hat(1);
hold on;
plot(time, y_hat,'cx');
plot(time(1:handles.plotDelta:end), y_hat(1:handles.plotDelta:end),'mo');
%plot(time(1:handles.plotDelta:end), flatend_filt(1:handles.plotDelta:end), 'c--');
%hold on;

%% --- fitting function
function y_hat = fun(b,x)

y_hat = b(1) + exp( -b(2)*x ) + exp( -b(3)*x.^3 );

if( isinf(sum(y_hat)) || isnan(sum(y_hat)))
    y_hat = 10000000.*ones(size(y_hat));
end


%% --- Executes on button press in chkFishingData.
function chkFishingData_Callback(hObject, eventdata, handles)
% hObject    handle to chkFishingData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkFishingData

%draw line for fishing level
handles.fishingFile = get(handles.chkFishingData,'Value');

%update current sweep and plot
loadCurrentSweep( hObject, eventdata, handles );
%get changes to handles made in loadCurrentSweep
handles = guidata( hObject );

% Update handles structure
guidata(hObject, handles);



%% --- Executes on button press in chkShowCursors.
function chkAnalyzeBounds_Callback(hObject, eventdata, handles)
% hObject    handle to chkShowCursors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkShowCursors


%draw line for fishing level
handles.limitAnalysisBetweenCursors = get(handles.chkShowCursors,'Value');

%update current sweep and plot
loadCurrentSweep( hObject, eventdata, handles );
%get changes to handles made in loadCurrentSweep
handles = guidata( hObject );

% Update handles structure
guidata(hObject, handles);


%% --- Executes on button press in chkShowCursors.
function chkShowCursors_Callback(hObject, eventdata, handles)
% hObject    handle to chkShowCursors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkShowCursors

%draw line for fishing level
handles.analyzeBounds = get(handles.chkShowCursors,'Value');

handles.analyzeUpperBound = 0.75*(handles.plotWidthMax-handles.plotWidthMin)+handles.plotWidthMin;
handles.analyzeLowerBound = 0.25*(handles.plotWidthMax-handles.plotWidthMin)+handles.plotWidthMin;

%update current sweep and plot
loadCurrentSweep( hObject, eventdata, handles );
%get changes to handles made in loadCurrentSweep
handles = guidata( hObject );

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in btnShowEventTrace.
function btnShowEventTrace_Callback(hObject, eventdata, handles)
% hObject    handle to btnShowEventTrace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

currentVisible = get( handles.axesCurrent, 'Visible');

if( strcmp(currentVisible, 'on') )
    
    set( handles.figMain, 'Pointer', 'arrow' );
    
     try
         %remove cursor line
         delete(handles.resizeLine1);
         delete(handles.resizeLine2);
     catch
         %ignore if doesn't exist
     end

    %clear plot
    cla(handles.axesCurrent);
    cla(handles.axesVoltage);

    %hide event plot and restore trace plot
    set(handles.axes6, 'Visible', 'on');
    set(handles.txtEventPlotTitle, 'Visible', 'on');
    set(handles.axesCurrent, 'Visible', 'off');
    set(handles.axesVoltage, 'Visible', 'off');
    set(handles.txtCurrentPlot, 'Visible', 'off');
    set(handles.txtVoltagePlot, 'Visible', 'off');
    
    handles.plotEventWidthMin  = handles.plotEventWidthMinDefault; 
    handles.plotEventWidthMax  = handles.plotEventWidthMaxDefault; 
    handles.plotEventHeightMin = handles.plotEventHeightMinDefault; 
    handles.plotEventHeightMax = handles.plotEventHeightMaxDefault; 
    
    set(handles.btnShowEventTrace, 'String', 'Show Current Trace', 'FontSize', 7.0);
    
    %update plot
    btnShowEventPlot_Callback(hObject, eventdata, handles)
else
    
    %clear plot
    cla(handles.axes6)
    
    %hide current plot and restore event plot
    set(handles.axes6, 'Visible', 'off');
    set(handles.txtEventPlotTitle, 'Visible', 'off');    
    set(handles.axesCurrent, 'Visible', 'on');
    set(handles.axesVoltage, 'Visible', 'on');
    set(handles.txtCurrentPlot, 'Visible', 'on');
    set(handles.txtVoltagePlot, 'Visible', 'on');

    
    set(handles.btnShowEventTrace, 'String', 'Show Event Plot', 'FontSize', 8.0);
    
    
    %update current sweep and plot
    lstSweeps_Callback(hObject, eventdata, handles);
    
    %loadCurrentSweep( hObject, eventdata, handles );
end

%get changes to handles made in loadCurrentSweep
handles = guidata( hObject );

% Update handles structure
guidata(hObject, handles);


%%
function edtBaseline_Callback(hObject, eventdata, handles)
% hObject    handle to edtBaseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtBaseline as text
%        str2double(get(hObject,'String')) returns contents of edtBaseline as a double
handles.baseline = str2double( get( handles.edtBaseline, 'String') );
%plot calculated baseline
if( handles.baseline > handles.plotCHeightMin && handles.baseline < handles.plotCHeightMax)
    try
	set( handles.lineBaseline, 'YData', [handles.baseline handles.baseline]);
    catch
        handles.lineBaseline = line( [handles.plotWidthMin handles.plotWidthMax], [handles.baseline handles.baseline], 'LineStyle', '--', 'Color', 'k' );

    end
end
set( handles.txtEstBaseline, 'String', ['Est. Baseline: ', num2str(handles.baseline) ,' pA'] );


% Update handles structure
guidata(hObject, handles);

%% --- Executes during object creation, after setting all properties.
function edtBaseline_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtBaseline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%%
% --- Executes on button press in chkKeepCutoffEvents.
function chkKeepCutoffEvents_Callback(hObject, eventdata, handles)
% hObject    handle to chkKeepCutoffEvents (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkKeepCutoffEvents

%draw line for fishing level
handles.keepCutoffEvents = get(handles.chkKeepCutoffEvents,'Value');

% Update handles structure
guidata(hObject, handles);



%%
function edtRMS_Callback(hObject, eventdata, handles)
% hObject    handle to edtRMS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtRMS as text
%        str2double(get(hObject,'String')) returns contents of edtRMS as a double

handles.rmsNoise = str2double( get( handles.edtRMS, 'String') );
%plot calculated baseline

set( handles.txtRMSNoise, 'String', ['RMS Noise: ', num2str(handles.rmsNoise) ,' pA RMS'] );

% Update handles structure
guidata(hObject, handles);

%%
function edtThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to edtThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtThreshold as text
%        str2double(get(hObject,'String')) returns contents of edtThreshold as a double

handles.rmsThreshold = str2double( get( handles.edtThreshold, 'String') );

% Update handles structure
guidata(hObject, handles);
%%
% --- Executes during object creation, after setting all properties.
function edtRMS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtRMS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%%
% --- Executes during object creation, after setting all properties.
function edtThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%%
% --- Executes on button press in chkFilterSignal.
function chkFilterSignal_Callback(hObject, eventdata, handles)
% hObject    handle to chkFilterSignal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkFilterSignal
% if( get(handles.chkFilterSignal,'Value') )
%     handles.plotCHeightMinDefault = -20;
%     handles.plotCHeightMaxDefault = 50;
%     handles.plotCHeightMin = -20;
%     handles.plotCHeightMax = 50;
%     loadCurrentSweep( hObject, eventdata, handles );
% else
%     handles.plotCHeightMinDefault = 100;
%     handles.plotCHeightMaxDefault = -50;
%     handles.plotCHeightMin = -50;
%     handles.plotCHeightMax = 100;
%     loadCurrentSweep( hObject, eventdata, handles );
% end

if( handles.firstLoad == 0) loadCurrentSweep( hObject, eventdata, handles ); end
%get changes to handles made in loadCurrentSweep
handles = guidata( hObject );

%%
% --- Executes on button press in btnLoadABF.
function btnLoadABF_Callback(hObject, eventdata, handles)
% hObject    handle to btnLoadABF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.currentFileIndex = length( get(handles.lstOpenFiles, 'String') );

%set default search path for open dialog
if( ~exist(handles.defaultPathname)  )
    handles.defaultPathname = [pwd '\'];
end
      
[filename, pathname] = uigetfile('*.abf','Open Trace Data', [handles.defaultPathname],'MultiSelect', 'on');

if isequal(filename,0) || isequal(pathname,0)
    %do nothing, user pressed cancel
    handles.defaultPathname = [pwd '\'];
else
    %update file info, append files to array of open files
    handles.defaultPathname = pathname;
      
    %check if one file selected and if duplicates
    
    if( iscell( filename) )
        lengthFilename = length(filename);
        listFiles = get(handles.lstOpenFiles, 'String');
        for x = 1:lengthFilename            
            duplicateFiles = sum(strcmp(filename{x}, listFiles));
            if( duplicateFiles > 0 ) %found a duplicate
                if( x == 1 )
                    filename = {filename{2:end}};
                    pathname = {pathname{2:end}};   
                 
                elseif( x == lengthFilename )
                    filename = {filename{1:end-1}};
                    pathname = {pathname{1:end-1}};
                 
                else
                    filename = {filename{1:x-1}, filename{x+1,end}};
                    pathname = {pathname{1:x-1}, filename{x+1,end}};
                 
                end
                %update length if files removed
                lengthFilename = lengthFilename - 1;
            end
        end
    else
        %only one file selected
        lengthFilename = 1;
%         duplicateFiles = sum(strcmp(filename, listFiles));
%         if(duplicateFiles > 0)
%             lengthFilename = 0;
%         end
    end
    
    %if any files left after duplicate checking
    if( lengthFilename > 0 )
        for i = 1:lengthFilename
            %catch if only one file selected            
            try
                handles.filename{length(handles.filename)+1} = filename{i};
            catch
                handles.filename{length(handles.filename)+1} = filename;
            end
            handles.pathname{length(handles.pathname)+1} = pathname;

            %get current contents of listbox
            listFiles = get(handles.lstOpenFiles, 'String');
            indexFile = get(handles.lstSweeps, 'Value');
            handles.currentFileIndex = handles.currentFileIndex + 1;
          
            %add new file to list            
            listFiles{length(listFiles)+1} = handles.filename{end};

            %update list
            set(handles.lstOpenFiles, 'String', sort(listFiles) );
            set(handles.lstOpenFiles, 'Value', 1 );                
          

            handles.numOpenFiles = length(listFiles);
            %fprintf('number of open files: %d\n', handles.numOpenFiles);
        end
    else
        handles.currentFileIndex = 0;
    end
   
end

% Update handles structure
guidata(hObject, handles);

%%
% --- Executes on button press in btnRemoveABF.
function btnRemoveABF_Callback(hObject, eventdata, handles)
% hObject    handle to btnRemoveABF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    listFiles = get(handles.lstOpenFiles, 'String');
    indexFile = get(handles.lstOpenFiles, 'Value');

    if( indexFile <= length( listFiles ) )
        selectedFile = listFiles{indexFile};
        
        %convert from list index to handles.filename index
        listIndex = indexFile;
        indexFile = find( strcmp( selectedFile, handles.filename ) == 1 );
     
        %loop over files after removed file and shift items down
        %to fill in gap
        %if first file selected
        if( indexFile == 1 )
            handles.filename = {handles.filename{2:end}};
            handles.pathname = {handles.pathname{2:end}};
            handles.sweepsAnalyzed = {handles.sweepsAnalyzed{2:end}};
            %if last file is selected
        elseif( indexFile == length( listFiles ) )
            handles.filename = {handles.filename{1:end-1}};
            handles.pathname = {handles.pathname{1:end-1}};       
            handles.sweepsAnalyzed = {handles.sweepsAnalyzed{1:end-1}};            
        else
            %if any other file selected
            handles.filename = {handles.filename{1:indexFile-1},handles.filename{indexFile+1:end}};
            handles.pathname = {handles.pathname{1:indexFile-1},handles.pathname{indexFile+1:end}};
            handles.sweepsAnalyzed = {handles.sweepsAnalyzed{1:indexFile-1},handles.sweepsAnalyzed{indexFile+1:end}};
        end
        
        if( indexFile > 1 )
            handles.currentFileIndex = indexFile-1;
        elseif( indexFile == 0 )
            handles.currentFileIndex = 0;
        else
            handles.currentFileIndex = 1;
        end
        
        %search and remove events from file in detectedEvents,
        %detectedEvents_ms, analyzedSweeps, etc
        handles.numEvents = handles.numEvents( find(handles.numEvents ~= indexFile ) );
        handles.numEventsTot = handles.numEventsTot( find(handles.numEventsTot ~= indexFile ) );
        handles.numDetectedEvents = handles.numDetectedEvents( find( handles.numDetectedEvents ~= indexFile ) );
        try
        handles.detectedEvents = handles.detectedEvents( find( handles.detectedEvents(:,10) ~= indexFile ), : );
        handles.detectedEvents_ms = handles.detectedEvents_ms( find( handles.detectedEvents_ms(:,12) ~= indexFile ), : );
        catch
            handles.detectedEvents = [];
            handles.detectedEvents_ms = [];
        end
                
        %update event count
        set( handles.txtNumEvents, 'String', ['Detected Events: ', num2str(sum(handles.numDetectedEvents))] );
                
        %remove element from list
        listFiles(listIndex) = [];

        handles.numOpenFiles = length(listFiles);
        if( handles.numOpenFiles == 0 )
            %zero out cell arrays if nothing left
            handles.filename = {};
            handles.pathname = {};
            listFiles = {};
            handles.currentFileIndex = 0;
        end
        
        %update list
        set(handles.lstOpenFiles, 'String', listFiles );

        %if last element is selected, select new last
        if( listIndex >= length( listFiles ) && listIndex > 1 )
            listIndex = listIndex - 1;
            set(handles.lstOpenFiles, 'Value', listIndex );
        end
        
        %fprintf('number of open files: %d\n', handles.numOpenFiles);
                    
    end

    % Update handles structure
    guidata(hObject, handles);
    
%%
% --- Executes on selection change in lstOpenFiles.
function lstOpenFiles_Callback(hObject, eventdata, handles)
% hObject    handle to lstOpenFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns lstOpenFiles contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstOpenFiles

%%
% --- Executes during object creation, after setting all properties.
function lstOpenFiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstOpenFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%
% --- Executes on button press in btnCloseLoadABF.
function btnCloseLoadABF_Callback(hObject, eventdata, handles)
% hObject    handle to btnCloseLoadABF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    %if files were selected, update eventDetector window
    set( handles.btnCloseLoadABF, 'Enable', 'off' );
 
    handles.currentSweep = 1;
    handles.firstLoad = 1;
    handles.firstMeanLoad = 1;
    handles.firstDiffLoad = 1;
    handles.detectedEvents_ms = [];
    handles.detectedEvents = [];
    handles.detectedEventsVoltages = [];
    handles.axes5 = 0;
    handles.axesTracePlot = 0;
    handles.analyzeEvents = 0;
    handles.numEvents = zeros(handles.numOpenFiles,1);
    handles.baseline = 0;
    handles.fishingBaseline = 0;
    handles.currentFileIndex = 1;

    %read current sweep from file and plot
    %re-allocate numEvents and numEventsTot
    
    handles.numEventsTot = zeros(handles.numOpenFiles,1);
    handles.numEvents = zeros(handles.numOpenFiles,1);
    handles.numDetectedEvents = zeros(handles.numOpenFiles,1);

    %sort handles.filename and handles.pathname
    [tempFile, tempIndices] = sort( handles.filename );
    handles.pathname = handles.pathname(tempIndices);
    handles.filename = handles.filename(tempIndices);
    
    %fprintf('numOpenFiles: %d\n', handles.numOpenFiles);
    if( handles.numOpenFiles > 0 )
        for j = 1:handles.numOpenFiles
            handles.currentFileIndex = j;
            loadTraceData( hObject, eventdata, handles );

            %get changes to handles made in loadCurrentSweep
            handles = guidata( hObject );
            handles.numEvents(handles.currentFileIndex) = handles.numEventsTot(handles.currentFileIndex);
            handles.sweepsAnalyzed{handles.currentFileIndex} = zeros(handles.numEventsTot(handles.currentFileIndex), 1);
        end


        %update file info information in window        
        windowFilename = makeWindowTitleString(hObject, eventdata, handles);
        
        set( handles.figMain, 'Name', ['Nanopore Event Detector - ' windowFilename] );
        set( handles.txtNumSweeps, 'String', ['Total Sweeps: ', num2str(sum(handles.numEvents))] );
        set( handles.txtSamplePeriod, 'String', ['Sample Period: ', num2str(handles.samplePeriod/(1e-6)), ' usec'] );
        set( handles.txtNumEvents, 'String', ['Detected Events: ', num2str(sum(handles.numDetectedEvents))] );


        %zero out lstSweepss
        set( handles.lstSweeps, 'String', {} );
        set( handles.lstSweepsExcluded, 'String', {} );
        sweepsCellList = get( handles.lstSweeps, 'String' );

        %populate sweep list with sweeps from all files
        for fileNo = 1:(handles.numOpenFiles)
            %fprintf('file no.: %d numEvents: %d\n', fileNo, handles.numEvents(fileNo));
            for sweepNo = 1:handles.numEvents(fileNo)
                %fprintf('%d, ', length(sweepsCellList)+1);
                sweepsCellList{sum(handles.numEvents(1:fileNo-1))+sweepNo} = ['Sweep ', num2str(fileNo), '.', num2str(sweepNo)];

            end
        end

        %set selected item in lstSweeps to first item
        set(handles.lstSweeps, 'String', sweepsCellList, 'Value', 1);

        %show current sweep
        handles.currentFileIndex = 1;
        handles.currentSweep = 1;
        loadCurrentSweep( hObject, eventdata, handles );

        % Update handles structure
        guidata(hObject, handles);

    else
        %zero out lstSweepss
        set( handles.lstSweeps, 'String', {} );
        set( handles.lstSweepsExcluded, 'String', {} );
        sweepsCellList = get( handles.lstSweeps, 'String' );
        %set selected item in lstSweeps to first item
        set(handles.lstSweeps, 'String', sweepsCellList, 'Value', 1);

        %update file info information in window
        set( handles.figMain, 'Name', 'Nanopore Event Detector' );
        set( handles.txtNumSweeps, 'String', 'Total Sweeps: ' );
        set( handles.txtSamplePeriod, 'String', 'Sample Period: ' );
        set( handles.txtNumEvents, 'String', 'Detected Events: ');

    end

    
    %show load ABF file window
    set(handles.pnlOpenABFFiles, 'Visible', 'off');
    set( handles.btnCloseLoadABF, 'Enable', 'on' );
    enableButtons(handles);

    % Update handles structure
    guidata(hObject, handles);

%%
function windowFilename = makeWindowTitleString(hObject, eventdata, handles)
%     fprintf('open files: %s\n', char(handles.filename));
    
    %update number of open files in case it's wrong
    handles.numOpenFiles = length(handles.filename);
%     fprintf('number of open files: %d\n', handles.numOpenFiles);
    
    %update file info information in window        
    windowFilename = '';    
%     try
    if( iscell(handles.filename) )
        for i = 1:handles.numOpenFiles
            try
                windowFilename = [windowFilename, ', ', handles.filename{i}];
            catch
                handles.numOpenFiles = i-1;
            end
        end
        windowFilename = windowFilename(3:end);
    else
        windowsFilename = handles.filename;
    end
%     catch
%         windowFilename = char(handles.filename);
%         %fix numOpenFiles
%         handles.numOpenFiles = 1;
%     end
    
    % Update handles structure
    guidata(hObject, handles);


% --- Executes on button press in chkRampData.
function chkRampData_Callback(hObject, eventdata, handles)
% hObject    handle to chkRampData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkRampData
%draw line for fishing level
handles.rampingFile = get(handles.chkRampData,'Value');

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in chkPrimerData.
function chkPrimerData_Callback(hObject, eventdata, handles)
% hObject    handle to chkPrimerData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkPrimerData

handles.primerFile = get(handles.chkPrimerData,'Value');

if( handles.primerFile == 1 )
    set(handles.edtProbingVoltage, 'Visible', 'on');
    set(handles.edtHoldingVoltage, 'Visible', 'on');
    set(handles.edtEjectVoltage, 'Visible', 'on');
    set(handles.txtProbingVoltage, 'Visible', 'on');
    set(handles.txtHoldingVoltage, 'Visible', 'on');
    set(handles.txtEjectVoltage, 'Visible', 'on');

else
    set(handles.edtProbingVoltage, 'Visible', 'off');
    set(handles.edtHoldingVoltage, 'Visible', 'off');
    set(handles.edtEjectVoltage, 'Visible', 'off');
    set(handles.txtProbingVoltage, 'Visible', 'off');
    set(handles.txtHoldingVoltage, 'Visible', 'off');
    set(handles.txtEjectVoltage, 'Visible', 'off');

end
    
handles.firstLoad = 1;

% Update handles structure
guidata(hObject, handles);



function edtHoldingVoltage_Callback(hObject, eventdata, handles)
% hObject    handle to edtHoldingVoltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtHoldingVoltage as text
%        str2double(get(hObject,'String')) returns contents of edtHoldingVoltage as a double


% --- Executes during object creation, after setting all properties.
function edtHoldingVoltage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtHoldingVoltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edtProbingVoltage_Callback(hObject, eventdata, handles)
% hObject    handle to edtProbingVoltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtProbingVoltage as text
%        str2double(get(hObject,'String')) returns contents of edtProbingVoltage as a double


% --- Executes during object creation, after setting all properties.
function edtProbingVoltage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtProbingVoltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edtEjectVoltage_Callback(hObject, eventdata, handles)
% hObject    handle to edtEjectVoltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtEjectVoltage as text
%        str2double(get(hObject,'String')) returns contents of edtEjectVoltage as a double


% --- Executes during object creation, after setting all properties.
function edtEjectVoltage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtEjectVoltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


