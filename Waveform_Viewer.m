%% Initilization of GUI

%% Waveform_Viewer
function varargout = Waveform_Viewer(varargin)
% WAVEFORM_VIEWER M-file for Waveform_Viewer.fig
%      WAVEFORM_VIEWER, by itself, creates a new WAVEFORM_VIEWER or raises the existing
%      singleton*.
%
%      H = WAVEFORM_VIEWER returns the handle to a new WAVEFORM_VIEWER or the handle to
%      the existing singleton*.
%
%      WAVEFORM_VIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WAVEFORM_VIEWER.M with the given input arguments.
%
%      WAVEFORM_VIEWER('Property','Value',...) creates a new WAVEFORM_VIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Waveform_Viewer_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Waveform_Viewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Waveform_Viewer

% Last Modified by GUIDE v2.5 12-Mar-2010 21:31:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Waveform_Viewer_OpeningFcn, ...
    'gui_OutputFcn',  @Waveform_Viewer_OutputFcn, ...
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


%% Waveform_Viewer_OpeningFcn
% --- Executes just before Waveform_Viewer is made visible.

function Waveform_Viewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Waveform_Viewer (see VARARGIN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------handles Structure
% data_type: which mode the program is in, either 'binary_ternary' or 'DNA_alone', this is a
%    a string and is set in the waveform_settings gui.
%
% zoom: Whether the 'zoom' feature is on or off
%
% trace_data: A cell array containing all of the raw events, since its a cell array each 'row'
%    has a unique number of 'columns'. You access an event for by using {}, and a set of ()
%    after to get a specific data points within that event. The first 'event' in trace_data is the
%    time vector for all of the events, its length is the length of the longest event. So for example,
%    if I wanted the first ten points in the 3rd event I would write, handles.trace_data{4}(1:10),
%    since the first 'event' is the time vector.
%
% time: same thing as trace_data{1}
%
% time_tick: the sample period
%
% num_sweeps: the number of events in the current data file
%
%  edt_exp_descrip: a string describing the current data files experiment, set in the waveform_settings
%
%  edt_exp_voltage: the voltage the experiment was run at (mV)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Choose default command line output for Waveform_Viewer
handles.MATfilename = [];
handles.MATpathname = [];
handles.cursweep = [];
handles.event_string = {};
handles.output = hObject;
handles.zoom = 0;
handles.processed_data = [];
handles.terminal_traces = [];
handles.filename_arc = [];
handles.sweep_data = [];
handles.processed_sweep = [];
handles.baseline_window = .75;
handles.alpha = 0.95;
handles.thresh_amp = -6.5;
handles.data_type = 'binary_ternary';
handles.filename = '';
handles.pathname = '';
handles.defaultPathname = '';
handles.orgfilename = '';
handles.orgpathname = '';
handles.baseline_type = ' ms';
handles.rem_trans = 0;
handles.curindex = 0;
handles.firstLoad = 1;
handles.plotWidthMin = 0;
handles.plotWidthMax = 1;
handles.termStep = 0;
handles.termStepLine = 0;
handles.trace_length = 0;
handles.event_trans = 0;
handles.edt_save_filename_value = 'Summary';
handles.plotCHeightMin = 0;
handles.plotCHeightMax = 1;
handles.buttonUp = [0, 0];
handles.buttonDown = [0, 0];
handles.trace_data = [];
handles.voltage_data = [];
handles.numOpenMATFiles = 0;
handles.min_current = 0;
handles.max_current = 50;
handles.preview_status = 'off';
handles.expCoefs = {};
handles.expCoefs_computed = {};
handles.maxTransientLength = {}; %take full length of event, up to this value in seconds.
handles.saturatedSampleTime = {};  %number of samples from end of event to when amp is out of saturation
%zero out open files listbox
set( handles.lstOpenFiles, 'String', {} );
handles.miktexdir = 'C:\Program Files\MiKTeX 2.7\miktex\bin';
handles.gsdir = 'C:\Program Files\gs\gs8.61\bin';
handles.onFigure = 0;
handles.transAddSamples = 0;
handles.trans_remove = 0;
handles.get_raw_data = 0;
handles.combined_ramps = [];

%set position of progress bar to right above main window
tempUnits = get(handles.figure1, 'Units');
set(handles.figure1, 'Units', 'pixels');
currentMainWindowPosition = get(handles.figure1,'Position');
screenSize = get(0, 'ScreenSize');
%normalize to be less than 1
handles.progressBarPos = [currentMainWindowPosition(1)/screenSize(3),(currentMainWindowPosition(2)+currentMainWindowPosition(4)+55)/screenSize(4)];
set(handles.figure1, 'Units', tempUnits);
set(handles.edt_sample_freq, 'String', {' '})

%load the Waveform Settings GUI
handles.settings = Waveform_Settings('Visible', 'off','handles',handles);

% Update handles structure
set(handles.figure1,'KeyPressFcn',@keyboardShortcuts);
set(handles.listbox_sweeps,'KeyPressFcn',@keyboardShortcuts);
set(handles.edtWinLen,'KeyPressFcn',@checkForCharacters);
set(handles.edtVar,'KeyPressFcn',@checkForCharacters);
set(handles.edtPVal,'KeyPressFcn',@checkForCharacters);
set(handles.edt_beginning_trim,'KeyPressFcn',@checkForCharacters);
set(handles.edt_end_trim,'KeyPressFcn',@checkForCharacters);
set(handles.edt_exp_voltage,'KeyPressFcn',@checkForCharacters);
guidata(hObject, handles);


% UIWAIT makes Waveform_Viewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


%% Waveform_Viewer_OutputFcn
% --- Outputs from this function are returned to the command line.

function varargout = Waveform_Viewer_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function Exit_GUI_Callback(hObject, eventdata, handles)
% hObject    handle to Exit_GUI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.settings)
delete(handles.figure1)

function listbox_sweeps_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_sweeps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox_sweeps contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_sweeps

switch handles.preview_status
    %If we aren't looking at the event plot
    case 'off'
        %get the index of the highlighted event
        handles.curindex = get(handles.listbox_sweeps, 'Value');

        %Makes the plot work ;)
        handles.firstLoad = 1;

        %Plot the event
        load_current_sweep(hObject, handles);

        %update from load_current_sweep
        handles = guidata(hObject);
        guidata(hObject, handles);
end

function btn_chk_term_Callback(hObject, eventdata, handles)
% hObject    handle to btn_chk_term (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch handles.data_type
    case 'binary_ternary'
        % Get currently highlighted event
        handles.curindex = get(handles.listbox_sweeps, 'Value');

        % Old terminal step detection
        % processedSweep = analyze_term(hObject,handles,'yes',0);

        % New One
        [processedSweep, processedSteps] = ChisquaredFit_Custom(hObject, handles, 'yes');

        handles = guidata(hObject);

        % Mark it as processed
        handles.processed_sweeps{handles.matIndex}(handles.event) = 1;

        % Store the steps multi/single
        handles.sweep_data{handles.matIndex}{handles.event,1} = processedSweep;
        handles.sweep_data{handles.matIndex}{handles.event,2} = processedSteps;

        %make the box red or green if all of the events have been analyzed
        if nnz(1-cell2mat(handles.processed_sweeps')) > 0
            set(handles.edt_events_analyzed, 'BackgroundColor', [1 0 0]);
        else
            set(handles.edt_events_analyzed, 'BackgroundColor', [0 1 0]);
        end

        % So you can keep on using the arrows to go through the events
        uicontrol(handles.listbox_sweeps);
        guidata(hObject,handles);

    case 'ramping'
        handles.curindex = get(handles.listbox_sweeps, 'Value');

        processedSweep = ChisquaredFit_Custom(hObject, handles, 'yes');

        handles.processed_sweeps{handles.matIndex}(handles.event) = 1;
        handles.sweep_data{handles.matIndex}{handles.event,1} = processedSweep;

        %make the box red or green if all of the events have been analyzed
        if nnz(1-cell2mat(handles.processed_sweeps')) > 0
            set(handles.edt_events_analyzed, 'BackgroundColor', [1 0 0]);
        else
            set(handles.edt_events_analyzed, 'BackgroundColor', [0 1 0]);
        end

        uicontrol(handles.listbox_sweeps);
        guidata(hObject,handles);

        [trace_data, time, voltage] = get_event(handles, hObject);

        try
            last10ms_min = min(voltage(end-round(5e-3/handles.cur_time_tick):end-round(.1e-3/handles.cur_time_tick)));
            kink_point = find(voltage(1:end-round(6e-3/handles.cur_time_tick)) < last10ms_min,1,'last');
        catch
            last10ms_min = mean(voltage);
            kink_point = find(voltage < last10ms_min, 1, 'last');
        end

        difference = trace_data./voltage;

        figure(10);
        subplot(3,2,1); plot(trace_data); subplot(3,2,2); plot(filter(handles.b,handles.a,trace_data));
        %                      subplot(3,2,3); plot(handles.combined_ramps); subplot(3,2,4); plot(filter(b,a,handles.combined_ramps));
        subplot(3,2,3); cla; plot(voltage); hold on; plot([kink_point kink_point], ylim, '-r'); subplot(3,2,4); plot(filter(handles.b,handles.a,voltage));
        subplot(3,2,5); plot(difference, 'b'); hold on; plot([round(processedSweep(11)/(handles.SampleInt*1e-6)) round(processedSweep(11)/(handles.SampleInt*1e-6))],  [processedSweep(2) processedSweep(2)] ,'-r');subplot(3,2,6); plot(filter(handles.b,handles.a,difference),'b')
    case 'primer'

        handles.curindex = get(handles.listbox_sweeps, 'Value');

        processedSweep = ChisquaredFit_Custom(hObject, handles, 'yes');

        handles.processed_sweeps{handles.matIndex}(handles.event) = 1;
        handles.sweep_data{handles.matIndex}{handles.event,1} = processedSweep;

        %make the box red or green if all of the events have been analyzed
        if nnz(1-cell2mat(handles.processed_sweeps')) > 0
            set(handles.edt_events_analyzed, 'BackgroundColor', [1 0 0]);
        else
            set(handles.edt_events_analyzed, 'BackgroundColor', [0 1 0]);
        end

        uicontrol(handles.listbox_sweeps);
        guidata(hObject,handles);
end

function btn_reset_zoom_Callback(hObject, eventdata, handles)
% hObject    handle to btn_reset_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Just plot like you just clicked the same one

hObject = handles.listbox_sweeps;
listbox_sweeps_Callback(hObject,[],handles);
handles = guidata(hObject);
guidata(hObject,handles);

function btn_proc_all_Callback(hObject, eventdata, handles)
% hObject    handle to btn_proc_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Process all of the events, just will 'click' each one individually
%calculate progress bar location above main window
tempUnits = get(handles.figure1, 'Units');
set(handles.figure1, 'Units', 'pixels');
currentMainWindowPosition = get(handles.figure1,'Position');
screenSize = get(0, 'ScreenSize');
%normalize to be less than 1
handles.progressBarPos = [currentMainWindowPosition(1)/screenSize(3),(currentMainWindowPosition(2)+currentMainWindowPosition(4)+55)/screenSize(4)];
set(handles.figure1, 'Units', tempUnits);

%get the first handles.matIndex & handles.event
handles.curindex = 1;
get_event(handles, hObject);
handles = guidata(hObject);
tic
for i=1:handles.num_events
    stopBar = progressbar(i/handles.num_events, handles.progressBarPos);
    if(stopBar)
        return
    end

    %grab the matIndex and event, trick from Noah ;)
    handles.curindex = i;
    tempEventString = handles.event_string{handles.curindex};
    selectedSweepInfo = sscanf(tempEventString, '%d.%d');
    handles.matIndex = selectedSweepInfo(1);
    handles.event = selectedSweepInfo(2);

    %if we've already processed it
    if (handles.processed_sweeps{handles.matIndex}(handles.event) == 0)
        try
            %                     processedSweep = analyze_term(hObject,handles,index,'no',0);
            [processedSweep, processedSteps] = ChisquaredFit_Custom(hObject, handles, 'no');
            handles = guidata(hObject);
            %                     terminal_step(ind,:) = processedSweep;

            % Save that the event has been processed
            handles.processed_sweeps{handles.matIndex}(handles.event) = 1;

            %For single step detection and general info
            handles.sweep_data{handles.matIndex}{handles.event,1} = processedSweep;

            %For multistep detection
            handles.sweep_data{handles.matIndex}{handles.event,2} = processedSteps;

            if( isempty(processedSteps) )
                %                 disp(['empty @: ' num2str(i)])
                %                 return
            end
        catch
            disp(['analyze term errored on waveform #: ' tempEventString]);
            disp(lasterror)
        end
    end

end
toc
beep;
set(handles.edt_events_analyzed, 'BackgroundColor', [0 1 0]);
guidata(hObject,handles);

function btn_clear_analyzed_Callback(hObject, eventdata, handles)
% hObject    handle to btn_clear_analyzed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%

% Clear all processed events

for i=1:length(handles.processed_sweeps)
    handles.processed_sweeps{i}(:) = 0;
    handles.sweep_data{i}(:,:) = [];
end

set(handles.edt_events_analyzed, 'BackgroundColor', [1 0 0]);

guidata(hObject,handles);


function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(handles.settings);
delete(hObject);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function checkForCharacters(src,evnt)
%keyPressFcn automatically takes in two inputs
%src is the object that was active when the keypress occurred
%evnt stores the data for the key pressed
%to get which box we're in get(src,'Tag')
handles = guidata(src);
set(src,'Enable','inactive')

k = evnt.Key;

% if ( (regexp(k, '[\d\b]') == 1) || strcmp(k,'backspace'))
if ( ~isempty(str2num(k)) || strcmp(k,'backspace') || strcmp(k,'period'))
    set(src,'Enable', 'on')
end

guidata(src,handles)
pause(.01)
set(src,'Enable','on')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function keyboardShortcuts(src,evnt)
%keyPressFcn automatically takes in two inputs
%src is the object that was active when the keypress occurred
%evnt stores the data for the key pressed

%brings in the handles structure in to the function
control = 0;
alt = 0;
shift =0;

%determine which modifiers have been pressed
for x=1:length(evnt.Modifier)
    switch(evnt.Modifier{x})
        case 'control'
            control = 1;
        case 'alt'
            alt = 1;
        case 'shift'
            shift = 1;
    end
end

handles = guidata(src);

k= evnt.Key; %k is the key that is pressed

if( strcmp( k, 's'))
    pause(0.01)
    handles.detected_events(get(handles.listbox_sweeps,'Value'),:)
end

if( strcmp( k, 'c' ) )
    pause(0.01) %allows time to update

    %define hObject as the object of the callback that we are going to use
    %in this case, we are mapping the enter key to the add_pushbutton
    %therefore, we define hObject as the add pushbutton
    %this is done mostly as an error precaution
    hObject = handles.btn_chk_term;

    %call the add pushbutton callback.
    %the middle argument is not used for this callback
    btn_chk_term_Callback(hObject, [], handles);
end

function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%For zooming in on the x and y axis
curPos = get( handles.sweep_plot, 'CurrentPoint' );

%if vertical current resize
if( (curPos(1,1) < handles.plotWidthMin) && (curPos(1,1) > (handles.plotWidthMin - (.07*abs(handles.plotWidthMax-handles.plotWidthMin))))&& (curPos(1,2) > handles.plotCHeightMin) && (curPos(1,2) < handles.plotCHeightMax) )
    set( handles.figure1, 'Pointer', 'top' );
    handles.onFigure = 1;

    %draw line denoting resize region
    handles.resizeLine1 = line( [handles.plotWidthMin handles.plotWidthMax], [curPos(1,2) curPos(1,2)], 'LineStyle', '--' );
    handles.resizeLine2 = line( [handles.plotWidthMin handles.plotWidthMax], [curPos(1,2) curPos(1,2)], 'LineStyle', '--'  );

    %if horizontal current resize
elseif( (curPos(1,2) < handles.plotCHeightMin) && (curPos(1,2) > (handles.plotCHeightMin - (.1*abs(handles.plotCHeightMax-handles.plotCHeightMin)))) && (curPos(1,1) > handles.plotWidthMin) && (curPos(1,1) < handles.plotWidthMax) )
    set( handles.figure1, 'Pointer', 'left' );
    handles.onFigure = 2;

    %draw line denoting resize region
    handles.resizeLine1 = line( [curPos(1,1) curPos(1,1)], [handles.plotCHeightMin handles.plotCHeightMax], 'LineStyle', '--'  );
    handles.resizeLine2 = line( [curPos(1,1) curPos(1,1)], [handles.plotCHeightMin handles.plotCHeightMax], 'LineStyle', '--'  );

else
    set( handles.figure1, 'Pointer', 'arrow' );
    handles.onFigure = 0;
end

handles.buttonDown = [ curPos(1,1) curPos(1,2) ];

% Update handles structure
guidata(hObject, handles);

% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%For the zoom in, the makes the blue lines move

curPos = get( handles.sweep_plot, 'CurrentPoint' );

if( (curPos(1,1) < handles.plotWidthMin) && (curPos(1,1) > (handles.plotWidthMin - (.07*abs(handles.plotWidthMax-handles.plotWidthMin))))&& (curPos(1,2) > handles.plotCHeightMin) && (curPos(1,2) < handles.plotCHeightMax) )
    set( handles.figure1, 'Pointer', 'top' );
    %move line to where mouse button is released
    if( handles.onFigure == 1)
        set( handles.resizeLine2, 'YData', [curPos(1,2), curPos(1,2)] );
    end

    %if horizontal current resize
elseif( (curPos(1,2) < handles.plotCHeightMin) && (curPos(1,2) > (handles.plotCHeightMin - (.1*abs(handles.plotCHeightMax-handles.plotCHeightMin)))) && (curPos(1,1) > handles.plotWidthMin) && (curPos(1,1) < handles.plotWidthMax) )
    set( handles.figure1, 'Pointer', 'left' );
    %move line to where mouse button is released
    if( handles.onFigure == 2 )
        set( handles.resizeLine2, 'XData', [curPos(1,1), curPos(1,1)] );
    end
    %if not in any hot zone
else
    set( handles.figure1, 'Pointer', 'arrow' );
end


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% For zoom in, once the user lets up on the mouse button then we resize the
% graph
curPos = get( handles.sweep_plot, 'CurrentPoint' );

handles.buttonUp = [ curPos(1,1) curPos(1,2) ];

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

%if vertical current resize
if( (handles.onFigure == 1 ) && (curPos(1,1) < handles.plotWidthMin) && (curPos(1,1) > (handles.plotWidthMin - (.07*abs(handles.plotWidthMax-handles.plotWidthMin))))&& (curPos(1,2) > handles.plotCHeightMin) && (curPos(1,2) < handles.plotCHeightMax) )
    set( handles.figure1, 'Pointer', 'top' );

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
    if( handles.firstLoad == 0)
        load_current_sweep( hObject, handles );
    end
    %get changes to handles made in loadCurrentSweep
    handles = guidata( hObject );

    %if horizontal current resize
elseif( (handles.onFigure == 2) && (curPos(1,2) < handles.plotCHeightMin) && (curPos(1,2) > (handles.plotCHeightMin - (.1*abs(handles.plotCHeightMax-handles.plotCHeightMin)))) && (curPos(1,1) > handles.plotWidthMin) && (curPos(1,1) < handles.plotWidthMax) )
    set( handles.figure1, 'Pointer', 'left' );

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
    if( handles.firstLoad == 0)
        load_current_sweep( hObject, handles );
    end
    %get changes to handles made in loadCurrentSweep
    handles = guidata( hObject );

    %not a resize event
else
    set( handles.figure1, 'Pointer', 'arrow' );

    try
        %remove cursor line
        delete(handles.resizeLine1);
        delete(handles.resizeLine2);
    catch
        %ignore if doesn't exist
    end


end

handles.onFigure = 0;

% Update handles structure
uicontrol(handles.listbox_sweeps);
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load current sweep
function load_current_sweep(hObject,handles)

handles.plotDelta = 1;
set(handles.txt_loading, 'Visible', 'on');

% set the current axis to our plot
axes(handles.sweep_plot)

pause(0.001);

%get event data
[data, time, voltage] = get_event(handles, hObject);

handles = guidata(hObject);

handles.timeVector = time;

if isempty(voltage)
    handles.numSignals = 1;
else
    handles.numSignals = 2;
end

% This will be the default behavior unless we are zooming
if( handles.firstLoad == 1 )

    %set dimensions of plot windows
    handles.plotWidthMinDefault = 0;
    handles.plotWidthMaxDefault = length(data)*handles.cur_time_tick;

    if( handles.primer_file{handles.matIndex} == 1 && handles.detected_events{handles.matIndex}(handles.event,8) ~= 0 )
        try
            %remove cursor line
            delete(handles.primer_transient_removal_line);
        catch
        end
        handles.plotCHeightMaxDefault = handles.holding_voltage{handles.matIndex}+10;
        handles.plotCHeightMinDefault = handles.eject_voltage{handles.matIndex}+20;
    else
        handles.plotCHeightMinDefault = -20;
        handles.plotCHeightMaxDefault = 40;
    end

    plotStartEndBuffer = 0.025;

    axis([time(1)-(time(end).*plotStartEndBuffer),time(end)+(time(end).*plotStartEndBuffer),min(data)-10,max(data)+10]);

    %set dimensions of plot windows
    xaxis = get(handles.sweep_plot,'XLim');
    yaxis = get(handles.sweep_plot,'YLim');
    handles.plotWidthMin = xaxis(1);
    handles.plotWidthMax = xaxis(2);
    handles.plotVHeightMin = -250;
    handles.plotVHeightMax = 250;
    handles.trace_length = time(end);

    if( handles.primer_file{handles.matIndex} == 1 && handles.detected_events{handles.matIndex}(handles.event,8) ~= 0 )
        handles.plotCHeightMax = handles.plotCHeightMaxDefault;
        handles.plotCHeightMin = handles.plotCHeightMinDefault;
    else
        handles.plotCHeightMin = yaxis(1);
        handles.plotCHeightMax = yaxis(2);
    end
    if( handles.numSignals == 2 )
        set(handles.txt_voltage, 'String', ['Event Voltage:' num2str(mean(voltage))]);
        set(handles.edt_exp_voltage, 'String', num2str(mean(voltage)));
    end

   

    handles.firstLoad = 0;

end

% % For transient removal
% if( get(handles.chkbox_trans_remove_enable, 'Value') == 1)
%     if(handles.maxTransientLength{handles.matIndex}(handles.detected_events{handles.matIndex}(handles.event,7)+1) == 0)
%         set(handles.edt_max_trans_length, 'String', '10')
%         set(handles.edt_saturation_sample_time, 'String', '0.8')
%     else
%         set(handles.edt_max_trans_length, 'String', num2str(handles.maxTransientLength{handles.matIndex}(handles.detected_events{handles.matIndex}(handles.event,7)+1)*1e3));
%         set(handles.edt_saturation_sample_time, 'String', num2str(handles.saturatedSampleTime{handles.matIndex}(handles.detected_events{handles.matIndex}(handles.event,7)+1)*double(handles.time_tick{handles.matIndex})*1e3));
%     end
% end

% %Fishing file specific
% if(handles.fishing_file{handles.matIndex} == 1)
%     set(handles.txt_molecule_number, 'String', num2str(handles.detected_events{handles.matIndex}(handles.event,7)))
% end
% 
% % Ramping file specific
% if(handles.ramping_file{handles.matIndex} == 1)
%     if(handles.openChanRamps{handles.matIndex}(handles.event) == 1)
%         set(handles.txt_ramp_is_unbound, 'BackgroundColor', [0 1 0]);
%     else
%         set(handles.txt_ramp_is_unbound, 'BackgroundColor', [1 0 0]);
%     end
% end

%this is a speedup borrowed from Noah's code when we are plotting large
%data sets

%adjust plotDelta based on how much the signal is zoomed
handles.plotDelta = round(((handles.plotWidthMax - handles.plotWidthMin)/handles.cur_time_tick)/70000);
if(handles.plotDelta <= 0) handles.plotDelta = 1; end
handles.plotIndexBegin = round(handles.plotWidthMin / handles.cur_time_tick);
handles.plotIndexEnd = round(handles.plotWidthMax / handles.cur_time_tick);

if( handles.plotIndexBegin <= 0) handles.plotIndexBegin = 1; end
if( isequal(handles.plotIndexEnd, 0) )
    handles.plotIndexEnd = 1;
elseif( handles.plotIndexEnd > length(data) )
    handles.plotIndexEnd = length(data);
end
plot(time(handles.plotIndexBegin:handles.plotDelta:handles.plotIndexEnd), data(handles.plotIndexBegin:handles.plotDelta:handles.plotIndexEnd,1), 'b-');
axis([handles.plotWidthMin, handles.plotWidthMax, handles.plotCHeightMin, handles.plotCHeightMax]);

if ( handles.processed_sweeps{handles.matIndex}(handles.event) == 1 && handles.primer_file{handles.matIndex} ~= 1 )
    if isempty(handles.sweep_data{handles.matIndex}{handles.event,2})
        handles.termStep = handles.sweep_data{handles.matIndex}{handles.event,1}(11);
        if handles.sweep_data{handles.matIndex}{handles.event,1}(1) ~= 0 && handles.termStep < handles.plotWidthMax && handles.termStep > handles.plotWidthMin
            xdata = [handles.termStep handles.termStep];
            hold on;
            handles.termStepLine = line(xdata, [handles.plotCHeightMin handles.plotCHeightMax],'color','r','linestyle','-');
            hold off;
        end
    else
        %         colormap_to_use = colormap('lines');
        %         colormap_to_use(7,:) = [];
        for i=1:size(handles.sweep_data{handles.matIndex}{handles.event,2},1)
            %JSJ: where the plotting happens! This data needs the first two
            %columns of processedSteps... maybe I can get away with using
            %the next two columns...
            if i==1
                hold on; plot([0 handles.sweep_data{handles.matIndex}{handles.event,2}(i,2)], [handles.sweep_data{handles.matIndex}{handles.event,2}(i,1) handles.sweep_data{handles.matIndex}{handles.event,2}(i,1)], 'LineStyle', '-', 'Color', 'r');
            else
                plot([sum(handles.sweep_data{handles.matIndex}{handles.event,2}(1:i-1,2)) sum(handles.sweep_data{handles.matIndex}{handles.event,2}(1:i,2))], [handles.sweep_data{handles.matIndex}{handles.event,2}(i,1) handles.sweep_data{handles.matIndex}{handles.event,2}(i,1)], 'LineStyle', '-', 'Color', 'r');
            end

        end
        hold off;
    end
end

%If we are a primer file and the event is a primer detection, then display
%where we changed voltages
if( handles.primer_file{handles.matIndex} == 1 && handles.detected_events{handles.matIndex}(handles.event,8) ~= 0)
    line([time(handles.detected_events_voltages{handles.matIndex}(handles.event,1)) time(handles.detected_events_voltages{handles.matIndex}(handles.event,1))], ylim(), 'color', 'r', 'linestyle', '-')
    line([time(handles.detected_events_voltages{handles.matIndex}(handles.event,2)) time(handles.detected_events_voltages{handles.matIndex}(handles.event,2))], ylim(), 'color', 'k', 'linestyle', '-')
    line([time(handles.detected_events_voltages{handles.matIndex}(handles.event,2)+str2num(get(handles.edt_primer_transient_trim,'String'))) time(handles.detected_events_voltages{handles.matIndex}(handles.event,2)+str2num(get(handles.edt_primer_transient_trim,'String')))], ylim(), 'color', 'k', 'linestyle', '--');
    %     handles.primer_transient_removal_line =
    %     set(handles.primer_transient_removal_line, 'XData', [time(handles.detected_events_voltages{handles.matIndex}(handles.event,2)+str2num(get(handles.edt_primer_transient_trim,'String'))) time(handles.detected_events_voltages{handles.matIndex}(handles.event,2)+str2num(get(handles.edt_primer_transient_trim,'String')))])
    if( handles.detected_events{handles.matIndex}(handles.event,8) == 11 )
        line([time(handles.detected_events_voltages{handles.matIndex}(handles.event,3)) time(handles.detected_events_voltages{handles.matIndex}(handles.event,3))], ylim(), 'color', 'g', 'linestyle', '-');
    end
end
%set the axis labels appropriately
if handles.plotWidthMax < 1e-003
    xlabel('1e-4 s');
else if handles.plotWidthMax < 1e-002
        xlabel('ms');
    else
        xlabel('s');
    end
end
ylabel('pA');

set(handles.txt_loading, 'Visible', 'off');
set(handles.sweep_plot, 'ButtonDownFcn', @sweep_plot_ButtonDownFcn );

guidata(hObject,handles);

% --- Executes on button press in btn_OpenMAT.
function btn_OpenMAT_Callback(hObject, eventdata, handles)
% hObject    handle to btn_OpenMAT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.sweep_plot);
cla;

%show load MAT file window
set(handles.pnlOpenMATFiles, 'Visible', 'on');

% --- Executes on button press in btn_LoadMAT.
function btn_LoadMAT_Callback(hObject, eventdata, handles)
% hObject    handle to btn_LoadMAT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.currentMATFileIndex = length( get(handles.lstOpenFiles, 'String') );

%set default search path for open dialog
if( ~exist(handles.orgpathname)  )
    handles.orgpathname = [pwd '\'];
end

[filename, pathname] = uigetfile('*.mat','Open Event Data', [handles.orgpathname handles.orgfilename],'MultiSelect', 'on');

if isequal(filename,0) || isequal(pathname,0)
    %do nothing, user pressed cancel
    handles.orgpathname = [pwd '\'];
else

    %check if one file selected and if duplicates

    if( iscell( filename) )
        %update file info, append files to array of open files
        if iscell( pathname )

            handles.orgpathname = pathname{1};
        else
            handles.orgpathname = pathname;
        end
        handles.orgfilename = filename{1};
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
        %update file info, append files to array of open files
        handles.orgpathname = pathname;
        handles.orgfilename = filename;
        %only one file selected
        lengthFilename = 1;
    end

    %if any files left after duplicate checking
    if( lengthFilename > 0 )
        for i = 1:lengthFilename
            %catch if only one file selected
            try
                handles.MATfilename{length(handles.MATfilename)+1} = filename{i};
            catch
                handles.MATfilename{length(handles.MATfilename)+1} = filename;
            end
            handles.MATpathname{length(handles.MATpathname)+1} = pathname;

            %get current contents of listbox
            listFiles = get(handles.lstOpenFiles, 'String');

            %add new file to list
            listFiles{length(listFiles)+1} = handles.MATfilename{end};

            %update list
            set(handles.lstOpenFiles, 'String', sort(listFiles) );
            set(handles.lstOpenFiles, 'Value', 1 );

            handles.numOpenMATFiles = length(listFiles);
        end
    end

end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in btn_RemoveMAT.
function btn_RemoveMAT_Callback(hObject, eventdata, handles)
% hObject    handle to btn_RemoveMAT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

listFiles = get(handles.lstOpenFiles, 'String');
indexFile = get(handles.lstOpenFiles, 'Value');

if( indexFile <= length( listFiles ) )
    selectedFile = listFiles{indexFile};

    %convert from list index to handles.filename index
    listIndex = indexFile;
    indexFile = find( strcmp( selectedFile, handles.MATfilename ) == 1 );

    %loop over files after removed file and shift items down
    %to fill in gap
    %if first file selected
    if( indexFile == 1 )
        handles.MATfilename = {handles.MATfilename{2:end}};
        handles.MATpathname = {handles.MATpathname{2:end}};
        handles.processed_sweeps = {handles.processed_sweeps{2:end}};
        handles.sweep_data = {handles.sweep_data{2:end}};
        %if last file is selected
    elseif( indexFile == length( listFiles ) )
        handles.MATfilename = {handles.MATfilename{1:end-1}};
        handles.MATpathname = {handles.MATpathname{1:end-1}};
        handles.processed_sweeps = {handles.processed_sweeps{1:end-1}};
        handles.sweep_data = {handles.sweep_data{1:end-1}};
    else
        %if any other file selected
        handles.MATfilename = {handles.MATfilename{1:indexFile-1},handles.MATfilename{indexFile+1:end}};
        handles.MATpathname = {handles.MATpathname{1:indexFile-1},handles.MATpathname{indexFile+1:end}};
        handles.processed_sweeps = {handles.processed_sweeps{1:indexFile-1},handles.processed_sweeps{indexFile+1:end}};
        handles.sweep_data = {handles.sweep_data{1:indexFile-1},handles.sweep_data{indexFile+1:end}};
    end

    if( indexFile > 1 )
        handles.currentFileIndex = indexFile-1;
    elseif( indexFile == 0 )
        handles.currentFileIndex = 0;
    else
        handles.currentFileIndex = 1;
    end

    %remove element from list
    listFiles(listIndex) = [];

    handles.numOpenMATFiles = length(listFiles);
    if( handles.numOpenMATFiles == 0 )
        %zero out cell arrays if nothing left
        disp('cleared')
        handles.MATfilename = {};
        handles.MATpathname = {};
        handles.processed_sweeps = {};
        handles.sweep_data = {};
        handles.cursweep = [];
        handles.filename = [];
        handles.pathname = [];
        handles.expCoefs = [];
        handles.expCoefs_computed = [];
        handles.allFilenames = {};
        handles.allPathnames = {};
        handles.num_events = 0;
        handles.voltage_data = [];
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

end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in btn_CloseLoadMAT.
function btn_CloseLoadMAT_Callback(hObject, eventdata, handles)
% hObject    handle to btn_CloseLoadMAT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%if files were selected, update eventDetector window
set( handles.btn_CloseLoadMAT, 'Enable', 'off' );

%Assume the files selected are a mix of the old and the new file types
%Lets construct our event string so that it takes the form x.y, where x
%is the mat file number and y is the event within that mat file.
%handles.detected_events will then be a 3 dimensional arrary where the 3 dim
%is the mat number and the first 2 dims will be the detected events
%matrix.  We will make allPathnames & allFilenames a cell array of cell
%arrays. To access the correct path and filenames for an event we will
%take handles.allFilenames{x}(handles.detected_events{x}(y,10)). The same for handles.detected_events,
%except handles.detected_events{x}(y,:)

%What will show in the listbox

if handles.numOpenMATFiles > 0
    %load the files
    handles.event_string = {};
    handles.processed_data = [];
    handles.processed_sweeps = {};
    handles.sweep_data = {};
    handles.detected_events = {};
    handles.detected_events_ms = {};
    handles.detected_events_voltages = {};
    handles.allFilenames = {};
    handles.allPathnames = {};
    handles.num_events = 0;
    handles.num_molecules = [];
    handles.trace_data = [];

    for i=1:handles.numOpenMATFiles

        %Load each of the MAT files
        load([handles.MATpathname{i} handles.MATfilename{i}]);

        %handles.allFilenames/allPathnames holds all the pathnames and
        %filesnames
        if iscell(poreEventData.filename)
            handles.allFilenames{i} = poreEventData.filename;
            handles.allPathnames{i} = poreEventData.pathname;
        else
            handles.allFilenames{i} = {poreEventData.filename};
            poreEventData.filename = {poreEventData.filename};
            handles.allPathnames{i} = {poreEventData.pathname};
            poreEventData.pathname = {poreEventData.pathname};
        end

        handles.time_tick{i} = poreEventData.samplePeriod;

        handles.time_tick_options{i} = cellstr([num2str(1./(poreEventData.samplePeriod*1e3.*[1:5])',4), repmat(' kHz', 5, 1)]);

        %detected events holds the event information and ptr to the
        %allFilenames
        %add current file index to detectedEvents arrays
        %if mat file is old might not have all the columns. just
        %add the proper number of columns to make it the right
        %width and then put filenumber on the end
        %width is 11 for handles.detectedEvents
        %and 12 for handles.detectedEvents_ms
        [m,n] = size(poreEventData.detectedEvents);

        if( n < 11 )
            handles.detected_events{i} = [poreEventData.detectedEvents, zeros(m,9-n), ones(m,1), zeros(m,1)];
        else
            handles.detected_events{i} = poreEventData.detectedEvents;
        end
        handles.detected_events_ms{i} = poreEventData.detectedEvents_ms;
        %handles.num_events contains the number of events
        temp_num_events = size(handles.detected_events{i},1);
        handles.num_events = handles.num_events + temp_num_events;
        handles.ramping_file{i} = 0;
        handles.primer_file{i} = 0;
        
%         handles.fishing_file{i} = poreEventData.fishingFile;
% 
%         try
%             handles.ramping_file{i} = poreEventData.rampingFile;
%         catch
%             handles.ramping_file{i} = 0;
%         end
% 
%         try
%             handles.primer_file{i} = poreEventData.primerFile;
%             handles.detected_events_voltages{i} = poreEventData.detectedEventsVoltages;
%             handles.holding_voltage{i} = poreEventData.holdingVoltage;
%             handles.probing_voltage{i} = poreEventData.probingVoltage;
%             handles.eject_voltage{i} = poreEventData.ejectVoltage;
%         catch
%             handles.primer_file{i} = 0;
%         end
% 
%         % For ramping file we have some different inputs that show up
%         if(handles.ramping_file{i} == 1)
%             set(handles.pnl_trans_ramp, 'Title', 'Ramp Functions');
%             set(handles.chkbox_trans_remove_enable, 'Visible', 'off');
%             set(handles.text49, 'Visible', 'off');
%             set(handles.text50, 'Visible', 'off');
%             set(handles.edt_saturation_sample_time, 'Visible', 'off');
%             set(handles.edt_max_trans_length, 'Visible', 'off');
%             set(handles.txt_have_coeff, 'Visible', 'off');
%             set(handles.btn_get_coeff, 'Visible', 'off');
%             set(handles.text51, 'Visible', 'off');
%             set(handles.text54, 'Visible', 'off');
%             set(handles.chkbox_analyze_probing_events, 'Visible', 'off');
%             set(handles.pnl_ramp_calc_open, 'Visible', 'on');
%             set(handles.pnl_ramp_save, 'Visible', 'on');
%             handles.data_type = 'ramping';
%         elseif( handles.primer_file{i} == 1 )
%             set(handles.pnl_primer_anneal, 'Visible', 'on')
%             set(handles.pnl_trans_ramp, 'Visible', 'off');
%             set(handles.chkbox_trans_remove_enable, 'Visible', 'off');
%             set(handles.text49, 'Visible', 'off');
%             set(handles.text50, 'Visible', 'off');
%             set(handles.edt_saturation_sample_time, 'Visible', 'off');
%             set(handles.edt_max_trans_length, 'Visible', 'off');
%             set(handles.txt_have_coeff, 'Visible', 'off');
%             set(handles.btn_get_coeff, 'Visible', 'off');
%             set(handles.text51, 'Visible', 'off');
%             set(handles.text54, 'Visible', 'off');
%             set(handles.chkbox_analyze_probing_events, 'Visible', 'off');
%             set(handles.pnl_ramp_calc_open, 'Visible', 'off');
%             set(handles.pnl_ramp_save, 'Visible', 'off');
%             handles.data_type = 'primer';
%         else
%             set(handles.pnl_trans_ramp, 'Title', 'Transient Removal/Fishing');
%             set(handles.pnl_trans_ramp, 'Visible', 'on');
%             set(handles.chkbox_trans_remove_enable, 'Visible', 'on');
%             set(handles.text49, 'Visible', 'on');
%             set(handles.text50, 'Visible', 'on');
%             set(handles.edt_saturation_sample_time, 'Visible', 'on');
%             set(handles.edt_max_trans_length, 'Visible', 'on');
%             set(handles.txt_have_coeff, 'Visible', 'on');
%             set(handles.btn_get_coeff, 'Visible', 'on');
%             set(handles.text51, 'Visible', 'on');
%             set(handles.text54, 'Visible', 'on');
%             set(handles.chkbox_analyze_probing_events, 'Visible', 'on');
%             set(handles.pnl_ramp_calc_open, 'Visible', 'off');
%             set(handles.pnl_ramp_save, 'Visible', 'off');
            handles.data_type = 'binary_ternary';
%         end

        if poreEventData.fishingFile == 1
            handles.num_molecules(i) = max(poreEventData.detectedEvents(:,7))+1;
        end

        %initialize processed_sweeps and sweep_data
        handles.processed_sweeps{i} = zeros(temp_num_events,1);
        handles.sweep_data{i} = cell(temp_num_events,2);
        handles.expCoefs_computed{i} = zeros(max(handles.detected_events{i}(:,7))+1,1);
        handles.maxTransientLength{i} = zeros(max(handles.detected_events{i}(:,7))+1,1);
        handles.saturatedSampleTime{i} = zeros(max(handles.detected_events{i}(:,7))+1,1);
        handles.openChanRamps{i} = zeros(temp_num_events,1);
        %             handles.

        %construct the string cell array that will appear in our list
        %box... don't ask how its contructed
        handles.event_string = [handles.event_string; ...
            cellstr([strjust(num2str(i*ones(temp_num_events,1)),'right'),...
            repmat('.',temp_num_events,1), strjust(num2str((1:temp_num_events)'),'left')])];

        %if the events in the mat file were filtered
        try
            handles.filtered(i) = poreEventData.filtered;
        catch
            handles.filtered(i) = 0;
        end

        %loop over files and locate them all
        if( iscell( poreEventData.filename ) )
            numFiles = length(poreEventData.filename);
        else
            numFiles = 1;
        end

        for currentFile = 1:numFiles
            handles.filename = handles.allFilenames{i}{currentFile};
            handles.pathname = handles.allPathnames{i}{currentFile};

            if( ~exist([handles.pathname handles.filename],'file') )

                correctFile = 0;
                while( correctFile == 0 || correctFile == -1)
                    %load sweep from disk
                    try
                        [Trace_Data, SampleInt] = import_abf( [handles.pathname handles.filename], 1);
                        poreEventData.pathname{currentFile} = handles.pathname;
                        handles.allPathnames{i}{currentFile} = handles.pathname;
                        if correctFile == -1
                            correctFile = 1;
                            save([handles.MATpathname{i} handles.MATfilename{i}], '-append', 'poreEventData');
                        end
                        correctFile = 1;
                    catch
                        if correctFile == -1
                            rethrow(lasterror)
                        end
                        correctFile = -1;

                        disp('file not found');

                        %set default search path for open dialog
                        if( ~exist(handles.defaultPathname)  )
                            handles.defaultPathname = [pwd '\'];
                        end

                        [filename, pathname] = uigetfile('*.abf',['Please find ', handles.filename], [handles.defaultPathname]);

                        if isequal(filename,0) || isequal(pathname,0)
                            %do nothing, user pressed cancel
                            handles.defaultPathname = [pwd '\'];
                            return;
                        else
                            %update file info
                            handles.defaultPathname = pathname;
                            %handles.filename = filename;
                            handles.pathname = pathname;
                        end
                    end
                end %while

            end
        end

    end

    % For filtered data
    if ( handles.filtered == 1 || handles.ramping_file{1} == 1)
        [Trace_Data, SampleInt] = import_abf( [handles.pathname handles.filename], 1);
        cutoff_freq = 5000/((1/(SampleInt*1e-6))/2);
        [z,p,k] = butter(8,cutoff_freq,'low');
        [s,g] = zp2sos(z,p,k);
        handles.Hb = dfilt.df2sos(s,g);

        cur_freq = 1/(SampleInt*1e-6);

        rp = 1;           % Passband ripple in dB
        rs = 50;          % Stopband ripple in Db
        fs = cur_freq;        % Sampling frequency
        f = [1000 3000];    % Cutoff frequencies
        handles.a = [1 0];        % Desired amplitudes
        % Compute deviations
        dev = [(10^(rp/20)-1)/(10^(rp/20)+1)  10^(-rs/20)];
        [n,fo,ao,w] = firpmord(f,handles.a,dev,fs);
        handles.b = firpm(n,fo,ao,w);
    end

    %handles.filename/pathname will hold the current
    %file to load
    %         handles.filename = poreEventData.filename{1};
    %         handles.pathname = poreEventData.pathname{1};
    set(handles.listbox_sweeps,'String',handles.event_string,...
        'Value',1);
    set(handles.pnlOpenMATFiles, 'Visible', 'off');
    set(handles.btn_CloseLoadMAT, 'Enable', 'on' );
    set(handles.figure1, 'Name', ['Waveform_Viewer - ' handles.MATfilename{1}]);
    temptok = regexp(handles.MATfilename{1}, '_(.*)\.mat', 'tokens');
    disp(temptok{:})
    set(handles.edt_exp_descrip, 'String', temptok{1})
    handles.curindex = 1;
    load_current_sweep(hObject, handles);
else
    %clear the lstbox
    handles.event_string = {};
    set(handles.listbox_sweeps,'String',handles.event_string,...
        'Value',1);
    set(handles.pnlOpenMATFiles, 'Visible', 'off');
    set(handles.btn_CloseLoadMAT, 'Enable', 'on' );
    set(handles.figure1, 'Name', 'Waveform_Viewer');
end

handles = guidata(hObject);

%change our indiciator
set(handles.edt_events_analyzed, 'BackgroundColor', [1 0 0]);

guidata(hObject, handles);

function edt_save_filename_Callback(hObject, eventdata, handles)
% hObject    handle to edt_save_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_save_filename as text
%        str2double(get(hObject,'String')) returns contents of edt_save_filename as a double

handles.edt_save_filename_value = get(hObject, 'String');
if( isempty(handles.edt_save_filename_value) )
    handles.edt_save_filename_value = 'Summary';
else
    handles.edt_save_filename_value = strrep(handles.edt_save_filename_value,' ','_');
end
set(handles.edt_save_filename, 'String', handles.edt_save_filename_value);

guidata(hObject, handles);

% --- Executes on selection change in mnu_data_type.
function mnu_data_type_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_data_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns mnu_data_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mnu_data_type

switch get(hObject,'Value')
    case 1
        handles.data_type = 'binary_ternary';
        set(handles.mnu_settings, 'enable', 'off');
    case 2
        handles.data_type = 'Fishing';
        set(handles.mnu_settings, 'enable', 'on');
    case 3
        handles.data_type = 'ramping';
        set(handles.mnu_settings, 'enable', 'off');
end

guidata(hObject, handles);

% --- Executes on button press in btn_preview.
function btn_preview_Callback(hObject, eventdata, handles)
% hObject    handle to btn_preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%This will display the event plot so far
switch handles.preview_status

    case 'off'
        set(hObject, 'String', 'Show Events');
        axes(handles.sweep_plot)
        cla;
        %         keyboard

        % Make a temp copy of the processed data, this will get rid of the
        % cell behavior of the data
        temp_processed_data = cell2mat(handles.sweep_data{:}(:,1));
        temp_processed_data_multi = cell2mat(handles.sweep_data{:}(:,2));

        disp(['hey!: ' num2str(length(temp_processed_data))])
        %find events with ts or not based on temp_processed_data(:,1) == 1 %and no ts with temp_processed_data(:,1) ~= 1
        %if the analzye only probing events is check, than we will get the
        %terminal step for only those types of events
        number_of_steps = str2num(get(handles.edtWinLen,'String'));

        % For typical ts data
        if( ~strcmp(handles.data_type, 'ramping') && number_of_steps == 1 )

            events_with_ts = find (temp_processed_data(:,1) == 1);

            events_wo_ts = find( (temp_processed_data(:,1) <= 0) ); %& (temp_processed_data(:,5) <= handles.max_current) & (temp_processed_data(:,5) >= handles.min_current) );

            number_with_ts = length( events_with_ts );
            number_without_ts = length( events_wo_ts );

            semilogx(temp_processed_data( events_wo_ts ,6).*1000,...
                temp_processed_data(events_wo_ts ,5),'.m','MarkerSize', 10);hold on; grid on;
            %we'll mark each of the molecules with a unique color
            if get(handles.chkbox_analyze_probing_events, 'Value') == 1
                for filenum = 1:handles.numOpenMATFiles;
                    for i=1:handles.num_molecules(filenum)
                        events_with_ts = find (handles.sweep_data{filenum}{:,1}(1) == 1 ...
                            & handles.detected_events{filenum}(:,8) == 3 & handles.detected_events{filenum}(:,7) == i);
                        semilogx(handles.sweep_data{filenum}{events_with_ts,1}(9).*1000,...
                            handles.sweep_data{filenum}{events_with_ts,1}(8),'.','MarkerSize', 10, 'Color', [rand rand rand]);
                        semilogx(temp_processed_data( events_with_ts ,3).*1000, temp_processed_data( events_with_ts ,2),'.c','MarkerSize', 10);
                    end
                end

                semilogx(temp_processed_data( events_with_ts ,3).*1000, temp_processed_data( events_with_ts ,2),'.c','MarkerSize', 10);
            else
                semilogx(temp_processed_data( events_with_ts ,9).*1000,...
                    temp_processed_data( events_with_ts ,8),'.k','MarkerSize', 10);
                semilogx(temp_processed_data( events_with_ts ,3).*1000, temp_processed_data( events_with_ts ,2),'.c','MarkerSize', 10);
            end
            axis([0.2 11000 2 30]);
            xlabel('Duration (ms)','FontSize', 10);
            ylabel('Amplitude (pA)','FontSize', 10);
            legend('\fontsize{10pt}Events "without" terminal steps','\fontsize{10pt}Extracted Pre-Terminal Steps', '\fontsize{10pt}Extracted Terminal Steps');
            handles.preview_status = 'on';

        elseif( strcmp(handles.data_type, 'ramping') )
            % For ramping data
            events_with_ts = find (temp_processed_data(:,1) == 1);

            events_wo_ts = find( (temp_processed_data(:,1) <= 0) ); %& (temp_processed_data(:,5) <= handles.max_current) & (temp_processed_data(:,5) >= handles.min_current) );

            number_with_ts = length( events_with_ts )
            number_without_ts = length( events_wo_ts );
            for filenum = 1:handles.numOpenMATFiles;
                events_with_ts_before_kink = find (handles.sweep_data{filenum}(:,1) == 1 ...
                    & handles.detected_events{filenum}(:,8) == 5 & handles.sweep_data{filenum}(:,12) == 0);
                events_with_ts_after_kink = find (handles.sweep_data{filenum}(:,1) == 1 ...
                    & handles.detected_events{filenum}(:,8) == 5 & handles.sweep_data{filenum}(:,12) == 1);
                semilogx(handles.sweep_data{filenum}( events_with_ts_before_kink ,9).*1000,...
                    handles.sweep_data{filenum}( events_with_ts_before_kink ,8),'.k','MarkerSize', 10);hold on; grid on;
                semilogx(handles.sweep_data{filenum}( events_with_ts_after_kink ,9).*1000,...
                    handles.sweep_data{filenum}( events_with_ts_after_kink ,8),'.r','MarkerSize', 10);
                %                     semilogx(temp_processed_data( events_with_ts ,3).*1000, temp_processed_data( events_with_ts ,2),'.c','MarkerSize', 10);
            end

            %             semilogx(temp_processed_data( events_with_ts ,3).*1000, temp_processed_data( events_with_ts ,2),'.c','MarkerSize', 10);
            axis([0.2 300 -5 30]);
            xlabel('Duration (ms)','FontSize', 10);
            ylabel('Amplitude (pA)','FontSize', 10);
            %             legend('\fontsize{10pt}Events "without" terminal steps','\fontsize{10pt}Extracted Pre-Terminal Steps', '\fontsize{10pt}Extracted Terminal Steps');
            handles.preview_status = 'on';
        else
            % For multiple step

            colormap_to_use = [0 0 0;1 0 .5; 0 0 1; 1 0 1;repmat([0 .95 0], 6,1)];

            num_steps_all = zeros(length(handles.sweep_data{:}(:,2)),1);

            for i=1:length(handles.sweep_data{:}(:,2))
                num_steps_all(i) =size(handles.sweep_data{:}{i,2}(:,2),1)-1;
                if (size(handles.sweep_data{:}{i,2},1) > 1)
                    tempH = semilogx(temp_processed_data(i,6).*1000,...
                        temp_processed_data(i,5),'.b','MarkerSize', 10);hold on; grid on;
                    if ~exist('hSGroup2')
                        disp(['creating hSGroup2'])
                        hSGroup2 = hggroup;
                    end
                    set(tempH, 'Parent',hSGroup2);

                else
                    tempH = semilogx(handles.sweep_data{:}{i,2}(:,2).*1000, handles.sweep_data{:}{i,2}(:,1),'.', 'MarkerSize', 10, 'Color', colormap_to_use(1,:));
                    hold on; grid on;
                    if ~exist('hSGroup1')
                        disp('creating hSGroup1')
                        hSGroup1 =hggroup;
                    end
                    set(tempH, 'Parent', hSGroup1);
                end
            end
            legend([hSGroup2;hSGroup1], 'Events "with" steps', 'Events without steps')

            handles.preview_status = 'on';
        end

    case 'on'

        % Plot the previous event
        axes(handles.sweep_plot)
        ylabel('pA');
        xlabel('');
        legend off;
        grid off;
        hold off;
        cla;
        set(handles.sweep_plot, 'XScale', 'linear');

        load_current_sweep(hObject, handles);
        handles = guidata(hObject);

        set(hObject, 'String', 'Preview');
        handles.preview_status = 'off';

end

guidata(hObject, handles)


% --- Executes on button press in chkbox_adv_options.
function chkbox_adv_options_Callback(hObject, eventdata, handles)
% hObject    handle to chkbox_adv_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkbox_adv_options

if get(hObject, 'Value') == 1
    set(handles.uipanel_points_to_trim, 'Visible', 'on')
    %set(handles.txt_data_type, 'Visible', 'on');
    %set(handles.mnu_data_type, 'Visible', 'on');
    set(handles.edt_sample_freq, 'Enable', 'on');
    %set(handles.chkbox_ignore_eject, 'Visible', 'on');
else
    set(handles.uipanel_points_to_trim, 'Visible', 'off')
    set(handles.txt_data_type, 'Visible', 'off');
    set(handles.mnu_data_type, 'Visible', 'off');
    set(handles.edt_sample_freq, 'Enable', 'off');
    set(handles.chkbox_ignore_eject, 'Visible', 'off');
end

guidata(hObject, handles);


% --- Executes on button press in chkbox_trans_remove_enable.
function chkbox_trans_remove_enable_Callback(hObject, eventdata, handles)
% hObject    handle to chkbox_trans_remove_enable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkbox_trans_remove_enable

if get(hObject,'Value') == 1
    set(handles.btn_get_coeff, 'Enable', 'on');
    set(handles.text49, 'Enable', 'on');
    set(handles.text50, 'Enable', 'on');
    set(handles.edt_saturation_sample_time, 'Enable', 'on');
    set(handles.edt_max_trans_length, 'Enable', 'on');
    set(handles.text54, 'Enable', 'on');
    set(handles.txt_have_coeff, 'Visible', 'on');
    set(handles.btn_get_coeff, 'Enable', 'on');
else
    set(handles.btn_get_coeff, 'Enable', 'off');
    set(handles.text49, 'Enable', 'off');
    set(handles.text50, 'Enable', 'off');
    set(handles.edt_saturation_sample_time, 'Enable', 'off');
    set(handles.edt_max_trans_length, 'Enable', 'off');
    set(handles.text54, 'Enable', 'off');
    set(handles.txt_have_coeff, 'Visible', 'off');
    set(handles.btn_get_coeff, 'Enable', 'off');
end

handles.trans_remove = get(hObject, 'Value');

guidata(hObject, handles);

% --- Executes on button press in btn_get_coeff.
function btn_get_coeff_Callback(hObject, eventdata, handles)
% hObject    handle to btn_get_coeff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%set some parameters
showPlots = 0;
maxIterations = 750;

% FITcoefs = [45.5519; -22.9977; 647.7805];
FITcoefs = [45.5519; -22.9977; 647.7805];
%initial guess for fit
x0 = FITcoefs;

%pre-allocate memory for exponential coefs
expCoefs = zeros(length(handles.detected_events{1}(:,1)), length(x0));

handles.get_raw_data = 1;

[traceData, time, voltage] = get_event(handles, hObject);

handles = guidata(hObject);

handles.get_raw_data = 0;

numEvents = handles.num_events;
typeAll_I = [1:1:numEvents];

t = handles.detected_events_ms{1}(typeAll_I,4);
a = handles.detected_events_ms{1}(typeAll_I,5);

% figure(5)
% semilogx(t,a,'b.')

if( handles.detected_events{handles.matIndex}(handles.event,8) == 3)

    sampleInt = double(handles.cur_time_tick);
    try
        handles.saturatedSampleTime{handles.matIndex}(handles.detected_events{handles.matIndex}(handles.event,7)+1) = round((str2num(get(handles.edt_saturation_sample_time, 'String'))*1e-3)/sampleInt);
    catch
        handles.saturatedSampleTime{handles.matIndex}(handles.detected_events{handles.matIndex}(handles.event,7)+1) = round((.8*1e-3)/sampleInt);
    end

    try
        handles.maxTransientLength{handles.matIndex}(handles.detected_events{handles.matIndex}(handles.event,7)+1) = str2num(get(handles.edt_max_trans_length, 'String'))*1e-3;
    catch
        handles.maxTransientLength{handles.matIndex}(handles.detected_events{handles.matIndex}(handles.event,7)+1) = 10e-3;
    end

    fprintf(' \n  Fit to transient will begin after %d samples (%f msec) \n Max Transient Length is %d samples (%f msec)\n',...
        handles.saturatedSampleTime{handles.matIndex}(handles.detected_events{handles.matIndex}(handles.event,7)+1),handles.saturatedSampleTime{handles.matIndex}(handles.detected_events{handles.matIndex}(handles.event,7)+1)*sampleInt*1e3, round(handles.maxTransientLength{handles.matIndex}(handles.detected_events{handles.matIndex}(handles.event,7)+1)/sampleInt), handles.maxTransientLength{handles.matIndex}(handles.detected_events{handles.matIndex}(handles.event,7)+1)*1e3);

    eventData = traceData(handles.saturatedSampleTime{handles.matIndex}(handles.detected_events{handles.matIndex}(handles.event,7)+1):end);

    time = 0:sampleInt:(length(eventData)-1)*sampleInt;

    % Correction to use probing event for fitting
    transStartIndex = handles.saturatedSampleTime{handles.matIndex}(handles.detected_events{handles.matIndex}(handles.event,7)+1);
    transientEndIndex = length(traceData);
    if( (transientEndIndex-transStartIndex+1) > (handles.maxTransientLength{handles.matIndex}(handles.detected_events{handles.matIndex}(handles.event,7)+1)/sampleInt) )
        transientEndIndex = round(handles.maxTransientLength{handles.matIndex}(handles.detected_events{handles.matIndex}(handles.event,7)+1)/sampleInt) + transStartIndex + 1;
    end

    transientData = traceData( transStartIndex:transientEndIndex, 1 );

    %Make sure it runs long enough to get a good fit 'display', 'iter',
    options = optimset('maxFunEvals', 1e16, 'maxIter', maxIterations);

    %Check out documentation, lsqcurvefit is just an application of lsqnonlin
    transientData = double(transientData);

    %     if( isfield(handles, 'expCoefs') )
    %         expFitCoefs = poreEventData.expCoefs(currentEvent, :);
    %     else
    %only fit to first transient
    %if( firstRun == 1 )

    fprintf(' fit to event %d \n\n',handles.event)
    if( showPlots == 1 )
        figure(12)
        plot(time(1:length(transientData)),transientData)
        hold on
        grid on
    end
    % pause

    [expFitCoefs, resnorm, RESIDUAL, EXITFLAG] = lsqcurvefit( ...
        @(x,t) exponentialFit(x,time(1:length(transientData))), x0, double( time(1:length(transientData)) ), transientData', [], [], options);
    %                     disp(EXITFLAG)

    % First Fit coeffs hardcoded from DNA alone data
    %expFitCoefs = FITcoefs;

    %save fit coefs
    disp(expFitCoefs)
    expCoefs(1,:) = expFitCoefs;

    %firstRun = 0;

    %update inital guess to previous fit
    x0 = expFitCoefs;

    %add inverse transient to compensate event
    fitData = exponentialFit( expFitCoefs,time(1:length(transientData)) );

    compensatedEvent = [-( fitData - fitData(end) )'+eventData(1:length(transientData)); eventData(length(transientData)+1:end)];

    %compensatedEvent = -fitData' + eventData + fitData(end);

    %calc mean
    meanAmplitude = mean( compensatedEvent);

    %%%%%%%%%%%%%%%%%%%%
    %plot results
    if( showPlots == 1 )
        %                     figure(12)
        %                     plot(time(1:length(transientData)),transientData)
        %                     hold on
        %                     grid on
        %                     %plot(eventData,'r')
        %plot(transientData-eventData,'k')


        figure(14)

        plot( time , eventData, 'r'); hold on, grid on
        % plot( [length(eventData):(length(eventData)+length(transientData)-1)], transientData, 'b');
        plot( time(1:length(transientData)), transientData, 'b');
        %plot(b, trans_decimate, 'go', time(1:length(transientData)), transSpline, 'm')
        plot( time(1:length(transientData)), exponentialFit(expFitCoefs,time(1:length(transientData))), 'r');


        figure(15);
        plot( time, compensatedEvent, 'k');hold on, grid on
        % axis([0 ((length(eventData)+length(transientData)-1)) -100 100]);
        axis([0 time(end) 0 60]);


        plot([0, length(eventData)-1], [meanAmplitude, meanAmplitude], 'g--');

        %         fprintf('eventData[%.4f:%.4f]s, transientData[%.4f,%.4f]s\n', ...
        %             handles.detected_events(handles.event,2)*sampleInt, ...
        %             handles.detected_events(handles.event,3)*sampleInt, ...
        %             (handles.detected_events(handles.event,3)+1)*sampleInt, ...
        %             transientEndIndex*sampleInt );

        % disp(size(expFitCoefs))

        hold off
    end
    handles.expCoefs{handles.matIndex}(handles.detected_events{handles.matIndex}(handles.event,7)+1,:) = expFitCoefs';
    handles.expCoefs_computed{handles.matIndex}(handles.detected_events{handles.matIndex}(handles.event,7)+1) = 1;
    set(handles.txt_have_coeff, 'BackgroundColor', [0 1 0]);
end

guidata(hObject, handles)

function btn_Open_Prev_Session_Callback(hObject, eventdata, handles)
% hObject    handle to btn_Open_Prev_Session (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Attempt to open up a previous session, probably not working right now

%set default search path for open dialog
if( ~exist(handles.orgpathname)  )
    handles.orgpathname = [pwd '\'];
end

[filename, pathname] = uigetfile('*.mat','Open Previous Session Data', [handles.orgpathname handles.orgfilename]);

if isequal(filename,0) || isequal(pathname,0)
    %do nothing, user pressed cancel
    handles.orgpathname = [pwd '\'];
else
    if ~isempty(findstr(filename, '_waveformViewerData.mat'))
        disp('found the right file')
        load([pathname filename])
        handles.trace_data = [];
        handles.voltage_data = [];
        handles.cursweep = [];
        handles.MATfilename = waveformViewerData.MATfilename;
        handles.MATpathname = waveformViewerData.MATpathname;
        handles.matIndex = waveformViewerData.matIndex;
        handles.event_string = waveformViewerData.event_string;
        handles.event = waveformViewerData.event;
        handles.curindex = waveformViewerData.curindex;
        handles.data_type = waveformViewerData.data_type;
        handles.filename = waveformViewerData.filename;
        handles.pathname = waveformViewerData.pathname;
        handles.orgfilename = waveformViewerData.orgfilename;
        handles.orgpathname = waveformViewerData.orgpathname;
        handles.defaultPathname = waveformViewerData.defaultPathname;
        handles.baseline_type = waveformViewerData.baseline_type;
        handles.rem_trans = waveformViewerData.rem_trans;
        handles.event_trans = waveformViewerData.event_trans;
        handles.numOpenMATFiles = waveformViewerData.numOpenMATFiles;
        handles.preview_status = 'off';
        handles.trans_remove = waveformViewerData.trans_remove;
        handles.expCoefs = waveformViewerData.expCoefs;
        handles.expCoefs_computed = waveformViewerData.expCoefs_computed;
        handles.maxTransientLength = waveformViewerData.maxTransientLength;
        handles.saturatedSampleTime = waveformViewerData.saturatedSampleTime;
        handles.processed_sweeps = waveformViewerData.processed_sweeps;
        handles.detected_events = waveformViewerData.detected_events;
        handles.detected_events_ms = waveformViewerData.detected_events_ms;
        handles.allFilenames = waveformViewerData.allFilenames;
        handles.allPathnames = waveformViewerData.allPathnames;
        handles.num_events = waveformViewerData.num_events;
        handles.num_molecules = waveformViewerData.num_molecules;
        handles.time_tick = waveformViewerData.time_tick;
        handles.time_tick_options = waveformViewerData.time_tick_options;
        handles.time_tick_options_value = waveformViewerData.time_tick_options_value;
        handles.fishing_file = waveformViewerData.fishing_file;
        try
            handles.ramping_file = waveformViewerData.ramping_file;
        catch
            handles.ramping_file = cell{length(handles.fishing_file),1};
            handles.ramping_file{:} = [0];
        end
        handles.filtered = waveformViewerData.filtered;
        handles.cur_time_tick = waveformViewerData.cur_time_tick;
        handles.SampleInt = waveformViewerData.SampleInt;
        handles.numSignals = waveformViewerData.numSignals;
        handles.edt_save_filename_value = waveformViewerData.edt_save_filename_value;
        handles.sweep_data = waveformViewerData.sweep_data;

        set(handles.edt_beginning_trim, 'String', waveformViewerData.edt_beginning_trim_value);
        set(handles.edt_end_trim, 'String', waveformViewerData.edt_end_trim_value);
        set(handles.edt_save_filename, 'String', handles.edt_save_filename_value);
        set(handles.edt_exp_stn, 'String', waveformViewerData.edt_exp_stn_value);
        set(handles.edt_exp_voltage, 'String', waveformViewerData.edt_exp_voltage_value);
        set(handles.edt_exp_descrip, 'String', waveformViewerData.edt_exp_descrip_value);
        set(handles.chkbox_trans_remove_enable, 'Value', handles.trans_remove);
        set(handles.listbox_sweeps,'String',handles.event_string,...
            'Value',handles.curindex);
        set(handles.edtPVal, 'String', waveformViewerData.MinStepLength_value);
        set(handles.edtVar, 'String', waveformViewerData.MinStepSize_value);
        set(handles.lstOpenFiles, 'String', waveformViewerData.lstOpenFiles_value)
        if waveformViewerData.adv_options_value == 1
            set(handles.uipanel_points_to_trim, 'Visible', 'on')
            set(handles.txt_data_type, 'Visible', 'on');
            set(handles.mnu_data_type, 'Visible', 'on');
            set(handles.edt_sample_freq, 'Enable', 'on');
            set(handles.chkbox_adv_options, 'Value', 1);
        else
            set(handles.uipanel_points_to_trim, 'Visible', 'off')
            set(handles.txt_data_type, 'Visible', 'off');
            set(handles.mnu_data_type, 'Visible', 'off');
            set(handles.edt_sample_freq, 'Enable', 'off');
            set(handles.chkbox_adv_options, 'Value', 0);
        end

        if waveformViewerData.trans_remove_enable_value == 1
            set(handles.btn_get_coeff, 'Enable', 'on');
            set(handles.text49, 'Enable', 'on');
            set(handles.text50, 'Enable', 'on');
            set(handles.edt_saturation_sample_time, 'Enable', 'on');
            set(handles.edt_max_trans_length, 'Enable', 'on');
            set(handles.text54, 'Enable', 'on');
            set(handles.txt_have_coeff, 'Visible', 'on');
            set(handles.btn_get_coeff, 'Enable', 'on');
            set(handles.chkbox_trans_remove_enable, 'Value', 1)
        else
            set(handles.btn_get_coeff, 'Enable', 'off');
            set(handles.text49, 'Enable', 'off');
            set(handles.text50, 'Enable', 'off');
            set(handles.edt_saturation_sample_time, 'Enable', 'off');
            set(handles.edt_max_trans_length, 'Enable', 'off');
            set(handles.text54, 'Enable', 'off');
            set(handles.txt_have_coeff, 'Visible', 'off');
            set(handles.btn_get_coeff, 'Enable', 'off');
            set(handles.chkbox_trans_remove_enable, 'Value', 0);
        end

        if waveformViewerData.analyze_probing_events_value == 1
            set(handles.chkbox_analyze_probing_events,'Value', 1);
        else
            set(handles.chkbox_analyze_probing_events,'Value', 0);
        end

        handles.firstLoad = 1;

        load_current_sweep(hObject,handles)
        handles = guidata(hObject);
        set(handles.edt_sample_freq, 'Value', handles.time_tick_options_value);
    else
        errordlg('Incorrect File Type. Please find file ending in "_waveformViewerData.mat"', 'Error: Wrong File Type');
    end
end
guidata(hObject,handles);

% --- Executes on button press in chkbox_ignore_eject.
function chkbox_ignore_eject_Callback(hObject, eventdata, handles)
% hObject    handle to chkbox_ignore_eject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkbox_ignore_eject

if (get(handles.chkbox_ignore_eject, 'Value') == 1 )
    set(handles.edt_max_event_length, 'Visible', 'on');
    set(handles.text56, 'Visible', 'on');
else
    set(handles.edt_max_event_length, 'Visible', 'off');
    set(handles.text56, 'Visible', 'off');
end

% --- Executes on button press in btn_ramp_add_unbound.
function btn_ramp_add_unbound_Callback(hObject, eventdata, handles)
% hObject    handle to btn_ramp_add_unbound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if( handles.detected_events{handles.matIndex}(handles.event, 8) == 5 )
    handles.openChanRamps{handles.matIndex}(handles.event) = 1;
    set(handles.txt_ramp_is_unbound, 'BackgroundColor', [0 1 0]);
end
uicontrol(handles.listbox_sweeps);
guidata(hObject, handles);

% --- Executes on button press in btn_ramp_remove_unbound.
function btn_ramp_remove_unbound_Callback(hObject, eventdata, handles)
% hObject    handle to btn_ramp_remove_unbound (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if( handles.detected_events{handles.matIndex}(handles.event, 8) == 5 )
    handles.openChanRamps{handles.matIndex}(handles.event) = 0;
    set(handles.txt_ramp_is_unbound, 'BackgroundColor', [1 0 0]);
end
uicontrol(handles.listbox_sweeps);
guidata(hObject, handles);

% --- Executes on button press in btn_ramp_calc_average.
function btn_ramp_calc_average_Callback(hObject, eventdata, handles)
% hObject    handle to btn_ramp_calc_average (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%BEGIN OLD SECTION FOR CALCULATING AVERAGE RAMPS
oldindex = handles.curindex;
combined_ramps = [];
combined_voltage = [];
number_of_ramps = 0;
for i=1:handles.num_events
    handles.curindex = i;
    tempEventString = handles.event_string{handles.curindex};
    selectedSweepInfo = sscanf(tempEventString, '%d.%d');
    tmatIndex = selectedSweepInfo(1);
    tevent = selectedSweepInfo(2);
    if( handles.openChanRamps{tmatIndex}(tevent) == 1)
        [traceData, time, voltage] = get_event(handles, hObject);
        handles = guidata(hObject);
        disp(class(traceData))
        if( number_of_ramps == 0)
            combined_ramps = double(traceData);
            combined_voltage = double(voltage);
        else

            %combined_ramps = [combined_ramps [traceData;zeros(length(combined_ramps)-length(traceData),1)]];
            %combined_voltage = [combined_voltage [voltage; zeros(length(combined_ramps)-length(voltage),1)]];

            if( length(combined_ramps) >= length(traceData) )
                combined_ramps = combined_ramps(1:length(traceData)) + double(traceData);
                combined_voltage = combined_voltage(1:length(voltage)) + double(voltage);
            else
                combined_ramps = combined_ramps + double(traceData(1:length(combined_ramps)));
                combined_voltage = combined_voltage + double(voltage(1:length(combined_voltage)));
            end

        end

        number_of_ramps = number_of_ramps + 1;

    end

end
whos
handles.combined_ramps = combined_ramps./number_of_ramps;
combined_voltage = combined_voltage./number_of_ramps;
handles.combined_voltage = combined_voltage;
assignin('base', 'traceData', handles.combined_ramps)
figure(11);plot(handles.combined_ramps)
figure(12);plot(combined_voltage)
handles.curindex = oldindex;
%%%% END CALC AVG RAMP SECTION

load_current_sweep(hObject, handles);
handles = guidata(hObject);
guidata(hObject, handles);

% --- Executes on button press in btn_ramp_save.
function btn_ramp_save_Callback(hObject, eventdata, handles)
% hObject    handle to btn_ramp_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%%%%%%% BEGIN SAVING OF DATA

temp_processed_data = cell2mat(handles.sweep_data{:}(:,1));

median_dur = median(temp_processed_data(:,6));
% events_with_ts_before_kink = find (temp_processed_data(:,1) == 1 & temp_processed_data(:,6) > (median_dur-0.005) & temp_processed_data(:,8) == 5 & temp_processed_data(:,12) == 0)
% events_with_ts_after_kink = find (temp_processed_data(:,1) == 1 & temp_processed_data(:,6) > (median_dur-0.005) & temp_processed_data(:,8) == 5 & temp_processed_data(:,12) == 1);
% events_wo_ts = find( (temp_processed_data(:,1) <= 0) &
% temp_processed_data(:,6) > (median_dur-0.005)  );
dissocation_time_before_kink = [];
dissocation_time_after_kink = [];
voltage_at_dissocation_before_kink = [];
voltage_at_dissocation_after_kink = [];

for filenum = 1:handles.numOpenMATFiles;
    temp_processed_data_2 = cell2mat(handles.sweep_data{filenum}(:,1));
    events_with_ts_before_kink = find (temp_processed_data_2(:,1) == 1 ...
        & handles.detected_events{filenum}(:,8) == 5 & temp_processed_data_2(:,12) == 0 & temp_processed_data_2(:,6) > (median_dur-0.005))
    events_with_ts_after_kink = find (temp_processed_data_2(:,1) == 1 ...
        & handles.detected_events{filenum}(:,8) == 5 & temp_processed_data_2(:,12) == 1 & temp_processed_data_2(:,6) > (median_dur-0.005))
    dissocation_time_before_kink = [dissocation_time_before_kink; temp_processed_data_2(events_with_ts_before_kink,9)];
    dissocation_time_after_kink = [dissocation_time_after_kink; temp_processed_data_2(events_with_ts_after_kink,9)];
    voltage_at_dissocation_before_kink = [voltage_at_dissocation_before_kink; temp_processed_data_2(events_with_ts_before_kink,13)];
    voltage_at_dissocation_after_kink = [voltage_at_dissocation_after_kink; temp_processed_data_2(events_with_ts_after_kink,13)];
end
if(~isempty(num2str(get(handles.edt_ramp_rate,'String'))))
    save([get(handles.edt_ramp_rate,'String') '_ramp_dissociation.mat'], 'dissocation_time_before_kink', 'dissocation_time_after_kink', 'voltage_at_dissocation_before_kink', 'voltage_at_dissocation_after_kink');
else
    save('0_ramp_dissociation.mat', 'dissocation_time_before_kink', 'dissocation_time_after_kink', 'voltage_at_dissocation_before_kink', 'voltage_at_dissocation_after_kink');
end
%%%%%% END SAVING OF DATA

% function edtWinLen_Callback(hObject, eventdata, handles)
% % hObject    handle to edtWinLen (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
%
% % Hints: get(hObject,'String') returns contents of edtWinLen as text
% %        str2double(get(hObject,'String')) returns contents of edtWinLen as a double
%
% try
%     if str2num(get(hObject,'String')) > 1
%         set(handles.mnu_multi_step_option, 'Visible', 'on', 'Enable', 'on');
%     else
%         set(handles.mnu_multi_step_option, 'Enable', 'off', 'Visible', 'off');
%     end
% catch
%     set(handles.mnu_multi_step_option, 'Visible', 'off', 'Enable', 'off');
% end

%% ----------------Create/button down/unsused callbacks--------------------------

function mnu_settings_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Not being used right now, but might be useful as another area to input
% parameters

set(handles.settings,'Visible','on');
uiwait(handles.settings);
temp_data = guidata(handles.settings);
handles.alpha = str2num(get(temp_data.edt_alpha,'String'));
handles.thresh_amp = str2num(get(temp_data.edt_thresh_amp,'String'));
handles.baseline_window = str2num(get(temp_data.edt_baseline_window,'String'));
handles.baseline_window_default = str2num(get(temp_data.edt_baseline_window,'String'));
% handles.edt_exp_descrip = get(temp_data.edt_exp_descrip,'String');
% handles.edt_exp_voltage = str2num(get(temp_data.edt_exp_voltage,'String'));
% handles.edt_exp_stn = get(temp_data.edt_exp_stn, 'String');
handles.rem_trans = get(temp_data.chkbox_rem_trans,'Value');
handles.transAddSamples = str2num( get(temp_data.edtTransAddSamples, 'String') );

% if( isempty(get(temp_data.edt_filename, 'String')) )
%     handles.edt_save_filename = 'Summary';
% else
%     handles.edt_save_filename = get(temp_data.edt_filename,'String');
%     handles.edt_save_filename = strrep(handles.edt_save_filename,' ','_');
% end
if( get(temp_data.chkbox_lim_cur,'Value') )
    handles.min_current = str2num(get(temp_data.edt_cur_min,'String'));
    handles.max_current = str2num(get(temp_data.edt_cur_max,'String'));
else
    handles.min_current = 0;
    handles.max_current = 50;
end
if( get(temp_data.chkbox_event_trans,'Value') )
    handles.event_trans = 1;
else
    handles.event_trans = 0;
end
% switch get(temp_data.popmnu_data_type,'Value')
%     case 1
%         handles.data_type = 'binary_ternary';
%     case 2
%         handles.data_type = 'DNA_alone';
% end

switch get(temp_data.popmnu_baseline_type,'Value')
    case 1
        handles.baseline_type = '\%';
    case 2
        handles.baseline_type = ' ms';
end
guidata(hObject,handles);

function edt_saturation_sample_time_Callback(hObject, eventdata, handles)
% hObject    handle to edt_saturation_sample_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_saturation_sample_time as text
%        str2double(get(hObject,'String')) returns contents of edt_saturation_sample_time as a double

function edt_max_trans_length_Callback(hObject, eventdata, handles)
% hObject    handle to edt_max_trans_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_max_trans_length as text
%        str2double(get(hObject,'String')) returns contents of edt_max_trans_length as a double

function chkbox_analyze_probing_events_Callback(hObject, eventdata, handles)
% hObject    handle to chkbox_analyze_probing_events (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of
% chkbox_analyze_probing_events

function edt_beginning_trim_Callback(hObject, eventdata, handles)
% hObject    handle to edt_beginning_trim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_beginning_trim as text
%        str2double(get(hObject,'String')) returns contents of edt_beginning_trim as a double

function edt_end_trim_Callback(hObject, eventdata, handles)
% hObject    handle to edt_end_trim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_end_trim as text
%        str2double(get(hObject,'String')) returns contents of edt_end_trim as a double

function edt_sample_freq_Callback(hObject, eventdata, handles)
% hObject    handle to edt_sample_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_sample_freq as text
%        str2double(get(hObject,'String')) returns contents of
%        edt_sample_freq as a double

function btn_add_ts_Callback(hObject, eventdata, handles)
% hObject    handle to btn_add_ts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Not used

function mnu_multi_step_option_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_multi_step_option (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function btn_proc_sel_Callback(hObject, eventdata, handles)
% hObject    handle to btn_proc_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function edt_ramp_rate_Callback(hObject, eventdata, handles)
% hObject    handle to edt_ramp_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_ramp_rate as text
%        str2double(get(hObject,'String')) returns contents of
%        edt_ramp_rate as a double

function btn_rem_ts_Callback(hObject, eventdata, handles)
% hObject    handle to btn_rem_ts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function btn_add_all_Callback(hObject, eventdata, handles)
% hObject    handle to btn_add_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function edt_max_event_length_Callback(hObject, eventdata, handles)
% hObject    handle to edt_max_event_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_max_event_length as text
%        str2double(get(hObject,'String')) returns contents of
%        edt_max_event_length as a double

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over btn_proc_all.
function btn_proc_all_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to btn_proc_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


function btn_proc_sel_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to btn_proc_sel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function btn_clear_analyzed_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to btn_clear_analyzed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function btn_xls_plot_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to btn_xls_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function btn_del_wv_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to btn_del_wv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function btn_calc_rms_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to btn_calc_rms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function btn_chk_term_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to btn_chk_term (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function btn_add_ts_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to btn_add_ts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function btn_rem_ts_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to btn_rem_ts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function btn_reset_zoom_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to btn_reset_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function btn_add_all_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to btn_add_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function chkbox_check_diff_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to chkbox_check_diff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function sweep_plot_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to sweep_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function edtMinStepSize_Callback(hObject, eventdata, handles)
% hObject    handle to edtVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtVar as text
%        str2double(get(hObject,'String')) returns contents of edtVar as a double



function lstOpenFiles_Callback(hObject, eventdata, handles)
% hObject    handle to lstOpenFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns lstOpenFiles contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstOpenFiles

function edt_exp_descrip_Callback(hObject, eventdata, handles)
% hObject    handle to edt_exp_descrip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_exp_descrip as text
%        str2double(get(hObject,'String')) returns contents of edt_exp_descrip as a double
% handles.edt_exp_descrip = get(hObject,'String');
% guidata(hObject, handles)

function edt_exp_voltage_Callback(hObject, eventdata, handles)
% hObject    handle to edt_exp_voltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_exp_voltage as text
%        str2double(get(hObject,'String')) returns contents of edt_exp_voltage as a double

% handles.edt_exp_voltage = str2num(get(hObject, 'String'));
% guidata(hObject, handles);

function edt_exp_stn_Callback(hObject, eventdata, handles)
% hObject    handle to edt_exp_stn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_exp_stn as text
%        str2double(get(hObject,'String')) returns contents of edt_exp_stn as a double

% handles.edt_exp_stn = get(hObject, 'String');
% guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edt_exp_descrip_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_exp_descrip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edt_beginning_trim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_beginning_trim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edt_sample_freq_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_sample_freq (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edt_exp_stn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_exp_stn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edt_exp_voltage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_exp_voltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function lstOpenFiles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstOpenFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edt_end_trim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_end_trim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function listbox_sweeps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_sweeps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edtMinStepSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edt_save_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_save_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function mnu_data_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mnu_data_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edt_saturation_sample_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_saturation_sample_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edt_max_trans_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_max_trans_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edt_max_event_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_max_event_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edt_ramp_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_ramp_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_primer_transient_trim_Callback(hObject, eventdata, handles)
% hObject    handle to edt_primer_transient_trim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_primer_transient_trim as text
%        str2double(get(hObject,'String')) returns contents of edt_primer_transient_trim as a double

load_current_sweep(hObject, handles);

%update from load_current_sweep
handles = guidata(hObject);

guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function edt_primer_transient_trim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_primer_transient_trim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_primer_hold_current_Callback(hObject, eventdata, handles)
% hObject    handle to edt_primer_hold_current (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_primer_hold_current as text
%        str2double(get(hObject,'String')) returns contents of edt_primer_hold_current as a double


% --- Executes during object creation, after setting all properties.
function edt_primer_hold_current_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_primer_hold_current (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edt_primer_open_channel_probe_Callback(hObject, eventdata, handles)
% hObject    handle to edt_primer_open_channel_probe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_primer_open_channel_probe as text
%        str2double(get(hObject,'String')) returns contents of edt_primer_open_channel_probe as a double


% --- Executes during object creation, after setting all properties.
function edt_primer_open_channel_probe_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_primer_open_channel_probe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in btn_primer_save_data.
function btn_primer_save_data_Callback(hObject, eventdata, handles)
% hObject    handle to btn_primer_save_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Plot/Display some of the results

temp_processed_data = [];

keyboard
for i=1:length(handles.sweep_data{:}(:,1))
    disp(i)
    temp_processed_data = [temp_processed_data; handles.sweep_data{:}{i,1}];
end

%This needs to be adjusted
%Will make a a graph to adjust stuff for outlier removal

%If our detectection time is under .1ms then we probably didn't reduce when
%we recgonize the DNA
if(~isempty(find(temp_processed_data(:,3) < 0.0001)))
    tempfig = figure;

    hist(temp_processed_data(find(temp_processed_data(:,2) < prctile(temp_processed_data(:,2),90)),2), 1000)
    [temp_cutoff_hard,y] = ginput(1)
    pause(0.01)
    close(tempfig)
    %       temp_cutoff_hard = 0.0076650;
    temp_processed_data_hard = temp_processed_data(find(temp_processed_data(:,2) < temp_cutoff_hard),:);
else
    %We'll just take the median detection time + %10
    temp_cutoff_hard = mode(temp_processed_data(:,3))*.1+mode(temp_processed_data(:,3))
    % temp_processed_data_soft = temp_processed_data(find(temp_processed_data(:,3) < temp_cutoff_soft),:);
    temp_processed_data_hard = temp_processed_data(find(temp_processed_data(:,3) < temp_cutoff_hard),:);
end

pause(0.01)
% temp_cutoff_hard = 2.6e-3;
% temp_cutoff_soft = 4.5e-3;

edges = -1.5:1:3.5;
% bins_soft = histc(temp_processed_data_soft(:,1),edges)./length(temp_processed_data_soft) .* 100;
% bins_soft(end) = [];

bins_hard = histc(temp_processed_data_hard(:,1),edges)./length(temp_processed_data_hard) .* 100;
bins_hard(end) = [];

%Look only at those that didn't dissociated but not eject, event type 3
% t_dur_data_soft = sort(temp_processed_data_soft(find(temp_processed_data_soft(:,1) == 3),4));
% dur_probs_soft = (length(t_dur_data_soft):-1:1)/length(t_dur_data_soft);

t_dur_data_hard = sort(temp_processed_data_hard(find(temp_processed_data_hard(:,1) == 3),4));
dur_probs_hard = (length(t_dur_data_hard):-1:1)/length(t_dur_data_hard);

handles.plotPrimer = figure;
maximize(gcf)
% subplot(2,3,1)
% barh(-1:3,bins_soft);
% set(gca,'YTickLabel',{'False Detect','DNA Escape','No Primer', 'Primer Eject','Primer Dissoc.'})
% xlabel('Percent of Events'); xlim([0 100]);
% title(['Event Types for Soft Cutoff, N=' num2str(length(temp_processed_data_soft)) ', Cutoff Diag. Time: ' num2str(temp_cutoff_soft.*1000) 'ms'])

subplot(2,3,1)
hist(temp_processed_data(find(temp_processed_data(:,3) < prctile(temp_processed_data(:,3),85)),3), 200);
hold on;
line([temp_cutoff_hard temp_cutoff_hard], ylim,'Color','r')
hold off;
title('DNA Detection Time, Line is Cutoff')
xlabel('Time (s)')
ylabel('Count')

subplot(2,3,4)
barh(-1:3,bins_hard);
set(gca,'YTickLabel',{'False Detect','DNA Escape','No Primer', 'Primer Eject','Primer Dissoc.'})
xlabel('Percent of Events'); xlim([0 100]);
title({['Event Types for Hard Cutoff, N=' num2str(length(temp_processed_data_hard))], [' Cutoff Diag. Time: ' num2str(temp_cutoff_hard.*1000) 'ms']})

% subplot(2,3,2)
% loglog(t_dur_data_soft,dur_probs_soft, '-rx', 'LineWidth', 2);
% xlabel('Duration')
% ylabel('Survival Probability')
% title('Survival Prob. for Primer Dissociation Events, Soft Cutoff')
% axis([0.0001 0.3 0.03 1])

subplot(2,3,5)
loglog(t_dur_data_hard,dur_probs_hard, '-rx', 'LineWidth', 2);
xlabel('Duration')
ylabel('Survival Probability')
title({'Survival Prob. for Primer Dissociation Events', 'Hard Cutoff'})
axis([0.0001 0.3 0.03 1])

% subplot(2,3,3)
% hist(t_dur_data_soft)
% ylabel('Number of Events')
% xlabel('Duration')
% title('Histogram of Dwell Time for Primer Dissociation Events, Soft Cutoff')
% xlim([0 0.3])

subplot(2,3,6)
edges_t_dur = linspace(0,.3, 200);
bar(edges_t_dur*1000,histc(t_dur_data_hard, edges_t_dur));
% hist(t_dur_data_hard)
xlim([-10 300])
title({'Histogram of Dwell Time for Primer Dissociation Events', 'Hard Cutoff'})
ylabel('Count')
xlabel('Dissociation Time (ms)')
%Need this for saving to eps
set(handles.plotPrimer, 'PaperPositionMode','auto')

primer_savefilename = ['Primer_' get(handles.edt_primer_hold_time,'String') 'ms, ' get(handles.edt_primer_hold_voltage,'String') 'mV, ' handles.allFilenames{1}{1}(1:end-4)]
if ( exist([primer_savefilename '.mat'],'file') )
    button = questdlg(['Do you want to overwrite ' primer_savefilename '?'],'Confirm Overwrite');
    switch button
        case 'Yes'
            %do nothing
        case 'No'
            return;
        case 'Cancel'
            return;
        otherwise
            return;
    end

end


saveas(handles.plotPrimer,[primer_savefilename '.eps'],'epsc2');

saveas(handles.plotPrimer,[primer_savefilename '.png'],'png');

% save(['Primer_' get(handles.edt_primer_hold_time,'String') 'ms, ' get(handles.edt_primer_hold_voltage,'String') 'mV, ' handles.allFilenames{1}{1}(1:end-4) '.mat'], ...
%     'temp_processed_data', 'temp_processed_data_soft', 'temp_processed_data_hard', 'temp_cutoff_soft', 'temp_cutoff_hard');

save([primer_savefilename '.mat'], ...
    'temp_processed_data','temp_processed_data_hard', 'temp_cutoff_hard');

function edt_primer_hold_time_Callback(hObject, eventdata, handles)
% hObject    handle to edt_primer_hold_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_primer_hold_time as text
%        str2double(get(hObject,'String')) returns contents of edt_primer_hold_time as a double


% --- Executes during object creation, after setting all properties.
function edt_primer_hold_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_primer_hold_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_primer_hold_voltage_Callback(hObject, eventdata, handles)
% hObject    handle to edt_primer_hold_voltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_primer_hold_voltage as text
%        str2double(get(hObject,'String')) returns contents of edt_primer_hold_voltage as a double


% --- Executes during object creation, after setting all properties.
function edt_primer_hold_voltage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_primer_hold_voltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in chkbox_primer_eject_working.
function chkbox_primer_eject_working_Callback(hObject, eventdata, handles)
% hObject    handle to chkbox_primer_eject_working (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkbox_primer_eject_working





function edtVar_Callback(hObject, eventdata, handles)
% hObject    handle to edtVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtVar as text
%        str2double(get(hObject,'String')) returns contents of edtVar as a double


% --- Executes during object creation, after setting all properties.
function edtVar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edtPVal_Callback(hObject, eventdata, handles)
% hObject    handle to edtPVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtPVal as text
%        str2double(get(hObject,'String')) returns contents of edtPVal as a double


% --- Executes during object creation, after setting all properties.
function edtPVal_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtPVal (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtWinLen_Callback(hObject, eventdata, handles)
% hObject    handle to edtWinLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtWinLen as text
%        str2double(get(hObject,'String')) returns contents of edtWinLen as a double


% --- Executes during object creation, after setting all properties.
function edtWinLen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtWinLen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in edt_label_steps.
function edt_label_steps_Callback(hObject, eventdata, handles)
% hObject    handle to edt_label_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edt_label_steps


% --- Executes on button press in edt_check_output.
function edt_check_output_Callback(hObject, eventdata, handles)
% hObject    handle to edt_check_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of edt_check_output


% --- Executes during object creation, after setting all properties.
function edt_label_steps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_label_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
    set(hObject,'Value',get(hObject,'Max')); % checked by default



% --- Executes during object creation, after setting all properties.
function edt_check_output_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_check_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
    set(hObject,'Value',get(hObject,'Max')); %checked by default


% --- Executes on button press in edt_ignore_term.
function edt_ignore_term_Callback(hObject, eventdata, handles)
% hObject    handle to edt_ignore_term (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 if get(hObject,'Value') == get(hObject,'Max')
    set(handles.term_algo,'Visible','on');
    contents = get(handles.term_algo,'String');
    val = contents{get(handles.term_algo,'Value')};
    if strcmp(val,'Simple')
        set(handles.termpa,'Visible','on');
        set(handles.termp,'Visible','on');
        set(handles.edt_term_amp,'Visible','on');
        set(handles.edt_term_p,'Visible','on');
    else
        set(handles.termpa,'Visible','off');
        set(handles.termp,'Visible','off');
        set(handles.edt_term_amp,'Visible','off');
        set(handles.edt_term_p,'Visible','off');
    end
 else
     set(handles.term_algo,'Visible','off');
     set(handles.termpa,'Visible','off');
     set(handles.termp,'Visible','off');
     set(handles.edt_term_amp,'Visible','off');        
     set(handles.edt_term_p,'Visible','off');
 end
% Hint: get(hObject,'Value') returns toggle state of edt_ignore_term


% --- Executes during object creation, after setting all properties.
function edt_ignore_term_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_ignore_term (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
    set(hObject,'Value',get(hObject,'Max'));


% --- Executes on selection change in term_algo.
function term_algo_Callback(hObject, eventdata, handles)
% hObject    handle to term_algo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns term_algo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from term_algo
    contents = get(hObject,'String');
    val = contents{get(hObject,'Value')};
    if strcmp(val,'Simple')
        set(handles.termpa,'Visible','on');
        set(handles.termp,'Visible','on');
        set(handles.edt_term_amp,'Visible','on');
        set(handles.edt_term_p,'Visible','on');
    else
        set(handles.termpa,'Visible','off');
        set(handles.termp,'Visible','off');
        set(handles.edt_term_amp,'Visible','off');
        set(handles.edt_term_p,'Visible','off');
    end


% --- Executes during object creation, after setting all properties.
function term_algo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to term_algo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_term_amp_Callback(hObject, eventdata, handles)
% hObject    handle to edt_term_amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_term_amp as text
%        str2double(get(hObject,'String')) returns contents of edt_term_amp as a double


% --- Executes during object creation, after setting all properties.
function edt_term_amp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_term_amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_term_p_Callback(hObject, eventdata, handles)
% hObject    handle to edt_term_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_term_p as text
%        str2double(get(hObject,'String')) returns contents of edt_term_p as a double


% --- Executes during object creation, after setting all properties.
function edt_term_p_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_term_p (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
