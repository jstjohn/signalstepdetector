function varargout = postProcessor(varargin)
% POSTPROCESSOR M-file for postProcessor.fig
%      POSTPROCESSOR, by itself, creates a new POSTPROCESSOR or raises the existing
%      singleton*.
%
%      H = POSTPROCESSOR returns the handle to a new POSTPROCESSOR or the handle to
%      the existing singleton*.
%
%      POSTPROCESSOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POSTPROCESSOR.M with the given input arguments.
%
%      POSTPROCESSOR('Property','Value',...) creates a new POSTPROCESSOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before postProcessor_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to postProcessor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help postProcessor

% Last Modified by GUIDE v2.5 06-May-2009 12:54:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @postProcessor_OpeningFcn, ...
                   'gui_OutputFcn',  @postProcessor_OutputFcn, ...
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


% --- Executes just before postProcessor is made visible.
function postProcessor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to postProcessor (see VARARGIN)

% Choose default command line output for postProcessor
handles.output = hObject;

warning off all

handles.defaultPathname = '';
handles.default_curaxis = [0.2 11000 0 50];
handles.curaxis = handles.default_curaxis;
handles.saved_data = [];
handles.rboxdata = [];
handles.rboxloc = [];
handles.rbox = [];
handles.onFigure = 0;

tempUnits = get(handles.figure1, 'Units');
set(handles.figure1, 'Units', 'normalized');
%normalize to be less than 1
tempPos = get(handles.figure1,'Position');
set(handles.figure1, 'Position', [0 0.6 tempPos(3) tempPos(4)]);
set(handles.figure1, 'Units', tempUnits);
set(handles.btn_fishing_down_molecule, 'Enable', 'off');
set(handles.btn_fishing_up_molecule, 'Enable', 'off');
set(handles.mnu_fishing_molecule_number, 'Enable', 'off');
set(handles.btn_exclude_molecule, 'Enable', 'off');
set(handles.txt_exclude, 'Enable', 'off');
set(handles.btn_exclude_preview, 'Enable', 'off');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes postProcessor wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = postProcessor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in btn_open_mat.
function btn_open_mat_Callback(hObject, eventdata, handles)
% hObject    handle to btn_open_mat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if( ~exist(handles.defaultPathname) )
    handles.defaultPathname = [pwd '\'];

end

[filename, pathname] = uigetfile('*.mat','Open Trace Data', [handles.defaultPathname]);


if isequal(filename,0) || isequal(pathname,0)
    %do nothing, user pressed cancel
    handles.defaultPathname = [pwd '\'];
else
    %update file info
    filetypeind = strfind(filename, '_waveformViewerData.mat');
    %in case they select the wrong file type
    if(~isempty(filetypeind))
        filename = [filename(1:filetypeind-1) '.mat'];
    end
    
    handles.defaultPathname = pathname;
    handles.pathname = pathname;
    handles.filename = filename;
    
    try
        load( [pathname filename] );

        if( exist('temp_proc', 'var') )
            try
                saved_data_cell = struct2cell(handles.saved_data);
                matching_data = find( strcmp(saved_data_cell(1,1,:),filename));
            catch
                matching_data = [];
            end
            
            try
                delete(handles.rbox)
            end
            
            handles.rboxdata = [];
            
            if ( isempty( matching_data ) )
                handles.curIndex = length(handles.saved_data) + 1;
                handles.saved_data(handles.curIndex).filename = filename;
                handles.saved_data(handles.curIndex).pathname = pathname;
                handles.saved_data(handles.curIndex).processed_data = temp_proc;
                handles.saved_data(handles.curIndex).chkbox_e_extract_pre_term_steps = get(handles.chkbox_e_extract_pre_term_steps, 'Value');
                handles.saved_data(handles.curIndex).chkbox_e_extract_term_steps = get(handles.chkbox_e_extract_term_steps, 'Value');
                handles.saved_data(handles.curIndex).chkbox_e_wo_term_steps = get(handles.chkbox_e_wo_term_steps, 'Value');
                handles.saved_data(handles.curIndex).chkbox_ue_w_term_steps = get(handles.chkbox_ue_w_term_steps, 'Value');
                handles.saved_data(handles.curIndex).chkbox_ue_wo_term_steps = get(handles.chkbox_ue_wo_term_steps, 'Value');
                handles.saved_data(handles.curIndex).chkbox_fishing_disp_indiv_mol = 0;
                
                handles.saved_data(handles.curIndex).rboxdata = [];
                handles.saved_data(handles.curIndex).rboxloc = [];
                handles.processed_data = temp_proc;

                if( exist([pathname filename(1:end-4) '_waveformViewerData.mat']) == 2)
                    load( [pathname filename(1:end-4) '_waveformViewerData.mat'] );
                    %reconstruct detected_events so that it is one giant
                    %matrix and the molecule numbers don't overlap

                    if (waveformViewerData.fishing_file{1} == 1)
                        for i=2:waveformViewerData.numOpenMATFiles
                            waveformViewerData.detected_events{i:end}(:,7) =  waveformViewerData.detected_events{i:end}(:,7) + waveformViewerData.num_molecules(i-1);
                        end
                    end
                    handles.saved_data(handles.curIndex).detected_events = cell2mat(waveformViewerData.detected_events(:));
                    try
                        handles.saved_data(handles.curIndex).num_molecules = length(unique(waveformViewerData.detected_events{:}(:,7)));
                    catch
                        handles.saved_data(handles.curIndex).num_molecules = 0;
                    end
                    handles.saved_data(handles.curIndex).fishing_data = waveformViewerData.fishing_file{1};
                    handles.saved_data(handles.curIndex).fishing_molecule_index = 1;
                    handles.detected_events = handles.saved_data(handles.curIndex).detected_events;
                    handles.saved_data(handles.curIndex).fishing_molecules_str = cellstr(num2str(unique(handles.detected_events(:,7))));
                    handles.num_molecules = handles.saved_data(handles.curIndex).num_molecules;
                    handles.fishing_data = waveformViewerData.fishing_file{1};
                else
                    handles.saved_data(handles.curIndex).detected_events = [];
                    handles.saved_data(handles.curIndex).num_molecules = 0;
                    handles.saved_data(handles.curIndex).fishing_data = 0;
                    handles.detected_events = [];
                    handles.num_molecules = 0;
                    handles.fishing_data = 0;
                end
            else
                handles.curIndex = matching_data;
                handles.processed_data = handles.saved_data(handles.curIndex).processed_data;
                handles.detected_events = handles.saved_data(handles.curIndex).detected_events;
                handles.num_molecules = handles.saved_data(handles.curIndex).num_molecules;
                handles.fishing_data = handles.saved_data(handles.curIndex).fishing_data;
                
                handles.rboxdata = handles.saved_data(handles.curIndex).rboxdata;
                handles.rboxloc = handles.saved_data(handles.curIndex).rboxloc;
                set(handles.chkbox_e_extract_pre_term_steps, 'Value', handles.saved_data(handles.curIndex).chkbox_e_extract_pre_term_steps);
                set(handles.chkbox_e_extract_term_steps, 'Value', handles.saved_data(handles.curIndex).chkbox_e_extract_term_steps);
                set(handles.chkbox_e_wo_term_steps, 'Value', handles.saved_data(handles.curIndex).chkbox_e_wo_term_steps);
                set(handles.chkbox_ue_w_term_steps, 'Value', handles.saved_data(handles.curIndex).chkbox_ue_w_term_steps);
                set(handles.chkbox_ue_wo_term_steps, 'Value',handles.saved_data(handles.curIndex).chkbox_ue_wo_term_steps);
                set(handles.chkbox_fishing_disp_indiv_mol, 'Value', handles.saved_data(handles.curIndex).chkbox_fishing_disp_indiv_mol);
            end

            set( handles.figure1, 'Name', ['Post Processor - ' filename] );

            handles.events_with_ts = find( (handles.processed_data(:,1) == 1) );
            handles.events_wo_ts = find( (handles.processed_data(:,1) <= 0) );

            load_axes(hObject,handles,'ue')
            handles = guidata( hObject );
            load_axes(hObject,handles,'e')
            handles = guidata( hObject );

            tempcell = struct2cell(handles.saved_data);

            set(handles.lstbox_saved_files, 'Enable', 'on')
            set(handles.lstbox_saved_files, 'String', tempcell(1,1,:))
            set(handles.btn_export_plots, 'Enable', 'on');
            set(handles.btn_del_file, 'Enable', 'on');
            set(handles.btn_clear_list, 'Enable', 'on');
            set(handles.btn_open_pdf, 'Enable', 'on')
            if(handles.fishing_data == 1)
                set(handles.uipanel_fishing_options, 'Visible', 'on');
            else
                set(handles.uipanel_fishing_options, 'Visible', 'off');
            end
            set(handles.lstbox_saved_files, 'Value', handles.curIndex );
        elseif(  exist('saved_data', 'var') )
            
            handles.saved_data = saved_data;
            handles.curIndex = 1;
            
            handles.processed_data = handles.saved_data(handles.curIndex).processed_data;
            handles.rboxdata = handles.saved_data(handles.curIndex).rboxdata;
            handles.rboxloc = handles.saved_data(handles.curIndex).rboxloc;
            handles.events_with_ts = find( (handles.processed_data(:,1) == 1) );
            handles.events_wo_ts = find( (handles.processed_data(:,1) <= 0) );
            
            tempcell = struct2cell(handles.saved_data);
            
            set(handles.lstbox_saved_files, 'Enable', 'on')
            set(handles.lstbox_saved_files, 'String', tempcell(1,1,:))
            set(handles.btn_export_plots, 'Enable', 'on');
            set(handles.btn_del_file, 'Enable', 'on');
            set(handles.btn_clear_list, 'Enable', 'on');
            set(handles.btn_open_pdf, 'Enable', 'on');
            
            set(handles.lstbox_saved_files, 'Value', handles.curIndex );
            set(handles.chkbox_e_extract_pre_term_steps, 'Value', handles.saved_data(handles.curIndex).chkbox_e_extract_pre_term_steps);
            set(handles.chkbox_e_extract_term_steps, 'Value', handles.saved_data(handles.curIndex).chkbox_e_extract_term_steps);
            set(handles.chkbox_e_wo_term_steps, 'Value', handles.saved_data(handles.curIndex).chkbox_e_wo_term_steps);
            set(handles.chkbox_ue_w_term_steps, 'Value', handles.saved_data(handles.curIndex).chkbox_ue_w_term_steps);
            set(handles.chkbox_ue_wo_term_steps, 'Value',handles.saved_data(handles.curIndex).chkbox_ue_wo_term_steps);
            
            set( handles.figure1, 'Name', ['Post Processor - ' handles.saved_data(handles.curIndex).filename] );

            load_axes(hObject,handles,'ue')
            handles = guidata( hObject );
            load_axes(hObject,handles,'e')
            handles = guidata( hObject );
            
        end
        
    catch
        rethrow(lasterror)
        set(handles.btn_export_plots, 'Enable', 'off');
        %don't do anything
    end
end

guidata(hObject, handles);

% --- Executes on button press in btn_del_file.
function btn_del_file_Callback(hObject, eventdata, handles)
% hObject    handle to btn_del_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.saved_data(handles.curIndex) = [];
tempStr = get(handles.lstbox_saved_files, 'String');
tempStr(handles.curIndex) = [];
set(handles.lstbox_saved_files, 'String', tempStr);
set(handles.lstbox_saved_files, 'Value', 1)
guidata(hObject, handles)
if ( ~isempty(handles.saved_data) )
    lstbox_saved_files_Callback(handles.lstbox_saved_files, [], handles)
else
    cla(handles.axes_unextracted)
    cla(handles.axes_extracted)
    handles.saved_data = [];
    handles.rboxdata = [];
    handles.rboxloc = [];
    handles.rbox = [];
    handles.onFigure = 0;
    handles.processed_data = [];
    handles.events_with_ts = [];
    handles.events_wo_ts = [];
    set(handles.btn_del_file, 'Enable', 'off');
    set(handles.btn_clear_list, 'Enable', 'off');
    set(handles.btn_open_pdf, 'Enable', 'off');
    set(handles.txt_median_dwell_time, 'String', [num2str([]) ' ms']);
    set(handles.txt_median_amp, 'String', [num2str(median([])) ' pA']);
    set(handles.txt_num_points, 'String', num2str(length([])));
    set(handles.txt_median_dwell_all, 'String', [num2str(median([])) ' ms']);
    set(handles.txt_mean_amp_all, 'String', [num2str(mean([])) ' pA']);
    set(handles.txt_std_dev_all, 'String', [num2str(std([])) ' pA']);
    set(handles.txt_num_points_all, 'String', num2str(sum([])));
    set( handles.figure1, 'Name', 'Post Processor' );
    
    guidata(hObject, handles);
end
    
% --- Executes on button press in btn_clear_list.
function btn_clear_list_Callback(hObject, eventdata, handles)
% hObject    handle to btn_clear_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
cla(handles.axes_unextracted)
cla(handles.axes_extracted)

handles.saved_data = [];
handles.rboxdata = [];
handles.rboxloc = [];
handles.rbox = [];
handles.onFigure = 0;
handles.processed_data = [];
handles.events_with_ts = [];
handles.events_wo_ts = [];
set(handles.lstbox_saved_files, 'String', '');
set(handles.lstbox_saved_files, 'Value', 0)
set(handles.btn_del_file, 'Enable', 'off');
set(handles.btn_clear_list, 'Enable', 'off');
set(handles.btn_open_pdf, 'Enable', 'off');
set(handles.txt_median_dwell_time, 'String', [num2str([]) ' ms']);
set(handles.txt_median_amp, 'String', [num2str(median([])) ' pA']);
set(handles.txt_num_points, 'String', num2str(length([])));
set(handles.txt_median_dwell_all, 'String', [num2str(median([])) ' ms']);
set(handles.txt_mean_amp_all, 'String', [num2str(mean([])) ' pA']);
set(handles.txt_std_dev_all, 'String', [num2str(std([])) ' pA']);
set(handles.txt_num_points_all, 'String', num2str(sum([])));
set(handles.figure1, 'Name', 'Post Processor' );

guidata(hObject, handles);

% --- Executes on button press in btn_open_pdf.
function btn_open_pdf_Callback(hObject, eventdata, handles)
% hObject    handle to btn_open_pdf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if( exist([ handles.saved_data(handles.curIndex).pathname handles.saved_data(handles.curIndex).filename(1:end-3) 'pdf']) )
    open([ handles.saved_data(handles.curIndex).pathname handles.saved_data(handles.curIndex).filename(1:end-3) 'pdf'])
end

% --- Executes on button press in btn_save.
function btn_save_Callback(hObject, eventdata, handles)
% hObject    handle to btn_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
saved_data = handles.saved_data;

uisave('saved_data')

% --- Executes on button press in chkbox_ue_wo_term_steps.
function chkbox_ue_wo_term_steps_Callback(hObject, eventdata, handles)
% hObject    handle to chkbox_ue_wo_term_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkbox_ue_wo_term_steps

if ( get(handles.chkbox_ue_wo_term_steps, 'Value') == 1 )
    set(handles.ph_ue_evnts_wo_ts, 'Visible', 'on');
    handles.saved_data(handles.curIndex).chkbox_ue_wo_term_steps = 1;
else
    set(handles.ph_ue_evnts_wo_ts, 'Visible', 'off');
    handles.saved_data(handles.curIndex).chkbox_ue_wo_term_steps = 0;
end

try 
    ydat = handles.rboxdata(2,:);
    xdat = handles.rboxdata(1,:);
    p1 = [xdat(1) ydat(1)];
    offset = [xdat(2) - xdat(1), ydat(3) - ydat(1)];
    if ( strcmp(handles.rboxloc, 'ue') ) 
        generate_stats(hObject, handles, [p1 offset],'ue')
    else
        generate_stats(hObject, handles, [p1 offset],'e')
    end

    handles = guidata(hObject);
end

guidata(hObject, handles);

% --- Executes on button press in chkbox_ue_w_term_steps.
function chkbox_ue_w_term_steps_Callback(hObject, eventdata, handles)
% hObject    handle to chkbox_ue_w_term_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkbox_ue_w_term_steps

if ( get(handles.chkbox_ue_w_term_steps, 'Value') == 1 )
    set(handles.ph_ue_evnts_with_ts, 'Visible', 'on');
    handles.saved_data(handles.curIndex).chkbox_ue_w_term_steps = 1;
else
    set(handles.ph_ue_evnts_with_ts, 'Visible', 'off');
    handles.saved_data(handles.curIndex).chkbox_ue_w_term_steps = 0;
end

try 
    ydat = handles.rboxdata(2,:);
    xdat = handles.rboxdata(1,:);
    p1 = [xdat(1) ydat(1)];
    offset = [xdat(2) - xdat(1), ydat(3) - ydat(1)];
    if ( strcmp(handles.rboxloc, 'ue') ) 
        generate_stats(hObject, handles, [p1 offset],'ue')
    else
        generate_stats(hObject, handles, [p1 offset],'e')
    end

    handles = guidata(hObject);
    
end

guidata(hObject, handles);

% --- Executes on button press in chkbox_e_extract_term_steps.
function chkbox_e_extract_term_steps_Callback(hObject, eventdata, handles)
% hObject    handle to chkbox_e_extract_term_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkbox_e_extract_term_steps

if ( get(handles.chkbox_e_extract_term_steps, 'Value') == 1 )
    set(handles.ph_e_term_steps, 'Visible', 'on');
    handles.saved_data(handles.curIndex).chkbox_e_extract_term_steps = 1;
else
    set(handles.ph_e_term_steps, 'Visible', 'off');
    handles.saved_data(handles.curIndex).chkbox_e_extract_term_steps = 0;
end

try 
    ydat = handles.rboxdata(2,:);
    xdat = handles.rboxdata(1,:);
    p1 = [xdat(1) ydat(1)];
    offset = [xdat(2) - xdat(1), ydat(3) - ydat(1)];
    if ( strcmp(handles.rboxloc, 'ue') ) 
        generate_stats(hObject, handles, [p1 offset],'ue')
    else
        generate_stats(hObject, handles, [p1 offset],'e')
    end

    handles = guidata(hObject);
    
end

guidata(hObject, handles);

% --- Executes on button press in chkbox_e_extract_pre_term_steps.
function chkbox_e_extract_pre_term_steps_Callback(hObject, eventdata, handles)
% hObject    handle to chkbox_e_extract_pre_term_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkbox_e_extract_pre_term_steps

if ( get(handles.chkbox_e_extract_pre_term_steps, 'Value') == 1 )
    set(handles.ph_e_pre_term_steps, 'Visible', 'on');
    set(handles.chkbox_disp_assoc, 'Visible', 'on');
    handles.saved_data(handles.curIndex).chkbox_e_extract_pre_term_steps = 1;
else
    set(handles.ph_e_pre_term_steps, 'Visible', 'off');
    set(handles.chkbox_disp_assoc, 'Value', 0);
    set(handles.chkbox_disp_assoc, 'Visible', 'off');
    handles.saved_data(handles.curIndex).chkbox_e_extract_pre_term_steps = 0;
end

try 
    ydat = handles.rboxdata(2,:);
    xdat = handles.rboxdata(1,:);
    p1 = [xdat(1) ydat(1)];
    offset = [xdat(2) - xdat(1), ydat(3) - ydat(1)];
    if ( strcmp(handles.rboxloc, 'ue') ) 
        generate_stats(hObject, handles, [p1 offset],'ue')
    else
        generate_stats(hObject, handles, [p1 offset],'e')
    end

    handles = guidata(hObject);
    
end

guidata(hObject, handles);

% --- Executes on button press in chkbox_e_wo_term_steps.
function chkbox_e_wo_term_steps_Callback(hObject, eventdata, handles)
% hObject    handle to chkbox_e_wo_term_steps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkbox_e_wo_term_steps

if ( get(handles.chkbox_e_wo_term_steps, 'Value') == 1 )
    set(handles.ph_e_evnts_wo_ts, 'Visible', 'on');
    handles.saved_data(handles.curIndex).chkbox_e_wo_term_steps = 1;
else
    set(handles.ph_e_evnts_wo_ts, 'Visible', 'off');
    handles.saved_data(handles.curIndex).chkbox_e_wo_term_steps = 0;
end

try 
    ydat = handles.rboxdata(2,:);
    xdat = handles.rboxdata(1,:);
    p1 = [xdat(1) ydat(1)];
    offset = [xdat(2) - xdat(1), ydat(3) - ydat(1)];
    if ( handles.axes_unextracted == gca) 
        generate_stats(hObject, handles, [p1 offset],'ue')
    else
        generate_stats(hObject, handles, [p1 offset],'e')
    end

    handles = guidata(hObject);
    
end

guidata(hObject, handles);

function load_axes(hObject, handles, axesToLoad)

if( strcmp( axesToLoad, 'ue') )
    
    %plot non-ts events
    axes(handles.axes_unextracted)
    cla(handles.axes_unextracted)
    hold on;
    grid on;

    if( get(handles.chkbox_fishing_disp_indiv_mol, 'Value') == 1 && handles.fishing_data == 1)
        h1 = semilogx( handles.processed_data(intersect(find( handles.detected_events(:,7) == handles.cur_molecule),  handles.events_wo_ts),6).*1000, ...
            handles.processed_data( intersect(find( handles.detected_events(:,7) == handles.cur_molecule), handles.events_wo_ts),5),'.m','MarkerSize', 10);
    else
        h1 = semilogx(handles.processed_data( handles.events_wo_ts,6).*1000, ...
            handles.processed_data( handles.events_wo_ts,5),'.m','MarkerSize', 10);
    end
    set(h1, 'ButtonDownFcn', {@plot_1_Callback, hObject});
    if ( get(handles.chkbox_ue_wo_term_steps, 'Value') == 0)
        set(h1, 'Visible', 'off');
    end
    handles.ph_ue_evnts_wo_ts = h1;
    
    if( get(handles.chkbox_fishing_disp_indiv_mol, 'Value') == 1 && handles.fishing_data == 1)
        h2 = semilogx( handles.processed_data(intersect(find( handles.detected_events(:,7) == handles.cur_molecule),  handles.events_with_ts),6).*1000, ...
            handles.processed_data( intersect(find( handles.detected_events(:,7) == handles.cur_molecule), handles.events_with_ts),5),'.b','MarkerSize', 10);
    else
        h2 = semilogx(handles.processed_data( handles.events_with_ts,6).*1000, ...
            handles.processed_data( handles.events_with_ts,5),'.b','MarkerSize', 10);
    end
    set(h2, 'ButtonDownFcn', {@plot_1_Callback, hObject});
    if ( get(handles.chkbox_ue_w_term_steps, 'Value') == 0)
        set(h2,'Visible','off');
    end
    handles.ph_ue_evnts_with_ts = h2;

    if ( ~isempty(handles.rboxdata) && strcmp(handles.rboxloc, 'ue') )
        h3 = plot(handles.rboxdata(1,:),handles.rboxdata(2,:),'-r','linewidth',3); 
        handles.rbox = h3;
    end
    
    axis(handles.curaxis);
    xlabel('Duration (ms)','FontSize', 10);
    ylabel('Amplitude (pA)','FontSize', 10);

elseif( strcmp( axesToLoad, 'e') )

    %plot ts events
    axes(handles.axes_extracted)
    cla(handles.axes_extracted)
    hold on;
    grid on;
    if( get(handles.chkbox_fishing_disp_indiv_mol, 'Value') == 1 && handles.fishing_data == 1)
        h1 = semilogx( handles.processed_data(intersect(find( handles.detected_events(:,7) == handles.cur_molecule),  handles.events_wo_ts),6).*1000, ...
            handles.processed_data( intersect(find( handles.detected_events(:,7) == handles.cur_molecule), handles.events_wo_ts),5),'.m','MarkerSize', 10);
    else
        h1 = semilogx(handles.processed_data( handles.events_wo_ts ,6).*1000,...
            handles.processed_data( handles.events_wo_ts ,5),'.m','MarkerSize', 10);
    end
    set(h1, 'ButtonDownFcn', {@plot_2_Callback, hObject});
    if ( get(handles.chkbox_e_wo_term_steps, 'Value') == 0)
        set(h1, 'Visible', 'off');
    end
    handles.ph_e_evnts_wo_ts = h1;
    
    if( get(handles.chkbox_fishing_disp_indiv_mol, 'Value') == 1 && handles.fishing_data == 1)
        h2 = semilogx( handles.processed_data(intersect(find( handles.detected_events(:,7) == handles.cur_molecule),  handles.events_with_ts),9).*1000, ...
            handles.processed_data( intersect(find( handles.detected_events(:,7) == handles.cur_molecule), handles.events_with_ts),8),'.k','MarkerSize', 10);
    else
        h2 = semilogx(handles.processed_data( handles.events_with_ts ,9).*1000,...
            handles.processed_data( handles.events_with_ts ,8),'.k','MarkerSize', 10);
    end
    set(h2, 'ButtonDownFcn', {@plot_2_Callback, hObject});
    if( get(handles.chkbox_e_extract_pre_term_steps, 'Value') == 0)
        set(h2, 'Visible', 'off');
    end
    handles.ph_e_pre_term_steps = h2;

    if( get(handles.chkbox_fishing_disp_indiv_mol, 'Value') == 1 && handles.fishing_data == 1)
        h3 = semilogx( handles.processed_data(intersect(find( handles.detected_events(:,7) == handles.cur_molecule),  handles.events_with_ts),3).*1000, ...
            handles.processed_data( intersect(find( handles.detected_events(:,7) == handles.cur_molecule), handles.events_with_ts),2),'.c','MarkerSize', 10);
    else
        h3 = semilogx(handles.processed_data( handles.events_with_ts ,3).*1000, ...
            handles.processed_data( handles.events_with_ts ,2),'.c','MarkerSize', 10);
    end
    set(h3, 'ButtonDownFcn', {@plot_2_Callback, hObject});
    if( get(handles.chkbox_e_extract_term_steps, 'Value') == 0)
        set(h3, 'Visible', 'off');
    end
    handles.ph_e_term_steps = h3;
    
    if ( ~isempty(handles.rboxdata) && strcmp(handles.rboxloc, 'e') )
        h3 = plot(handles.rboxdata(1,:),handles.rboxdata(2,:),'-r','linewidth',3); 
        handles.rbox = h3;
    end
    axis(handles.curaxis);
    xlabel('Duration (ms)','FontSize', 10);
    ylabel('Amplitude (pA)','FontSize', 10);
    
end

guidata(hObject,handles);

% --- Executes on mouse press over axes background.
function axes_unextracted_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_unextracted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ( handles.onFigure == 0 )
    handles.onFigure = 9;
    if( ~isempty(handles.rboxdata)  && get(handles.btn_zoom, 'Value') == 0 )
        set(handles.rbox, 'Visible', 'off');
    end
    set( handles.figure1, 'Pointer', 'cross' );
    point1 = get(handles.axes_unextracted,'CurrentPoint');    % button down detected

    curAxisSize = axis(handles.axes_unextracted);
    axesPosition = get(handles.axes_unextracted,'Position');

    xPointNorm = (log(point1(1,1))-log(curAxisSize(1)))/(log(curAxisSize(2))-log(curAxisSize(1)));
    xPointChar = axesPosition(3)*xPointNorm + axesPosition(1);

    yPointNorm = (point1(1,2) - curAxisSize(3))/(curAxisSize(4) - curAxisSize(3));
    yPointChar = axesPosition(4)*yPointNorm + axesPosition(2);

%     load_axes(hObject,handles,'ue')

    handles = guidata(hObject);
    axes(handles.axes_unextracted)

    finalRect = rbbox([xPointChar yPointChar 0 0]);                   % return figure units
    point2 = get(handles.axes_unextracted,'CurrentPoint');    % button up detected
    set( handles.figure1, 'Pointer', 'arrow' );

    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    p1 = min(point1,point2);             % calculate locations
    offset = abs(point1-point2);         % and dimensions
    x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
    y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
    hold on
    axis manual
    h3 = plot(x,y,'-r','linewidth',3);                            % redraw in dataspace units
    handles.onFigure = 0;
    
    if( get(handles.btn_zoom, 'Value') == 1)
        delete(h3);
        handles.curaxis = [p1(1) p1(1)+offset(1) p1(2) p1(2)+offset(2)];
        axis(handles.curaxis)
        guidata(hObject,handles);
    else
        handles.rbox = h3;
        handles.rboxloc = 'ue';
        handles.rboxdata = [x;y];
        handles.saved_data(handles.curIndex).rboxdata = [x;y];
        handles.saved_data(handles.curIndex).rboxloc = 'ue';
        guidata(hObject,handles);
        generate_stats(hObject, handles, [p1 offset],'ue')
    end
end

function plot_1_Callback(handle,evt,hObject)

handles = guidata( hObject );

if ( ~isempty(handles.rboxdata) && get(handles.btn_zoom, 'Value') == 0 )
	delete(handles.rbox)
end

point1 = get(handles.axes_unextracted,'CurrentPoint');    % button down detected

curAxisSize = axis(handles.axes_unextracted);
axesPosition = get(handles.axes_unextracted,'Position');

xPointNorm = (log(point1(1,1))-log(curAxisSize(1)))/(log(curAxisSize(2))-log(curAxisSize(1)));
xPointChar = axesPosition(3)*xPointNorm + axesPosition(1);

yPointNorm = (point1(1,2) - curAxisSize(3))/(curAxisSize(4) - curAxisSize(3));
yPointChar = axesPosition(4)*yPointNorm + axesPosition(2);

% load_axes(hObject,handles,'ue')

handles = guidata(hObject);

axes(handles.axes_unextracted)

finalRect = rbbox([xPointChar yPointChar 0 0]);                   % return figure units
point2 = get(handles.axes_unextracted,'CurrentPoint');    % button up detected
set( handles.figure1, 'Pointer', 'arrow' );

point1 = point1(1,1:2);              % extract x and y
point2 = point2(1,1:2);
p1 = min(point1,point2);             % calculate locations
offset = abs(point1-point2);         % and dimensions
x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
hold on
axis manual
h4 = plot(x,y,'-r','linewidth',3);
handles.onFigure = 0;

if( get(handles.btn_zoom, 'Value') == 1)
    delete(h4);
    handles.curaxis = [p1(1) p1(1)+offset(1) p1(2) p1(2)+offset(2)];
    axis(handles.curaxis);
    guidata(hObject,handles);
else
    handles.rbox = h4;
    handles.rboxloc = 'ue';
	handles.rboxdata = [x;y];
    handles.saved_data(handles.curIndex).rboxdata = [x;y];
    handles.saved_data(handles.curIndex).rboxloc = 'ue';
    guidata(hObject,handles);
    generate_stats(hObject, handles, [p1 offset],'ue')
end


% --- Executes on mouse press over axes background.
function axes_extracted_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_extracted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ( handles.onFigure == 0 )
    if ( ~isempty(handles.rboxdata)  && get(handles.btn_zoom, 'Value') == 0 )
        delete(handles.rbox)
    end
    handles.onFigure = 9;
    set( handles.figure1, 'Pointer', 'cross' );
    point1 = get(handles.axes_extracted,'CurrentPoint');    % button down detected

    curAxisSize = axis(handles.axes_extracted);
    axesPosition = get(handles.axes_extracted,'Position');

    %return the 'percentage' along the xaxis
    xPointNorm = (log(point1(1,1))-log(curAxisSize(1)))/(log(curAxisSize(2))-log(curAxisSize(1)));
    xPointChar = axesPosition(3)*xPointNorm + axesPosition(1);

    yPointNorm = (point1(1,2) - curAxisSize(3))/(curAxisSize(4) - curAxisSize(3));
    yPointChar = axesPosition(4)*yPointNorm + axesPosition(2);
    
%     load_axes(hObject,handles,'e')
    
    handles = guidata(hObject);
    
    axes(handles.axes_extracted)

    finalRect = rbbox([xPointChar yPointChar 0 0]);                   % return figure units
    point2 = get(handles.axes_extracted,'CurrentPoint');    % button up detected
    set( handles.figure1, 'Pointer', 'arrow' );

    point1 = point1(1,1:2);              % extract x and y
    point2 = point2(1,1:2);
    p1 = min(point1,point2);             % calculate locations
    offset = abs(point1-point2);         % and dimensions
    x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
    y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
    hold on
    axis manual
    h3 = plot(x,y,'-r','linewidth',3);                            % redraw in dataspace units
    handles.onFigure = 0;
    
    if( get(handles.btn_zoom, 'Value') == 1)
        delete(h3);
        handles.curaxis = [p1(1) p1(1)+offset(1) p1(2) p1(2)+offset(2)];
        axis(handles.curaxis)
        guidata(hObject,handles);
    else
        handles.rbox = h3;
        handles.rboxloc = 'e';
        handles.rboxdata = [x;y];
        handles.saved_data(handles.curIndex).rboxdata = [x;y];
        handles.saved_data(handles.curIndex).rboxloc = 'e';
        guidata(hObject,handles);
        generate_stats(hObject, handles, [p1 offset],'e')
    end

end


function plot_2_Callback(handle,evt,hObject)

handles = guidata( hObject );

if ( ~isempty(handles.rboxdata)  && get(handles.btn_zoom, 'Value') == 0 )
	delete(handles.rbox)
end

point1 = get(handles.axes_extracted,'CurrentPoint');    % button down detected

curAxisSize = axis(handles.axes_extracted);
axesPosition = get(handles.axes_extracted,'Position');

xPointNorm = (log(point1(1,1))-log(curAxisSize(1)))/(log(curAxisSize(2))-log(curAxisSize(1)));
xPointChar = axesPosition(3)*xPointNorm + axesPosition(1);

yPointNorm = (point1(1,2) - curAxisSize(3))/(curAxisSize(4) - curAxisSize(3));
yPointChar = axesPosition(4)*yPointNorm + axesPosition(2);

% load_axes(hObject,handles,'e')

handles = guidata(hObject);

axes(handles.axes_extracted)

finalRect = rbbox([xPointChar yPointChar 0 0]);                   % return figure units
point2 = get(handles.axes_extracted,'CurrentPoint');    % button up detected
set( handles.figure1, 'Pointer', 'arrow' );

point1 = point1(1,1:2);              % extract x and y
point2 = point2(1,1:2);
p1 = min(point1,point2);             % calculate locations
offset = abs(point1-point2);         % and dimensions
x = [p1(1) p1(1)+offset(1) p1(1)+offset(1) p1(1) p1(1)];
y = [p1(2) p1(2) p1(2)+offset(2) p1(2)+offset(2) p1(2)];
hold on
axis manual
h4 = plot(x,y,'-r','linewidth',3);
set(h4, 'ButtonDownFcn', {@plot_2_Callback, hObject});

if( get(handles.btn_zoom, 'Value') == 1)
    delete(h4);
    handles.curaxis = [p1(1) p1(1)+offset(1) p1(2) p1(2)+offset(2)];
    axis(handles.curaxis)
    guidata(hObject,handles);
else
    handles.rbox = h4;
    handles.rboxloc = 'e';
    handles.rboxdata = [x;y];
    handles.saved_data(handles.curIndex).rboxdata = [x;y];
    handles.saved_data(handles.curIndex).rboxloc = 'e';
    guidata(hObject,handles);
    generate_stats(hObject, handles, [p1 offset],'e')
end

%---------------------------------------------------------------
function generate_stats(hObject, handles, boxCoords, curAxes)

if( strcmp( curAxes, 'ue'))
    
    events_in_box_wo_ts = [];
    events_in_box_with_ts = [];
    
    if ( get(handles.chkbox_ue_wo_term_steps, 'Value') == 1)
        events_in_box_wo_ts = find( handles.processed_data(:,6).*1000 >= boxCoords(1) & handles.processed_data(:,6).*1000 <= ( boxCoords(1)+boxCoords(3)) & handles.processed_data(:,5) >= boxCoords(2) & handles.processed_data(:,5) <= (boxCoords(2)+boxCoords(4)) & handles.processed_data(:,1) <= 0);
    end
    if ( get(handles.chkbox_ue_w_term_steps, 'Value') == 1)
        events_in_box_with_ts = find( handles.processed_data(:,6).*1000 >= boxCoords(1) & handles.processed_data(:,6).*1000 <= ( boxCoords(1)+boxCoords(3)) & handles.processed_data(:,5) >= boxCoords(2) & handles.processed_data(:,5) <= (boxCoords(2)+boxCoords(4)) & handles.processed_data(:,1) == 1);
    end

    temp = [events_in_box_wo_ts(:); events_in_box_with_ts];

    handles.saved_data(handles.curIndex).ampInBox = handles.processed_data(temp,5);
    handles.saved_data(handles.curIndex).durInBox = handles.processed_data(temp,6).*1000;
    handles.saved_data(handles.curIndex).ptsInBox = length(temp);
    
%     disp(median(handles.processed_data([events_in_box_wo_ts; events_in_box_with_ts],6)))
    set(handles.txt_median_dwell_time, 'String', [num2str(median(handles.processed_data(temp,6).*1000)) ' ms']);
    set(handles.txt_median_amp, 'String', [num2str(mean(handles.processed_data(temp,5))) ' pA']);
    set(handles.txt_num_points, 'String', num2str(length(temp)));
elseif( strcmp( curAxes, 'e'))
    
    events_in_box_wo_ts = [];
    events_in_box_term_step = [];
    events_in_box_pre_term_step = [];
    
    if ( get(handles.chkbox_e_wo_term_steps, 'Value') == 1)
        if( get(handles.chkbox_fishing_disp_indiv_mol, 'Value') == 1 && handles.fishing_data == 1)
            events_in_box_wo_ts = find( handles.processed_data(:,6).*1000 >= boxCoords(1) & ...
                handles.processed_data(:,6).*1000 <= ( boxCoords(1)+boxCoords(3)) & ...
                handles.processed_data(:,5) >= boxCoords(2) & ...
                handles.processed_data(:,5) <= (boxCoords(2)+boxCoords(4)) & ...
                handles.processed_data(:,1) <= 0); 
            events_in_box_wo_ts = intersect(find( handles.detected_events(:,7) == handles.cur_molecule), events_in_box_wo_ts);
        else
            events_in_box_wo_ts = find( handles.processed_data(:,6).*1000 >= boxCoords(1) & ...
                handles.processed_data(:,6).*1000 <= ( boxCoords(1)+boxCoords(3)) & ...
                handles.processed_data(:,5) >= boxCoords(2) & ...
                handles.processed_data(:,5) <= (boxCoords(2)+boxCoords(4)) & ...
                handles.processed_data(:,1) <= 0);
        end
    end
    if ( get(handles.chkbox_e_extract_term_steps, 'Value') == 1)
        if( get(handles.chkbox_fishing_disp_indiv_mol, 'Value') == 1 && handles.fishing_data == 1)
            events_in_box_term_step = find( handles.processed_data(:,3).*1000 >= boxCoords(1) & handles.processed_data(:,3).*1000 <= ( boxCoords(1)+boxCoords(3)) & handles.processed_data(:,2) >= boxCoords(2) & handles.processed_data(:,2) <= (boxCoords(2)+boxCoords(4)) & handles.processed_data(:,1) == 1);
            events_in_box_term_step = intersect(find( handles.detected_events(:,7) == handles.cur_molecule), events_in_box_term_step);
        else
            events_in_box_term_step = find( handles.processed_data(:,3).*1000 >= boxCoords(1) & handles.processed_data(:,3).*1000 <= ( boxCoords(1)+boxCoords(3)) & handles.processed_data(:,2) >= boxCoords(2) & handles.processed_data(:,2) <= (boxCoords(2)+boxCoords(4)) & handles.processed_data(:,1) == 1);
        end
    end
    if ( get(handles.chkbox_e_extract_pre_term_steps, 'Value') == 1)
        if( get(handles.chkbox_fishing_disp_indiv_mol, 'Value') == 1 && handles.fishing_data == 1)
            events_in_box_pre_term_step = find( handles.processed_data(:,9).*1000 >= boxCoords(1) & handles.processed_data(:,9).*1000 <= ( boxCoords(1)+boxCoords(3)) & handles.processed_data(:,8) >= boxCoords(2) & handles.processed_data(:,8) <= (boxCoords(2)+boxCoords(4)) & handles.processed_data(:,1) == 1); 
            events_in_box_pre_term_step = intersect(find( handles.detected_events(:,7) == handles.cur_molecule), events_in_box_pre_term_step);
        else
            events_in_box_pre_term_step = find( handles.processed_data(:,9).*1000 >= boxCoords(1) & handles.processed_data(:,9).*1000 <= ( boxCoords(1)+boxCoords(3)) & handles.processed_data(:,8) >= boxCoords(2) & handles.processed_data(:,8) <= (boxCoords(2)+boxCoords(4)) & handles.processed_data(:,1) == 1); 
        end
    end

    temp = [handles.processed_data(events_in_box_wo_ts(:),6).*1000; handles.processed_data(events_in_box_term_step(:),3).*1000; handles.processed_data(events_in_box_pre_term_step(:),9).*1000];
%     temp = [handles.processed_data(events_in_box_pre_term_step(:),9).*1000];
    temp2 = [handles.processed_data(events_in_box_wo_ts(:),5); handles.processed_data(events_in_box_term_step(:),2); handles.processed_data(events_in_box_pre_term_step(:),8)]; 
	temp3 = [handles.processed_data( events_in_box_pre_term_step(:) ,3).*1000];
    temp4 = [handles.processed_data( events_in_box_pre_term_step(:) ,2)];
    %temp stats finding the precentage of events pasts a line in the sand
    lineInSand = 10;  %in ms
    
    try
        set(handles.edt_temp_stats, 'Max', 2);
        set(handles.edt_temp_stats, 'String', {[num2str((length(find(temp > lineInSand))/length(temp))*100) ' %']; [num2str(median(temp(find(temp > lineInSand)))) ' ms']; ['Mean Amp gt LIS: ' num2str(mean(temp2(find(temp > lineInSand)))) ' pA']; ['Mean Amp ls LIS: ' num2str(mean(temp2(find(temp <= lineInSand)))) ' pA']; ['Num events gt LIS: ' num2str(length(find(temp > lineInSand)))]; ['Median Dwell: ' num2str(median(temp3))]})
    catch
        disp(lasterr)
    end
    set(handles.txt_median_dwell_time, 'String', [num2str(median(temp)) ' ms']);
    set(handles.txt_median_amp, 'String', [num2str(mean(temp2)) ' pA']);
    set(handles.txt_num_points, 'String', num2str(length(temp)));
    
    handles.saved_data(handles.curIndex).ptsInBox = length(temp);
    
    handles.saved_data(handles.curIndex).ampInBox = temp2;
    handles.saved_data(handles.curIndex).durInBox = temp;
    
    
    if ( ( get(handles.chkbox_disp_assoc, 'Value') == 1 ) && ( get(handles.chkbox_e_extract_pre_term_steps, 'Value') == 1 ) )
        h4 = semilogx(handles.processed_data( events_in_box_pre_term_step(:) ,3).*1000, ...
            handles.processed_data( events_in_box_pre_term_step(:) ,2),'.b','MarkerSize', 10);
        set(h4, 'ButtonDownFcn', {@plot_2_Callback, hObject});
    end

end

% size(handles.saved_data(1).durInBox)

set(handles.txt_median_dwell_all, 'String', [num2str(median(cat(1,handles.saved_data(:).durInBox))) ' ms']);
set(handles.txt_mean_amp_all, 'String', [num2str(mean(cat(1,handles.saved_data(:).ampInBox))) ' pA']);
set(handles.txt_std_dev_all, 'String', [num2str(std(cat(1,handles.saved_data(:).ampInBox))) ' pA']);
set(handles.txt_num_points_all, 'String', num2str(sum(cat(1,handles.saved_data(:).ptsInBox))));


guidata(hObject, handles)

% --- Executes on button press in btn_zoom.
function btn_zoom_Callback(hObject, eventdata, handles)
% hObject    handle to btn_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in btn_reset_zoom.
function btn_reset_zoom_Callback(hObject, eventdata, handles)
% hObject    handle to btn_reset_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.curaxis = handles.default_curaxis;
axis(handles.default_curaxis);
guidata(hObject,handles);

% --- Executes on button press in btn_export_plots.
function btn_export_plots_Callback(hObject, eventdata, handles)
% hObject    handle to btn_export_plots (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%update baseline from edit box
handles.axes5 = figure(1);

%init figure
set(handles.axes5, 'Position', [200 100 1000 450], 'Name', 'Unextracted Events Export');

if ( get(handles.chkbox_ue_wo_term_steps, 'Value') == 1)
    semilogx(handles.processed_data( handles.events_wo_ts,6).*1000, ...
        handles.processed_data( handles.events_wo_ts,5),'.m','MarkerSize', 10);
end
hold on;
if ( get(handles.chkbox_ue_w_term_steps, 'Value') == 1)
    semilogx(handles.processed_data( handles.events_with_ts,6).*1000, ...
        handles.processed_data( handles.events_with_ts,5),'.b','MarkerSize', 10);
end

set(gca, 'XScale', 'log');
set(gca, 'FontName', 'Arial', 'FontSize', 16)
% grid on;
axis([0.2 11000 0 50]);
xlabel('Duration (ms)','FontSize', 16, 'FontName', 'Arial');
ylabel('Amplitude (pA)','FontSize', 16, 'FontName', 'Arial');


%save figure as eps
maximize(gcf);

saveFilename = [handles.filename, ' - unextracted '];
% exportfig(gcf, [saveFilename, '.eps'], 'bounds', 'tight', 'FontSize', 1.8, 'Color', 'rgb');

handles.axes6 = figure(2);

%init figure
set(handles.axes6, 'Position', [200 100 1000 450], 'Name', 'Unextracted Events Export');

% if ( get(handles.chkbox_e_wo_term_steps, 'Value') == 1)
%     semilogx(handles.processed_data( handles.events_wo_ts ,6).*1000,...
%         handles.processed_data( handles.events_wo_ts ,5),'ow','MarkerSize', 5, 'Linewidth', 0.2, 'MarkerEdgeColor', 'k',...
%         'MarkerFaceColor', 'w');
% end
% hold on;
% if( get(handles.chkbox_e_extract_pre_term_steps, 'Value') == 1)
%     semilogx(handles.processed_data( handles.events_with_ts ,9).*1000,...
%         handles.processed_data( handles.events_with_ts ,8),'ok','MarkerSize', 5, 'Linewidth', 0.2, 'MarkerEdgeColor', 'k',...
%         'MarkerFaceColor', 'k');
% end
% if( get(handles.chkbox_e_extract_term_steps, 'Value') == 1)
%     semilogx(handles.processed_data( handles.events_with_ts ,3).*1000, ...
%         handles.processed_data( handles.events_with_ts ,2),'ok','MarkerSize', 5, 'Linewidth', 0.2, 'MarkerEdgeColor', 'k',...
%         'MarkerFaceColor', [.6 .6 .6]);
% end

if ( get(handles.chkbox_e_wo_term_steps, 'Value') == 1)
    semilogx(handles.processed_data( handles.events_wo_ts ,6).*1000,...
        handles.processed_data( handles.events_wo_ts ,5),'ow','MarkerSize', 5, 'Linewidth', 0.2, 'MarkerEdgeColor', 'm',...
        'MarkerFaceColor', 'm');
end
hold on;
if( get(handles.chkbox_e_extract_pre_term_steps, 'Value') == 1)
    semilogx(handles.processed_data( handles.events_with_ts ,9).*1000,...
        handles.processed_data( handles.events_with_ts ,8),'ok','MarkerSize', 5, 'Linewidth', 0.2, 'MarkerEdgeColor', 'k',...
        'MarkerFaceColor', 'k');
end
if( get(handles.chkbox_e_extract_term_steps, 'Value') == 1)
    semilogx(handles.processed_data( handles.events_with_ts ,3).*1000, ...
        handles.processed_data( handles.events_with_ts ,2),'ok','MarkerSize', 5, 'Linewidth', 0.2, 'MarkerEdgeColor', 'c',...
        'MarkerFaceColor', 'c');
end

set(gca, 'XScale', 'log');
set(gca, 'FontName', 'Arial', 'FontSize', 16)
% grid on;
axis([0.2 11000 10 50]);
xlabel('Duration (ms)','FontSize', 16, 'FontName', 'Arial');
ylabel('Amplitude (pA)','FontSize', 16, 'FontName', 'Arial');

%save figure as eps
maximize(gcf);

saveFilename = [handles.filename, ' - extracted '];
% exportfig(gcf, [saveFilename, '.eps'], 'bounds', 'tight', 'FontSize', 1.8, 'Color', 'rgb');

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in chkbox_disp_assoc.
function chkbox_disp_assoc_Callback(hObject, eventdata, handles)
% hObject    handle to chkbox_disp_assoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkbox_disp_assoc




% --- Executes on mouse motion over figure - except title and menu.
function figure1_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% disp(handles.onFigure)
if ( handles.onFigure ~= 9 )
    try
        ydat = handles.rboxdata(2,:);
        xdat = handles.rboxdata(1,:);
        curPos = get( gca, 'CurrentPoint' );
        curAxisSize = axis(gca);

        %Create a 'detection zone' for the current box

        percent_y = 0.025*(curAxisSize(4)-curAxisSize(3));
        percent_min_x = ((log(xdat(1))-log(curAxisSize(1)))/(log(curAxisSize(2))-log(curAxisSize(1))));
        percent_max_x = ((log(xdat(2))-log(curAxisSize(1)))/(log(curAxisSize(2))-log(curAxisSize(1))));


        outer_box_x(1) = exp(((percent_min_x - 0.01)*(log(curAxisSize(2)) - log(curAxisSize(1))) + log(curAxisSize(1))));
        outer_box_x(2) = exp(((percent_max_x + 0.01)*(log(curAxisSize(2)) - log(curAxisSize(1))) + log(curAxisSize(1))));

        inner_box_x(1) = exp(((percent_min_x + 0.01)*(log(curAxisSize(2)) - log(curAxisSize(1))) + log(curAxisSize(1))));
        inner_box_x(2) = exp(((percent_max_x - 0.01)*(log(curAxisSize(2)) - log(curAxisSize(1))) + log(curAxisSize(1))));

        outer_box_y(1) = ydat(1) - percent_y;
        outer_box_y(2) = ydat(3) + percent_y;

        inner_box_y(1) = ydat(1) + percent_y;
        inner_box_y(2) = ydat(3) - percent_y;


        bottom_box_x = zeros(1,5);
        bottom_box_y = zeros(1,5);

        right_box_x = zeros(1,5);
        right_box_y = zeros(1,5);

        top_box_x = zeros(1,5);
        top_box_y = zeros(1,5);

        left_box_x = zeros(1,5);
        left_box_y = zeros(1,5);

        bottom_box_x = [outer_box_x(1) outer_box_x(2) outer_box_x(2) outer_box_x(1) outer_box_x(1)];
        bottom_box_y = [outer_box_y(1) outer_box_y(1) inner_box_y(1) inner_box_y(1) outer_box_y(1)];

        right_box_x = [inner_box_x(2) outer_box_x(2) outer_box_x(2) inner_box_x(2) inner_box_x(2)];
        right_box_y = [outer_box_y(1) outer_box_y(1) outer_box_y(2) outer_box_y(2) outer_box_y(1)];

        top_box_x = [outer_box_x(1) outer_box_x(2) outer_box_x(2) outer_box_x(1) outer_box_x(1)];
        top_box_y = [inner_box_y(2) inner_box_y(2) outer_box_y(2) outer_box_y(2) inner_box_y(2)];

        left_box_x = [outer_box_x(1) inner_box_x(1) inner_box_x(1) outer_box_x(1) outer_box_x(1)];
        left_box_y = [outer_box_y(1) outer_box_y(1) outer_box_y(2) outer_box_y(2) outer_box_y(1)];

        if ( inpolygon(curPos(1,1), curPos(1,2), bottom_box_x, bottom_box_y) && ~(handles.onFigure == 1 || handles.onFigure == 2 || handles.onFigure == 3 ) )
            if(  inpolygon(curPos(1,1), curPos(1,2), right_box_x, right_box_y) )
                set(handles.txt_curPos, 'String', 'Yes');
                set(handles.figure1, 'Pointer', 'botr');
                if ( handles.onFigure == 1 )
                    set(handles.rbox, 'YData', [curPos(1,2) curPos(1,2) ydat(3) ydat(4) curPos(1,2)]);
                    set(handles.rbox, 'XData', [xdat(1) curPos(1,1) curPos(1,1) xdat(4) xdat(5)]);
                end
            elseif( inpolygon(curPos(1,1), curPos(1,2), left_box_x, left_box_y) )
                set(handles.txt_curPos, 'String', 'Yes')
                set(handles.figure1, 'Pointer', 'botl');
                if ( handles.onFigure == 2 )
                    set(handles.rbox, 'YData', [curPos(1,2) curPos(1,2) ydat(3) ydat(4) curPos(1,2)]);
                    set(handles.rbox, 'XData', [curPos(1,1) xdat(2) xdat(3) curPos(1,1) curPos(1,1)]);
                end
            else
                set(handles.txt_curPos, 'String', 'Yes');
                set( handles.figure1, 'Pointer', 'bottom' );
                if ( handles.onFigure == 3 )
                    set(handles.rbox, 'YData', [curPos(1,2) curPos(1,2) ydat(3) ydat(4) curPos(1,2)]);
                end
            end
        elseif ( handles.onFigure == 1 )
            set(handles.rbox, 'YData', [curPos(1,2) curPos(1,2) ydat(3) ydat(4) curPos(1,2)]);
            set(handles.rbox, 'XData', [xdat(1) curPos(1,1) curPos(1,1) xdat(4) xdat(5)]);
        elseif ( handles.onFigure == 2 )
            set(handles.rbox, 'YData', [curPos(1,2) curPos(1,2) ydat(3) ydat(4) curPos(1,2)]);
            set(handles.rbox, 'XData', [curPos(1,1) xdat(2) xdat(3) curPos(1,1) curPos(1,1)]);
        elseif ( handles.onFigure == 3 )
            set(handles.rbox, 'YData', [curPos(1,2) curPos(1,2) ydat(3) ydat(4) curPos(1,2)]);
        elseif ( inpolygon(curPos(1,1), curPos(1,2), top_box_x, top_box_y)  && ~(handles.onFigure == 4 || handles.onFigure == 5 || handles.onFigure == 6 ) )
            if(  inpolygon(curPos(1,1), curPos(1,2), right_box_x, right_box_y) )
                set(handles.txt_curPos, 'String', 'Yes');
                set(handles.figure1, 'Pointer', 'topr');
                if ( handles.onFigure == 4 )
                    set(handles.rbox, 'YData', [ydat(1) ydat(2) curPos(1,2) curPos(1,2) ydat(5)]);
                    set(handles.rbox, 'XData', [xdat(1) curPos(1,1) curPos(1,1) xdat(4) xdat(5)]);
                end
            elseif( inpolygon(curPos(1,1), curPos(1,2), left_box_x, left_box_y) )
                set(handles.txt_curPos, 'String', 'Yes')
                set(handles.figure1, 'Pointer', 'topl');
                if ( handles.onFigure == 5 )
                    set(handles.rbox, 'YData', [ydat(1) ydat(2) curPos(1,2) curPos(1,2) ydat(5)]);
                    set(handles.rbox, 'XData', [curPos(1,1) xdat(2) xdat(3) curPos(1,1) curPos(1,1)]);
                end
            else
                set(handles.txt_curPos, 'String', 'Yes');
                set( handles.figure1, 'Pointer', 'top' );
                if ( handles.onFigure == 6 )
                    set(handles.rbox, 'YData', [ydat(1) ydat(2) curPos(1,2) curPos(1,2) ydat(5)]);
                end
            end
        elseif ( handles.onFigure == 4 )
            set(handles.rbox, 'YData', [ydat(1) ydat(2) curPos(1,2) curPos(1,2) ydat(5)]);
            set(handles.rbox, 'XData', [xdat(1) curPos(1,1) curPos(1,1) xdat(4) xdat(5)]);
        elseif ( handles.onFigure == 5 )
            set(handles.rbox, 'YData', [ydat(1) ydat(2) curPos(1,2) curPos(1,2) ydat(5)])
            set(handles.rbox, 'XData', [curPos(1,1) xdat(2) xdat(3) curPos(1,1) curPos(1,1)]);
        elseif ( handles.onFigure == 6 )
            set(handles.rbox, 'YData', [ydat(1) ydat(2) curPos(1,2) curPos(1,2) ydat(5)]);
        elseif(  inpolygon(curPos(1,1), curPos(1,2), right_box_x, right_box_y) )
            set(handles.txt_curPos, 'String', 'Yes');
            set( handles.figure1, 'Pointer', 'right' );
            if ( handles.onFigure == 7 )
                set(handles.rbox, 'XData', [xdat(1) curPos(1,1) curPos(1,1) xdat(4) xdat(5)]);
            end
        elseif( inpolygon(curPos(1,1), curPos(1,2), left_box_x, left_box_y) )
            set(handles.txt_curPos, 'String', 'Yes');
            set( handles.figure1, 'Pointer', 'left' );
            if ( handles.onFigure == 8 )
                set(handles.rbox, 'XData', [curPos(1,1) xdat(2) xdat(3) curPos(1,1) curPos(1,1)]);
            end
        elseif ( handles.onFigure == 7 )
            set(handles.rbox, 'XData', [xdat(1) curPos(1,1) curPos(1,1) xdat(4) xdat(5)]);
        elseif ( handles.onFigure == 8 )
            set(handles.rbox, 'XData', [curPos(1,1) xdat(2) xdat(3) curPos(1,1) curPos(1,1)]);
        else
            set(handles.txt_curPos, 'String', 'No');
            set(handles.figure1, 'Pointer', 'arrow');
        end

        handles.saved_data(handles.curIndex).rboxdata = [get(handles.rbox, 'XData');get(handles.rbox, 'YData')];
        handles.rboxdata = [get(handles.rbox, 'XData');get(handles.rbox, 'YData')];
    catch
%         rethrow(lasterror)
        set(handles.txt_curPos, 'String', 'No')
        set(handles.figure1, 'Pointer', 'arrow');
    end

    guidata(hObject, handles);
end

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try
    ydat = handles.rboxdata(2,:);
    xdat = handles.rboxdata(1,:);
    curPos = get( gca, 'CurrentPoint' );
    curAxisSize = axis(gca);

    %Create a 'detection zone' for the current box
    
    percent_y = 0.025*(curAxisSize(4)-curAxisSize(3));
    percent_min_x = ((log(xdat(1))-log(curAxisSize(1)))/(log(curAxisSize(2))-log(curAxisSize(1))));
    percent_max_x = ((log(xdat(2))-log(curAxisSize(1)))/(log(curAxisSize(2))-log(curAxisSize(1))));
    
    
    outer_box_x(1) = exp(((percent_min_x - 0.01)*(log(curAxisSize(2)) - log(curAxisSize(1))) + log(curAxisSize(1))));
    outer_box_x(2) = exp(((percent_max_x + 0.01)*(log(curAxisSize(2)) - log(curAxisSize(1))) + log(curAxisSize(1))));
    
    inner_box_x(1) = exp(((percent_min_x + 0.01)*(log(curAxisSize(2)) - log(curAxisSize(1))) + log(curAxisSize(1))));
    inner_box_x(2) = exp(((percent_max_x - 0.01)*(log(curAxisSize(2)) - log(curAxisSize(1))) + log(curAxisSize(1))));
    
    outer_box_y(1) = ydat(1) - percent_y;
    outer_box_y(2) = ydat(3) + percent_y;
    
    inner_box_y(1) = ydat(1) + percent_y;
    inner_box_y(2) = ydat(3) - percent_y;

    
    bottom_box_x = zeros(1,5);
    bottom_box_y = zeros(1,5);
    
    right_box_x = zeros(1,5);
    right_box_y = zeros(1,5);
    
    top_box_x = zeros(1,5);
    top_box_y = zeros(1,5); 
    
    left_box_x = zeros(1,5);
    left_box_y = zeros(1,5);
    
    bottom_box_x = [outer_box_x(1) outer_box_x(2) outer_box_x(2) outer_box_x(1) outer_box_x(1)];
    bottom_box_y = [outer_box_y(1) outer_box_y(1) inner_box_y(1) inner_box_y(1) outer_box_y(1)];
    
    right_box_x = [inner_box_x(2) outer_box_x(2) outer_box_x(2) inner_box_x(2) inner_box_x(2)];
    right_box_y = [outer_box_y(1) outer_box_y(1) outer_box_y(2) outer_box_y(2) outer_box_y(1)];
    
    top_box_x = [outer_box_x(1) outer_box_x(2) outer_box_x(2) outer_box_x(1) outer_box_x(1)];
    top_box_y = [inner_box_y(2) inner_box_y(2) outer_box_y(2) outer_box_y(2) inner_box_y(2)]; 
    
    left_box_x = [outer_box_x(1) inner_box_x(1) inner_box_x(1) outer_box_x(1) outer_box_x(1)];
    left_box_y = [outer_box_y(1) outer_box_y(1) outer_box_y(2) outer_box_y(2) outer_box_y(1)];

    if ( inpolygon(curPos(1,1), curPos(1,2), bottom_box_x, bottom_box_y) )
        if(  inpolygon(curPos(1,1), curPos(1,2), right_box_x, right_box_y) )
            handles.onFigure = 1;
            set(handles.txt_curPos, 'String', 'Yes');
            set(handles.figure1, 'Pointer', 'botr');
            set(handles.rbox, 'YData', [curPos(1,2) curPos(1,2) ydat(3) ydat(4) curPos(1,2)]);
            set(handles.rbox, 'XData', [xdat(1) curPos(1,1) curPos(1,1) xdat(4) xdat(5)]);
        elseif( inpolygon(curPos(1,1), curPos(1,2), left_box_x, left_box_y) )
            handles.onFigure = 2;
            set(handles.txt_curPos, 'String', 'Yes')
            set(handles.figure1, 'Pointer', 'botl');
            set(handles.rbox, 'YData', [curPos(1,2) curPos(1,2) ydat(3) ydat(4) curPos(1,2)]);
            set(handles.rbox, 'XData', [curPos(1,1) xdat(2) xdat(3) curPos(1,1) curPos(1,1)]);
        else
            handles.onFigure = 3;
            set(handles.txt_curPos, 'String', 'Yes');
            set( handles.figure1, 'Pointer', 'bottom' );
            set(handles.rbox, 'YData', [curPos(1,2) curPos(1,2) ydat(3) ydat(4) curPos(1,2)]);
        end
    elseif ( inpolygon(curPos(1,1), curPos(1,2), top_box_x, top_box_y) )
        if(  inpolygon(curPos(1,1), curPos(1,2), right_box_x, right_box_y) )
            handles.onFigure = 4;
            set(handles.txt_curPos, 'String', 'Yes');
            set(handles.figure1, 'Pointer', 'topr');
            set(handles.rbox, 'YData', [ydat(1) ydat(2) curPos(1,2) curPos(1,2) ydat(5)]);
            set(handles.rbox, 'XData', [xdat(1) curPos(1,1) curPos(1,1) xdat(4) xdat(5)]);
        elseif( inpolygon(curPos(1,1), curPos(1,2), left_box_x, left_box_y) )
            handles.onFigure = 5;
            set(handles.txt_curPos, 'String', 'Yes')
            set(handles.figure1, 'Pointer', 'topl');
            set(handles.rbox, 'YData', [ydat(1) ydat(2) curPos(1,2) curPos(1,2) ydat(5)])
            set(handles.rbox, 'XData', [curPos(1,1) xdat(2) xdat(3) curPos(1,1) curPos(1,1)]);
        else
            handles.onFigure = 6;
            set(handles.txt_curPos, 'String', 'Yes');
            set( handles.figure1, 'Pointer', 'top' );
        end
    elseif(  inpolygon(curPos(1,1), curPos(1,2), right_box_x, right_box_y) )
        handles.onFigure = 7;
        set(handles.txt_curPos, 'String', 'Yes');
        set( handles.figure1, 'Pointer', 'right' );
        set(handles.rbox, 'XData', [xdat(1) curPos(1,1) curPos(1,1) xdat(4) xdat(5)]);
    elseif( inpolygon(curPos(1,1), curPos(1,2), left_box_x, left_box_y) )
        handles.onFigure = 8;
        set(handles.txt_curPos, 'String', 'Yes');
        set( handles.figure1, 'Pointer', 'left' );
        set(handles.rbox, 'XData', [curPos(1,1) xdat(2) xdat(3) curPos(1,1) curPos(1,1)]);
    else
        handles.onFigure = 0;
        set(handles.txt_curPos, 'String', 'No');
        set(handles.figure1, 'Pointer', 'arrow');
    end

     handles.saved_data(handles.curIndex).rboxdata =  [get(handles.rbox, 'XData');get(handles.rbox, 'YData')];
     handles.rboxdata =  [get(handles.rbox, 'XData');get(handles.rbox, 'YData')];

catch
    %     rethrow(lasterror)
    set(handles.txt_curPos, 'String', 'No')
    set(handles.figure1, 'Pointer', 'arrow');
end

guidata(hObject, handles);

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

try 
    ydat = handles.rboxdata(2,:);
    xdat = handles.rboxdata(1,:);
    p1 = [xdat(1) ydat(1)];
    offset = [xdat(2) - xdat(1), ydat(3) - ydat(1)];
    if ( handles.axes_unextracted == gca) 
        generate_stats(hObject, handles, [p1 offset],'ue')
    else
        generate_stats(hObject, handles, [p1 offset],'e')
    end

    handles = guidata(hObject);
    
end

handles.onFigure = 0;

guidata(hObject, handles)

% --- Executes on selection change in lstbox_saved_files.
function lstbox_saved_files_Callback(hObject, eventdata, handles)
% hObject    handle to lstbox_saved_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns lstbox_saved_files contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lstbox_saved_files
    
handles.curIndex = get(handles.lstbox_saved_files,'Value');

set( handles.figure1, 'Name', ['Post Processor - ' handles.saved_data(handles.curIndex).filename] );

handles.processed_data = handles.saved_data(handles.curIndex).processed_data;
handles.detected_events = handles.saved_data(handles.curIndex).detected_events;
handles.num_molecules = handles.saved_data(handles.curIndex).num_molecules;
handles.fishing_data = handles.saved_data(handles.curIndex).fishing_data;

handles.rboxdata = handles.saved_data(handles.curIndex).rboxdata;
handles.rboxloc = handles.saved_data(handles.curIndex).rboxloc;

handles.events_with_ts = find( (handles.processed_data(:,1) == 1) );
handles.events_wo_ts = find( (handles.processed_data(:,1) <= 0) );

set(handles.chkbox_e_extract_pre_term_steps, 'Value', handles.saved_data(handles.curIndex).chkbox_e_extract_pre_term_steps);
set(handles.chkbox_e_extract_term_steps, 'Value', handles.saved_data(handles.curIndex).chkbox_e_extract_term_steps);
set(handles.chkbox_e_wo_term_steps, 'Value', handles.saved_data(handles.curIndex).chkbox_e_wo_term_steps);
set(handles.chkbox_ue_w_term_steps, 'Value', handles.saved_data(handles.curIndex).chkbox_ue_w_term_steps);
set(handles.chkbox_ue_wo_term_steps, 'Value',handles.saved_data(handles.curIndex).chkbox_ue_wo_term_steps);
set(handles.chkbox_fishing_disp_indiv_mol, 'Value', handles.saved_data(handles.curIndex).chkbox_fishing_disp_indiv_mol);

if(handles.fishing_data == 1)
    set(handles.uipanel_fishing_options, 'Visible', 'on');
    set(handles.mnu_fishing_molecule_number, 'String', handles.saved_data(handles.curIndex).fishing_molecules_str);
    set(handles.mnu_fishing_molecule_number, 'Value', handles.saved_data(handles.curIndex).fishing_molecule_index);
    if(handles.saved_data(handles.curIndex).chkbox_fishing_disp_indiv_mol == 1)
        set(handles.btn_fishing_down_molecule, 'Enable', 'on');
        set(handles.btn_fishing_up_molecule, 'Enable', 'on');
        set(handles.mnu_fishing_molecule_number, 'Enable', 'on');
        set(handles.btn_exclude_molecule, 'Enable', 'on');
        set(handles.txt_exclude, 'Enable', 'on');
        set(handles.btn_exclude_preview, 'Enable', 'on');
        handles.cur_molecule = str2num(strtok(handles.saved_data(handles.curIndex).fishing_molecules_str{handles.saved_data(handles.curIndex).fishing_molecule_index},'*'));  
    else
        set(handles.btn_fishing_down_molecule, 'Enable', 'off');
        set(handles.btn_fishing_up_molecule, 'Enable', 'off');
        set(handles.mnu_fishing_molecule_number, 'Enable', 'off');
        set(handles.btn_exclude_molecule, 'Enable', 'off');
        set(handles.txt_exclude, 'Enable', 'off');
        set(handles.btn_exclude_preview, 'Enable', 'off');
    end
else
    set(handles.uipanel_fishing_options, 'Visible', 'off');
end

load_axes(hObject,handles,'ue')
handles = guidata( hObject );
load_axes(hObject,handles,'e')
handles = guidata( hObject );

if( ~isempty(handles.rboxdata) )
    xdat = handles.rboxdata(1,:);
    ydat = handles.rboxdata(2,:);
    p1 = [xdat(1) ydat(1)];
    offset = [xdat(2) - xdat(1), ydat(3) - ydat(1)];

    if ( strcmp(handles.rboxloc, 'ue') )
        axes(handles.axes_unextracted)
        generate_stats(hObject, handles, [p1 offset],'ue')
    else
        axes(handles.axes_extracted)
        generate_stats(hObject, handles, [p1 offset], 'e')
    end
    
    handles = guidata(hObject);
    
end

guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function lstbox_saved_files_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lstbox_saved_files (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_temp_stats_Callback(hObject, eventdata, handles)
% hObject    handle to edt_temp_stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_temp_stats as text
%        str2double(get(hObject,'String')) returns contents of edt_temp_stats as a double


% --- Executes during object creation, after setting all properties.
function edt_temp_stats_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_temp_stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in btn_custom_function.
function btn_custom_function_Callback(hObject, eventdata, handles)
% hObject    handle to btn_custom_function (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

lstbox_str = get(handles.lstbox_saved_files, 'String');

global t_dur_data_really_all;
global dur_probs_really_all;
global fity;
global voltage_str;

t_dur_data_really_all = {};
dur_probs_really_all = {};
fity = {};
voltage_str = {};
combine_probs = 1;
t_dur_data = [];
t_dur_data_all = [];

for file_inds = 1:length(lstbox_str)

    disp(file_inds);
    set(handles.lstbox_saved_files, 'Value', file_inds);

    lstbox_saved_files_Callback(hObject,[],handles);

    handles = guidata(hObject);

    drawnow

    if combine_probs == 0
        temptok = regexp(handles.saved_data(handles.curIndex).filename, '.*_Oligo_(.*)_.*_(.*)mV', 'tokens');
        voltage_str{file_inds} = ['Oligo ' temptok{1}{1} ' @ ' temptok{1}{2} 'mV'];

        if (get(handles.chkbox_fishing_disp_indiv_mol, 'Value') ~= 1)
            set(handles.chkbox_fishing_disp_indiv_mol, 'Value', 1);
            chkbox_fishing_disp_indiv_mol_Callback(hObject, [], handles);
            handles = guidata(hObject);
            drawnow
        end

        tstr = get(handles.mnu_fishing_molecule_number, 'String');
        cmap = colormap(lines);

        t_dur_data = [];
        t_dur_data_all = [];
        for i = find(~strncmpi('*', tstr, 1))'

            handles.cur_molecule =str2num(tstr{i});

            try
                ydat = handles.rboxdata(2,:);
                xdat = handles.rboxdata(1,:);
                p1 = [xdat(1) ydat(1)];
                offset = [xdat(2) - xdat(1), ydat(3) - ydat(1)];

                generate_stats(hObject, handles, [p1 offset],'e')

                handles = guidata(hObject);
                drawnow

                if (~isempty(handles.saved_data(handles.curIndex).durInBox))
                    t_dur_data = sort(handles.saved_data(handles.curIndex).durInBox);
                    disp(size(t_dur_data))
                    t_dur_data_all = [t_dur_data_all; handles.saved_data(handles.curIndex).durInBox];
                    dur_probs = (length(t_dur_data):-1:1)/length(t_dur_data);

                    %                 tcolor = [rand rand rand];
                    %                 loglog(t_dur_data,dur_probs, '.-', 'Color', tcolor); hold on;
                end

            end

        end

        % set(handles.chkbox_fishing_disp_indiv_mol, 'Value', 0);
        % generate_stats(hObject, handles, [p1 offset],'e');
        % handles = guidata(hObject);
        if (~isempty(t_dur_data_all))
            t_dur_data_all = sort(t_dur_data_all);
            t_dur_data_really_all{file_inds} = t_dur_data_all;
            dur_probs = (length(t_dur_data_all):-1:1)/length(t_dur_data_all);
            dur_probs_really_all{file_inds} = dur_probs;
        end
    else
        
        try
            ydat = handles.rboxdata(2,:);
            xdat = handles.rboxdata(1,:);
            p1 = [xdat(1) ydat(1)];
            offset = [xdat(2) - xdat(1), ydat(3) - ydat(1)];

            generate_stats(hObject, handles, [p1 offset],'e')

            handles = guidata(hObject);
            drawnow

            if (~isempty(handles.saved_data(handles.curIndex).durInBox))
                t_dur_data = sort(handles.saved_data(handles.curIndex).durInBox);
                disp(size(t_dur_data))
                t_dur_data_all = [t_dur_data_all; handles.saved_data(handles.curIndex).durInBox];
                dur_probs = (length(t_dur_data):-1:1)/length(t_dur_data);

                %                 tcolor = [rand rand rand];
                %                 loglog(t_dur_data,dur_probs, '.-', 'Color', tcolor); hold on;
            end
        end
    end
%     cfun = fit(t_dur_data_all, dur_probs', 'exp2');
%     fity{file_inds} = feval(cfun, t_dur_data_all);

    guidata(hObject, handles)
end
figure(20);
% keyboard

if combine_probs == 0
    for j = 1:length(t_dur_data_really_all)
        loglog(t_dur_data_really_all{j},dur_probs_really_all{j}, 'Color', cmap(j,:), 'LineWidth', 2);hold on;

    %     plot(t_dur_data_all, fity, '-b', 'LineWidth', 2);

    end
    legend(voltage_str)

else
    t_dur_data_all = sort(t_dur_data_all);
    dur_probs_all = (length(t_dur_data_all):-1:1)/length(t_dur_data_all);
    loglog(t_dur_data_all, dur_probs_all, 'LineWidth', 2);
end

xlabel('Time (ms)')
ylabel('P_{survival}')
axis([0.9 1000 0.002 1.1])

save('Survival Probs\22mer_biz.mat', 't_dur_data_all', 'dur_probs_all')
% figure(handles.curIndex+1);
% 
% t_dur_data = sort(handles.saved_data(handles.curIndex).durInBox./1e3);
% dur_probs = (length(t_dur_data):-1:1)/length(t_dur_data);
% loglog(t_dur_data,dur_probs, 'o');

% figure(handles.curIndex+1);
% [N,X] = hist(handles.saved_data(handles.curIndex).ampInBox, 30);
% TotEvents = sum(N);
% N = N./TotEvents;
% bar(X,N,'hist')

% bins = -.5:.2:4;
% figure(1);
% for i = 1:length(handles.saved_data)
%     
%   subplot(length(handles.saved_data), 1, i)
%   hist(log10(handles.saved_data(i).processed_data(find( handles.saved_data(i).processed_data(:,1) == 1), 9).*1000), bins)
% end

% --- Executes on button press in chkbox_fishing_disp_indiv_mol.
function chkbox_fishing_disp_indiv_mol_Callback(hObject, eventdata, handles)
% hObject    handle to chkbox_fishing_disp_indiv_mol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkbox_fishing_disp_indiv_mol

if ( get(handles.chkbox_fishing_disp_indiv_mol, 'Value') == 1 )
    set(handles.btn_fishing_down_molecule, 'Enable', 'on');
    set(handles.btn_fishing_up_molecule, 'Enable', 'on');
    set(handles.btn_exclude_molecule, 'Enable', 'on');
    set(handles.mnu_fishing_molecule_number, 'Enable', 'on');
    set(handles.txt_exclude, 'Enable', 'on');
    set(handles.btn_exclude_preview, 'Enable', 'on');
    set(handles.mnu_fishing_molecule_number, 'String', handles.saved_data(handles.curIndex).fishing_molecules_str);
    handles.saved_data(handles.curIndex).chkbox_fishing_disp_indiv_mol = 1;
    tstr = get(handles.mnu_fishing_molecule_number, 'String');
    handles.cur_molecule = str2num(tstr{1});
else
    set(handles.btn_fishing_down_molecule, 'Enable', 'off');
    set(handles.btn_fishing_up_molecule, 'Enable', 'off');
    set(handles.mnu_fishing_molecule_number, 'Enable', 'off');
    set(handles.btn_exclude_molecule, 'Enable', 'off');
    set(handles.txt_exclude, 'Enable', 'off');
    set(handles.btn_exclude_preview, 'Enable', 'off');
    handles.saved_data(handles.curIndex).chkbox_fishing_disp_indiv_mol = 0;
end

load_axes(hObject,handles,'ue')
handles = guidata( hObject );
load_axes(hObject,handles,'e')
handles = guidata( hObject );
try 
    ydat = handles.rboxdata(2,:);
    xdat = handles.rboxdata(1,:);
    p1 = [xdat(1) ydat(1)];
    offset = [xdat(2) - xdat(1), ydat(3) - ydat(1)];
    if ( handles.axes_unextracted == gca) 
        generate_stats(hObject, handles, [p1 offset],'ue')
    else
        generate_stats(hObject, handles, [p1 offset],'e')
    end

    handles = guidata(hObject);
    
end
guidata(hObject, handles);
% --- Executes on button press in btn_fishing_down_molecule.
function btn_fishing_down_molecule_Callback(hObject, eventdata, handles)
% hObject    handle to btn_fishing_down_molecule (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if( get(handles.mnu_fishing_molecule_number,'Value') ~= 1 )
    set(handles.mnu_fishing_molecule_number,'Value', get(handles.mnu_fishing_molecule_number, 'Value')-1);
    handles.saved_data(handles.curIndex).fishing_molecule_index = get(handles.mnu_fishing_molecule_number,'Value');
    tstr = get(handles.mnu_fishing_molecule_number, 'String');
    handles.cur_molecule = str2num(strtok(tstr{get(handles.mnu_fishing_molecule_number, 'Value')}, '*'));
    if ( isempty(strfind(tstr{get(handles.mnu_fishing_molecule_number, 'Value')}, '*')))
        set(handles.btn_exclude_molecule, 'String', 'Exclude')
    else
        set(handles.btn_exclude_molecule, 'String', 'Include')
    end
end
load_axes(hObject,handles,'ue')
handles = guidata( hObject );
load_axes(hObject,handles,'e')
handles = guidata( hObject );
try 
    ydat = handles.rboxdata(2,:);
    xdat = handles.rboxdata(1,:);
    p1 = [xdat(1) ydat(1)];
    offset = [xdat(2) - xdat(1), ydat(3) - ydat(1)];
    if ( handles.axes_unextracted == gca) 
        generate_stats(hObject, handles, [p1 offset],'ue')
    else
        generate_stats(hObject, handles, [p1 offset],'e')
    end

    handles = guidata(hObject);
    
end
guidata(hObject, handles);
% --- Executes on selection change in mnu_fishing_molecule_number.
function mnu_fishing_molecule_number_Callback(hObject, eventdata, handles)
% hObject    handle to mnu_fishing_molecule_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns mnu_fishing_molecule_number contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mnu_fishing_molecule_number

tstr = get(handles.mnu_fishing_molecule_number, 'String');
handles.cur_molecule = str2num(strtok(tstr{get(handles.mnu_fishing_molecule_number, 'Value')}, '*'));
handles.saved_data(handles.curIndex).fishing_molecule_index = get(handles.mnu_fishing_molecule_number,'Value');
if ( isempty(strfind(tstr{get(handles.mnu_fishing_molecule_number, 'Value')}, '*')))
    set(handles.btn_exclude_molecule, 'String', 'Exclude')
else
    set(handles.btn_exclude_molecule, 'String', 'Include')
end
load_axes(hObject,handles,'ue')
handles = guidata( hObject );
load_axes(hObject,handles,'e')
handles = guidata( hObject );
try 
    ydat = handles.rboxdata(2,:);
    xdat = handles.rboxdata(1,:);
    p1 = [xdat(1) ydat(1)];
    offset = [xdat(2) - xdat(1), ydat(3) - ydat(1)];
    if ( handles.axes_unextracted == gca) 
        generate_stats(hObject, handles, [p1 offset],'ue')
    else
        generate_stats(hObject, handles, [p1 offset],'e')
    end

    handles = guidata(hObject);
    
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function mnu_fishing_molecule_number_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mnu_fishing_molecule_number (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btn_fishing_up_molecule.
function btn_fishing_up_molecule_Callback(hObject, eventdata, handles)
% hObject    handle to btn_fishing_up_molecule (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if( get(handles.mnu_fishing_molecule_number,'Value') ~= handles.num_molecules )
    set(handles.mnu_fishing_molecule_number,'Value', get(handles.mnu_fishing_molecule_number,'Value') + 1)
    handles.saved_data(handles.curIndex).fishing_molecule_index = get(handles.mnu_fishing_molecule_number,'Value');
    tstr = get(handles.mnu_fishing_molecule_number, 'String');
    handles.cur_molecule = str2num(strtok(tstr{get(handles.mnu_fishing_molecule_number, 'Value')}, '*'));
    if ( isempty(strfind(tstr{get(handles.mnu_fishing_molecule_number, 'Value')}, '*')))
        set(handles.btn_exclude_molecule, 'String', 'Exclude')
    else
        set(handles.btn_exclude_molecule, 'String', 'Include')
    end
end
load_axes(hObject,handles,'ue')
handles = guidata( hObject );
load_axes(hObject,handles,'e')
handles = guidata( hObject );

try 
    ydat = handles.rboxdata(2,:);
    xdat = handles.rboxdata(1,:);
    p1 = [xdat(1) ydat(1)];
    offset = [xdat(2) - xdat(1), ydat(3) - ydat(1)];
    if ( handles.axes_unextracted == gca) 
        generate_stats(hObject, handles, [p1 offset],'ue')
    else
        generate_stats(hObject, handles, [p1 offset],'e')
    end

    handles = guidata(hObject);
    
end

guidata(hObject, handles);


% --- Executes on button press in btn_felix.
function btn_felix_Callback(hObject, eventdata, handles)
% hObject    handle to btn_felix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%calculate number of events before first bound event

disp(unique(handles.detected_events(:,7)'))
first_occurance = [];
for i = unique(handles.detected_events(:,7))'
    
    cur_molecule_inds = find( handles.detected_events(:,7) == i);
    ai = find( (handles.processed_data(cur_molecule_inds,1) == 1), 1 );
    
    if isempty(ai)
        first_occurance = [first_occurance 0];
    else
        first_occurance = [first_occurance ai(1)];  
    end
    
end
set(handles.txt_felix_nnz, 'String', num2str(nnz(first_occurance)))
figure(2); plot(unique(handles.detected_events(:,7))', first_occurance,'-xb');xlabel('Molecule Number'); ylabel('Number of probes till first step')
first_occurance(find(first_occurance == 0)) = [];
set(handles.txt_felix_mean, 'String', mean(mean(first_occurance)))
set(handles.txt_felix_90, 'String', prctile(first_occurance, 90))


% --- Executes on button press in btn_surv_prob.
function btn_surv_prob_Callback(hObject, eventdata, handles)
% hObject    handle to btn_surv_prob (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%plot non-ts events
axes(handles.axes_unextracted)
cla(handles.axes_unextracted)
hold on;
grid on;

if (get(handles.chkbox_fishing_disp_indiv_mol, 'Value') ~= 1)
    set(handles.chkbox_fishing_disp_indiv_mol, 'Value', 1);
    pause(0.1)
    chkbox_fishing_disp_indiv_mol_Callback(hObject, [], handles);
    handles = guidata(hObject);
    pause(0.1)
end

tstr = get(handles.mnu_fishing_molecule_number, 'String');
cmap = colormap;
figure(handles.curIndex*2);

t_dur_data_all = [];

for i = find(~strncmpi('*', tstr, 1))'

    handles.cur_molecule =str2num(tstr{i});

    try
        ydat = handles.rboxdata(2,:);
        xdat = handles.rboxdata(1,:);
        p1 = [xdat(1) ydat(1)];
        offset = [xdat(2) - xdat(1), ydat(3) - ydat(1)];

        generate_stats(hObject, handles, [p1 offset],'e')

        handles = guidata(hObject);
        
        if (~isempty(handles.saved_data(handles.curIndex).durInBox))
            t_dur_data = sort(handles.saved_data(handles.curIndex).durInBox);
            t_dur_data_all = [t_dur_data_all; handles.saved_data(handles.curIndex).durInBox];
            dur_probs = (length(t_dur_data):-1:1)/length(t_dur_data);

            tcolor = [rand rand rand];
            loglog(t_dur_data,dur_probs, '.-', 'Color', tcolor); hold on;
        end

    end

end

if (~isempty(t_dur_data_all))
    t_dur_data_all = sort(t_dur_data_all);
    dur_probs = (length(t_dur_data_all):-1:1)/length(t_dur_data_all);

    loglog(t_dur_data_all,dur_probs, '-r', 'LineWidth', 5);
end

cfun = fit(t_dur_data_all, dur_probs', 'exp2');
fity = feval(cfun, t_dur_data_all);
plot(t_dur_data_all, fity, '-b', 'LineWidth', 2);
xlabel('Time (ms)')
ylabel('P_{survival}')
axis([0.9 1000 0.002 1.1])
% assignin('base', 't_dur_data', t_dur_data)
% assignin('base', 'dur_probs', dur_probs)



% --- Executes on button press in btn_exclude_molecule.
function btn_exclude_molecule_Callback(hObject, eventdata, handles)
% hObject    handle to btn_exclude_molecule (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tstr = get(handles.mnu_fishing_molecule_number, 'String');

cur_molecule_str = tstr{get(handles.mnu_fishing_molecule_number, 'Value')};

if ( isempty(strfind(cur_molecule_str, '*')) )
    tstr{get(handles.mnu_fishing_molecule_number, 'Value')} = ['*', cur_molecule_str];
    set(handles.mnu_fishing_molecule_number, 'String', tstr);
    handles.saved_data(handles.curIndex).fishing_molecules_str = tstr;
    set(handles.btn_exclude_molecule, 'String', 'Include');
else
    tstr{get(handles.mnu_fishing_molecule_number, 'Value')} = strtok(cur_molecule_str,'*');
    set(handles.mnu_fishing_molecule_number, 'String', tstr);
    handles.saved_data(handles.curIndex).fishing_molecules_str = tstr;
    set(handles.btn_exclude_molecule, 'String', 'Exclude');
end
guidata(hObject, handles)

% --- Executes on button press in btn_exclude_preview.
function btn_exclude_preview_Callback(hObject, eventdata, handles)
% hObject    handle to btn_exclude_preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tstr = get(handles.mnu_fishing_molecule_number, 'String');

%plot non-ts events
axes(handles.axes_unextracted)
cla(handles.axes_unextracted)
hold on;
grid on;

for i=find(~strncmpi('*', tstr, 1))'

    cur_molecule = str2num(tstr{i});

    if( get(handles.chkbox_fishing_disp_indiv_mol, 'Value') == 1 && handles.fishing_data == 1)
        h1 = semilogx( handles.processed_data(intersect(find( handles.detected_events(:,7) == cur_molecule),  handles.events_wo_ts),6).*1000, ...
            handles.processed_data( intersect(find( handles.detected_events(:,7) == cur_molecule), handles.events_wo_ts),5),'.m','MarkerSize', 10);
    end
    set(h1, 'ButtonDownFcn', {@plot_1_Callback, hObject});
    if ( get(handles.chkbox_ue_wo_term_steps, 'Value') == 0)
        set(h1, 'Visible', 'off');
    end
    handles.ph_ue_evnts_wo_ts = h1;
    
    if( get(handles.chkbox_fishing_disp_indiv_mol, 'Value') == 1 && handles.fishing_data == 1)
        h2 = semilogx( handles.processed_data(intersect(find( handles.detected_events(:,7) == cur_molecule),  handles.events_with_ts),6).*1000, ...
            handles.processed_data( intersect(find( handles.detected_events(:,7) == cur_molecule), handles.events_with_ts),5),'.b','MarkerSize', 10);
    end
    set(h2, 'ButtonDownFcn', {@plot_1_Callback, hObject});
    if ( get(handles.chkbox_ue_w_term_steps, 'Value') == 0)
        set(h2,'Visible','off');
    end
    handles.ph_ue_evnts_with_ts = h2;

    if ( ~isempty(handles.rboxdata) && strcmp(handles.rboxloc, 'ue') )
        h3 = plot(handles.rboxdata(1,:),handles.rboxdata(2,:),'-r','linewidth',3); 
        handles.rbox = h3;
    end
    
    axis(handles.curaxis);
    xlabel('Duration (ms)','FontSize', 10);
    ylabel('Amplitude (pA)','FontSize', 10);
    
end  
%plot ts events
axes(handles.axes_extracted)
cla(handles.axes_extracted)
hold on;
grid on;
    
for i=find(~strncmpi('*', tstr, 1))'

    cur_molecule = str2num(tstr{i});

    h1 = semilogx( handles.processed_data(intersect(find( handles.detected_events(:,7) == cur_molecule),  handles.events_wo_ts),6).*1000, ...
        handles.processed_data( intersect(find( handles.detected_events(:,7) == cur_molecule), handles.events_wo_ts),5),'.m','MarkerSize', 10);

    set(h1, 'ButtonDownFcn', {@plot_2_Callback, hObject});
    if ( get(handles.chkbox_e_wo_term_steps, 'Value') == 0)
        set(h1, 'Visible', 'off');
    end
    handles.ph_e_evnts_wo_ts = h1;

    h2 = semilogx( handles.processed_data(intersect(find( handles.detected_events(:,7) == cur_molecule),  handles.events_with_ts),9).*1000, ...
        handles.processed_data( intersect(find( handles.detected_events(:,7) == cur_molecule), handles.events_with_ts),8),'.k','MarkerSize', 10);

    set(h2, 'ButtonDownFcn', {@plot_2_Callback, hObject});
    if( get(handles.chkbox_e_extract_pre_term_steps, 'Value') == 0)
        set(h2, 'Visible', 'off');
    end
    handles.ph_e_pre_term_steps = h2;

    h3 = semilogx( handles.processed_data(intersect(find( handles.detected_events(:,7) == cur_molecule),  handles.events_with_ts),3).*1000, ...
        handles.processed_data( intersect(find( handles.detected_events(:,7) == cur_molecule), handles.events_with_ts),2),'.c','MarkerSize', 10);

    set(h3, 'ButtonDownFcn', {@plot_2_Callback, hObject});
    if( get(handles.chkbox_e_extract_term_steps, 'Value') == 0)
        set(h3, 'Visible', 'off');
    end
    handles.ph_e_term_steps = h3;

    if ( ~isempty(handles.rboxdata) && strcmp(handles.rboxloc, 'e') )
        h3 = plot(handles.rboxdata(1,:),handles.rboxdata(2,:),'-r','linewidth',3);
        handles.rbox = h3;
    end
    axis(handles.curaxis);
    xlabel('Duration (ms)','FontSize', 10);
    ylabel('Amplitude (pA)','FontSize', 10);

end
