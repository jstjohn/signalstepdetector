function varargout = Waveform_Settings(varargin)
% WAVEFORM_SETTINGS M-file for Waveform_Settings.fig
%      WAVEFORM_SETTINGS, by itself, creates a new WAVEFORM_SETTINGS or raises the existing
%      singleton*.
%
%      H = WAVEFORM_SETTINGS returns the handle to a new WAVEFORM_SETTINGS or the handle to
%      the existing singleton*.
%
%      WAVEFORM_SETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WAVEFORM_SETTINGS.M with the given input arguments.
%
%      WAVEFORM_SETTINGS('Property','Value',...) creates a new WAVEFORM_SETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Waveform_Settings_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Waveform_Settings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Waveform_Settings

% Last Modified by GUIDE v2.5 01-Aug-2008 12:37:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Waveform_Settings_OpeningFcn, ...
                   'gui_OutputFcn',  @Waveform_Settings_OutputFcn, ...
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes just before Waveform_Settings is made visible.
function Waveform_Settings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Waveform_Settings (see VARARGIN)

% Choose default command line output for Waveform_Settings
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if(nargin > 3)
    for index = 1:1:(nargin-3),
        switch lower(varargin{index})
            case 'handles'
                set(handles.edt_alpha,'String',num2str(varargin{index+1}.alpha));
                set(handles.edt_thresh_amp,'String',num2str(varargin{index+1}.thresh_amp));
                set(handles.edt_baseline_window,'String',num2str(varargin{index+1}.baseline_window));
                %set(handles.edt_exp_descrip, 'String', num2str(varargin{index+1}.exp_descrip ) );
                %set(handles.edt_filename, 'String', num2str(varargin{index+1}.save_filename ) );
                break;
            otherwise
        end
    end
end


% UIWAIT makes Waveform_Settings wait for user response (see UIRESUME)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Outputs from this function are returned to the command line.
function varargout = Waveform_Settings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in popmnu_filtering.
function popmnu_filtering_Callback(hObject, eventdata, handles)
% hObject    handle to popmnu_filtering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popmnu_filtering contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popmnu_filtering



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function popmnu_filtering_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popmnu_filtering (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edt_alpha_Callback(hObject, eventdata, handles)
% hObject    handle to edt_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_alpha as text
%        str2double(get(hObject,'String')) returns contents of edt_alpha as a double



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function edt_alpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_alpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edt_thresh_amp_Callback(hObject, eventdata, handles)
% hObject    handle to edt_thresh_amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_thresh_amp as text
%        str2double(get(hObject,'String')) returns contents of edt_thresh_amp as a double



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function edt_thresh_amp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_thresh_amp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);
set(handles.figure1,'Visible','off');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edt_baseline_window_Callback(hObject, eventdata, handles)
% hObject    handle to edt_baseline_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_baseline_window as text
%        str2double(get(hObject,'String')) returns contents of edt_baseline_window as a double



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function edt_baseline_window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_baseline_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edt_exp_descrip_Callback(hObject, eventdata, handles)
% hObject    handle to edt_exp_descrip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_exp_descrip as text
%        str2double(get(hObject,'String')) returns contents of edt_exp_descrip as a double



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function edt_exp_voltage_Callback(hObject, eventdata, handles)
% hObject    handle to edt_exp_voltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_exp_voltage as text
%        str2double(get(hObject,'String')) returns contents of edt_exp_voltage as a double



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function edt_exp_voltage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_exp_voltage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on selection change in popmnu_data_type.
function popmnu_data_type_Callback(hObject, eventdata, handles)
% hObject    handle to popmnu_data_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popmnu_data_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popmnu_data_type



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes during object creation, after setting all properties.
function popmnu_data_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popmnu_data_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popmnu_baseline_type.
function popmnu_baseline_type_Callback(hObject, eventdata, handles)
% hObject    handle to popmnu_baseline_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popmnu_baseline_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popmnu_baseline_type

contents = get(hObject,'String');

switch contents{get(hObject,'Value')}
    case 'Fraction'
        set(handles.txt_baseline_fraction,'Visible','on');
        set(handles.txt_baseline_time, 'Visible', 'off');
    case 'Time'
        set(handles.txt_baseline_fraction,'Visible','off');
        set(handles.txt_baseline_time, 'Visible', 'on');
end

% --- Executes during object creation, after setting all properties.
function popmnu_baseline_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popmnu_baseline_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in chkbox_rem_trans.
function chkbox_rem_trans_Callback(hObject, eventdata, handles)
% hObject    handle to chkbox_rem_trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkbox_rem_trans


% --- Executes on button press in chkbox_event_trans.
function chkbox_event_trans_Callback(hObject, eventdata, handles)
% hObject    handle to chkbox_event_trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkbox_event_trans
if( get(hObject,'Value') )
    set(handles.edtTransAddSamples, 'Visible', 'on');
    set(handles.txtTransAdd, 'Visible', 'on');
else
    set(handles.edtTransAddSamples, 'Visible', 'off');
    set(handles.txtTransAdd, 'Visible', 'off');
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in chkbox_lim_cur.
function chkbox_lim_cur_Callback(hObject, eventdata, handles)
% hObject    handle to chkbox_lim_cur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkbox_lim_cur

if( get(hObject,'Value') )
    set(handles.edt_cur_min, 'Visible', 'on');
    set(handles.edt_cur_max, 'Visible', 'on');
    set(handles.txt_min, 'Visible', 'on');
    set(handles.txt_max, 'Visible', 'on');
else
    set(handles.edt_cur_min, 'Visible', 'off');
    set(handles.edt_cur_max, 'Visible', 'off');
    set(handles.txt_min, 'Visible', 'off');
    set(handles.txt_max, 'Visible', 'off');
end


function edt_cur_min_Callback(hObject, eventdata, handles)
% hObject    handle to edt_cur_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_cur_min as text
%        str2double(get(hObject,'String')) returns contents of edt_cur_min as a double


% --- Executes during object creation, after setting all properties.
function edt_cur_min_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_cur_min (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_cur_max_Callback(hObject, eventdata, handles)
% hObject    handle to edt_cur_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_cur_max as text
%        str2double(get(hObject,'String')) returns contents of edt_cur_max as a double


% --- Executes during object creation, after setting all properties.
function edt_cur_max_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_cur_max (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edt_filename_Callback(hObject, eventdata, handles)
% hObject    handle to edt_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_filename as text
%        str2double(get(hObject,'String')) returns contents of edt_filename as a double


% --- Executes during object creation, after setting all properties.
function edt_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtTransAddSamples_Callback(hObject, eventdata, handles)
% hObject    handle to edtTransAddSamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtTransAddSamples as text
%        str2double(get(hObject,'String')) returns contents of edtTransAddSamples as a double
handles.transAddSamples = get( handles.edtTransAddSamples, 'Value');

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edtTransAddSamples_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtTransAddSamples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%set( handles.edtTranAddSamples, 'String', num2str( handles.transAddSamples) );

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_exp_stn_Callback(hObject, eventdata, handles)
% hObject    handle to edt_exp_stn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_exp_stn as text
%        str2double(get(hObject,'String')) returns contents of edt_exp_stn as a double


% --- Executes during object creation, after setting all properties.
function edt_exp_stn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_exp_stn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


