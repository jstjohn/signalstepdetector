function varargout = FSMRampCreator(varargin)
% FSMRAMPCREATOR M-file for FSMRampCreator.fig
%      FSMRAMPCREATOR, by itself, creates a new FSMRAMPCREATOR or raises the existing
%      singleton*.
%
%      H = FSMRAMPCREATOR returns the handle to a new FSMRAMPCREATOR or the handle to
%      the existing singleton*.
%
%      FSMRAMPCREATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FSMRAMPCREATOR.M with the given input arguments.
%
%      FSMRAMPCREATOR('Property','Value',...) creates a new FSMRAMPCREATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FSMRampCreator_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FSMRampCreator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help FSMRampCreator

% Last Modified by GUIDE v2.5 19-Sep-2008 18:16:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FSMRampCreator_OpeningFcn, ...
                   'gui_OutputFcn',  @FSMRampCreator_OutputFcn, ...
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


% --- Executes just before FSMRampCreator is made visible.
function FSMRampCreator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FSMRampCreator (see VARARGIN)

% Choose default command line output for FSMRampCreator
handles.output = hObject;

%Variables
handles.samplePeriodOutput = (5.3e-6/4);
handles.samplePeriodInput = 5.3e-6;
handles.bits_to_output_volts = 20e-3*(10/32768);
handles.step_size = 1;
handles.step_length = 1;
handles.start_voltage = -50e-3;
handles.end_voltage = 150e-3;
handles.time = 1000;

guidata(hObject, handles);

pause(0.01)

display_step(hObject, handles)

handles = guidata(hObject);

calculate_values(hObject,handles)

handles = guidata(hObject);

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes FSMRampCreator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FSMRampCreator_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edt_ramp_rate_Callback(hObject, eventdata, handles)
% hObject    handle to edt_ramp_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_ramp_rate as text
%        str2double(get(hObject,'String')) returns contents of edt_ramp_rate as a double


% --- Executes during object creation, after setting all properties.
function edt_ramp_rate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_ramp_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_time_Callback(hObject, eventdata, handles)
% hObject    handle to edt_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_time as text
%        str2double(get(hObject,'String')) returns contents of edt_time as a double

handles.time = round(str2num(get(handles.edt_time,'String')));

if ( handles.time < 1)
    handles.time = 1;
end

set(handles.edt_time, 'String', num2str(handles.time))

display_step(hObject, handles)

calculate_values(hObject,handles)

handles = guidata(hObject);

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edt_time_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_step_size_Callback(hObject, eventdata, handles)
% hObject    handle to edt_step_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_step_size as text
%        str2double(get(hObject,'String')) returns contents of edt_step_size as a double

handles.step_size = round(str2num(get(handles.edt_step_size,'String')));

if ( handles.step_size < 1)
    handles.step_size = 1;
end

set(handles.edt_step_size, 'String', num2str(handles.step_size))

display_step(hObject, handles)

calculate_values(hObject,handles)

handles = guidata(hObject);

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edt_step_size_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_step_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_step_length_Callback(hObject, eventdata, handles)
% hObject    handle to edt_step_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_step_length as text
%        str2double(get(hObject,'String')) returns contents of edt_step_length as a double

handles.step_length = round(str2num(get(handles.edt_step_length,'String')));

if ( handles.step_length < 1)
    handles.step_length = 1;
end

set(handles.edt_step_length, 'String', num2str(handles.step_length))

display_step(hObject, handles)

calculate_values(hObject,handles)

handles = guidata(hObject);

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edt_step_length_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_step_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_step_length_actual_Callback(hObject, eventdata, handles)
% hObject    handle to edt_step_length_actual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_step_length_actual as text
%        str2double(get(hObject,'String')) returns contents of edt_step_length_actual as a double


% --- Executes during object creation, after setting all properties.
function edt_step_length_actual_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_step_length_actual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edt_step_size_actual_Callback(hObject, eventdata, handles)
% hObject    handle to edt_step_size_actual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_step_size_actual as text
%        str2double(get(hObject,'String')) returns contents of edt_step_size_actual as a double


% --- Executes during object creation, after setting all properties.
function edt_step_size_actual_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_step_size_actual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in chkbox_time.
function chkbox_time_Callback(hObject, eventdata, handles)
% hObject    handle to chkbox_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkbox_time


% --- Executes on button press in chkbox_ramp_rate.
function chkbox_ramp_rate_Callback(hObject, eventdata, handles)
% hObject    handle to chkbox_ramp_rate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkbox_ramp_rate


% --- Executes on button press in chkbox_step_length.
function chkbox_step_length_Callback(hObject, eventdata, handles)
% hObject    handle to chkbox_step_length (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkbox_step_length


% --- Executes on button press in chkbox_step_size.
function chkbox_step_size_Callback(hObject, eventdata, handles)
% hObject    handle to chkbox_step_size (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of chkbox_step_size


% --- Executes on button press in btn_open_fsm.
function btn_open_fsm_Callback(hObject, eventdata, handles)
% hObject    handle to btn_open_fsm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edt_start_mv_Callback(hObject, eventdata, handles)
% hObject    handle to edt_start_mv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_start_mv as text
%        str2double(get(hObject,'String')) returns contents of edt_start_mv as a double

handles.start_voltage = str2num(get(handles.edt_start_mv,'String'))*1e-3;

display_step(hObject, handles)

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edt_start_mv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_start_mv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function edt_end_mv_Callback(hObject, eventdata, handles)
% hObject    handle to edt_end_mv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edt_end_mv as text
%        str2double(get(hObject,'String')) returns contents of edt_end_mv as a double

handles.end_voltage = str2num(get(handles.edt_end_mv,'String'))*1e-3;

display_step(hObject, handles)

guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edt_end_mv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edt_end_mv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%Responsible for calculating the values
function calculate_values(hObject,handles)

ramp_rate = (handles.step_size*handles.bits_to_output_volts)/(handles.step_length*handles.samplePeriodOutput);

set(handles.edt_ramp_rate, 'String', num2str(ramp_rate))
set(handles.txt_step_length_actual, 'String', handles.step_length*(handles.samplePeriodOutput)*1e3);
set(handles.txt_step_size_actual, 'String', handles.step_size*handles.bits_to_output_volts*1e3);
set(handles.txt_time_actual, 'String', handles.time*handles.samplePeriodInput*1e3);

guidata(hObject, handles);

%Plots a sample ramp with the given parameters
function display_step(hObject, handles)

%Sample Period is 5.3us/4 == 1.325us
%We'll show 2ms before and after ramp

total_cycles_output = 3000+handles.time*4;

plot_time = 1:total_cycles_output;

plot_time = plot_time.*handles.samplePeriodOutput;
plot_voltage = ones(1500,1)*handles.start_voltage;

ramp_index = 1;
ramp_voltage = [];
next_voltage = handles.start_voltage;

while ramp_index < handles.time*4
    ramp_voltage = [ramp_voltage; ones(handles.step_length,1).*next_voltage];
    next_voltage = next_voltage + handles.step_size*handles.bits_to_output_volts;
    ramp_index = ramp_index + handles.step_length;
end

if ( length(ramp_voltage) < handles.time*4)
    ramp_voltage = [ramp_voltage; ones((handles.time*4-length(ramp_voltage)),1)*ramp_voltage(end)];
else
    ramp_voltage = ramp_voltage(1:handles.time*4);
end

set(handles.txt_error_mv,'String', num2str(((ramp_voltage(end))-handles.end_voltage)*1e3))
disp(ramp_voltage(end)*1e3)

plot_voltage = [plot_voltage; ramp_voltage; ones(1500,1)*handles.end_voltage];

axes(handles.axes_sample_trace)
plot(plot_time, plot_voltage);

temp_axis = axis;

if ( get(handles.tglbtn_zoom,'Value') == 0 )
    axis([temp_axis(1)-0.2e-3 temp_axis(2)+0.2e-3 temp_axis(3)-10e-3 temp_axis(4)+10e-3])
else
    axis([plot_time(1500+handles.time*4-50) plot_time(1500+handles.time*4+50) plot_voltage(1500+handles.time*4-50) max(plot_voltage(1500+handles.time*4-50:1500+handles.time*4+50)+1e-3)])
end


% --- Executes on button press in tglbtn_zoom.
function tglbtn_zoom_Callback(hObject, eventdata, handles)
% hObject    handle to tglbtn_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tglbtn_zoom

display_step(hObject, handles)



% --- Executes on button press in btn_guess.
function btn_guess_Callback(hObject, eventdata, handles)
% hObject    handle to btn_guess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

voltage_change = (handles.end_voltage - handles.start_voltage);

approx_rate = voltage_change/(handles.time*handles.samplePeriodInput);

%ramp_rate = (handles.step_size*handles.bits_to_output_volts)/(handles.step_length*handles.samplePeriod);

approx_step_size = floor((approx_rate*(handles.step_length*handles.samplePeriodOutput))/handles.bits_to_output_volts)

handles.step_size = approx_step_size

set(handles.edt_step_size, 'String', num2str(handles.step_size));

guidata(hObject, handles)

display_step(hObject, handles)

handles = guidata(hObject);

notdone = 1;

preverr = str2num(get(handles.txt_error_mv,'String'));

while notdone
    handles.time = handles.time+1;
    set(handles.edt_time, 'String', num2str(handles.time));
    guidata(hObject, handles)
    display_step(hObject, handles)
    handles = guidata(hObject);
    curerr = str2num(get(handles.txt_error_mv,'String'))
    if curerr > 0
        if abs(curerr) > abs(preverr)
            handles.time = handles.time - 1;
            set(handles.edt_time, 'String', num2str(handles.time))
        end
        notdone = 0;
    end
    preverr = curerr;
end

guidata(hObject, handles)

display_step(hObject, handles)

handles = guidata(hObject);

guidata(hObject, handles)
