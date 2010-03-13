%% get_event
% --- This will return the event waveform specified in the call
%     It will truncate the event to the appropriate length
function [data, time, voltage] = get_event(handles, hObject)
% handles          handles of GUI
% event            the event of interest, which is the index of
%                  handles.event_string

% disp('---------------')
% disp([handles.pathname handles.filename])
% disp(handles.detected_events(event,1))
% disp('---------------')

voltage = [];

tempEventString = handles.event_string{handles.curindex};
selectedSweepInfo = sscanf(tempEventString, '%d.%d');
handles.matIndex = selectedSweepInfo(1);
handles.event = selectedSweepInfo(2);
detected_events = handles.detected_events{handles.matIndex};
handles.cur_time_tick = handles.time_tick{handles.matIndex};

tempEventType = detected_events(handles.event,8);
tempDNANum = detected_events(handles.event,7);
% fprintf('Event Type: %d\n DNA Number: %d\n', tempEventType, tempDNANum);
set(handles.txt_event_type, 'String', num2str(tempEventType));


try
    beginning_points_to_trim = round(str2num(get(handles.edt_beginning_trim, 'String')));
catch
    beginning_points_to_trim = 0;
    disp(lasterr)
end

try
    end_points_to_trim = round(str2num(get(handles.edt_end_trim, 'String')));
catch
    end_points_to_trim = 0;
    disp(lasterr)
end

if( detected_events(handles.event,10) > 0 )
    if (~strcmp(handles.filename,handles.allFilenames{handles.matIndex}{detected_events(handles.event,10)}))
        handles.cursweep = [];
    end
    handles.pathname = handles.allPathnames{handles.matIndex}{detected_events(handles.event,10)};
    handles.filename = handles.allFilenames{handles.matIndex}{detected_events(handles.event,10)};
end

start_sweep = detected_events(handles.event,4);
end_sweep = detected_events(handles.event,5);
%     cur_sweep = start_sweep;

sweeps_to_load = start_sweep:1:end_sweep;
sweeps_to_del_cur = [];
sweeps_to_del_to_load = [];

%handles.cursweep contains the sweeps currently loaded
%cursweep indexes trace_data

%first find out which sweeps we can delete trace_data
for i=1:length(handles.cursweep)
    searchIndex = find(handles.cursweep(i) == sweeps_to_load);

    if isempty(searchIndex)
        sweeps_to_del_cur = [sweeps_to_del_cur;i];
    else
        sweeps_to_del_to_load = [sweeps_to_del_to_load; searchIndex];
    end

end

handles.trace_data(sweeps_to_del_cur,:) = [];

if( ~isempty(handles.voltage_data) )
    handles.voltage_data(sweeps_to_del_cur,:) = [];
end
handles.cursweep(sweeps_to_del_cur) = [];
sweeps_to_load(sweeps_to_del_to_load) = [];

%handles.trace_data will contain the loaded traces for the current
%sweeps

if ( isempty(handles.cursweep) )
    for i=1:length(sweeps_to_load)
        [Trace_Data, SampleInt] = import_abf([handles.pathname handles.filename],sweeps_to_load(i));
        handles.SampleInt = SampleInt;
        handles.trace_data(i,:) = Trace_Data(:,1)';
        if ( size(Trace_Data,2) == 2)
            handles.voltage_data(i,:) = Trace_Data(:,2)';
        end
    end
    handles.cursweep = sweeps_to_load;
elseif ( isempty(sweeps_to_load) )
    %do nothing
else
    %Take out the sweeps that we have from sweeps_to_load

    if (sweeps_to_load(1) < handles.cursweep(1))
        %Put them at the beginning of the trace_data
        for i=length(sweeps_to_load):-1:1
            [Trace_Data, SampleInt] = import_abf([handles.pathname handles.filename],sweeps_to_load(i));
            handles.trace_data = [Trace_Data(:,1)'; handles.trace_data];
            if ( size(Trace_Data,2) == 2)
                handles.voltage_data = [Trace_Data(:,2)'; handles.voltage_data];
            end
        end
        handles.cursweep = [sweeps_to_load handles.cursweep];
    else
        %put them at the end
        for i=1:1:length(sweeps_to_load)
            [Trace_Data, SampleInt] = import_abf([handles.pathname handles.filename],sweeps_to_load(i));
            handles.trace_data = [handles.trace_data; Trace_Data(:,1)'];
            if ( size(Trace_Data,2) == 2)
                handles.voltage_data = [handles.voltage_data; Trace_Data(:,2)'];
            end
        end
        handles.cursweep = [handles.cursweep sweeps_to_load];
    end

end

switch( end_sweep-start_sweep )
    case 0
        %The event is in only one sweep
        if ( handles.filtered == 1 )
            %add on a bit to the front so it can filter
            if ( detected_events(handles.event,2) > round(1000/handles.SampleInt) + 1 )
                data = handles.trace_data(1,(detected_events(handles.event,2)-round(1000/handles.SampleInt)):detected_events(handles.event,3));
            else
                data = [handles.trace_data(1,detected_event(handles.event,2))*ones(1, round(1000/handles.SampleInt)); handles.trace_data(1,detected_events(handles.event,2):detected_events(handles.event,3))];
            end
        else
            data = handles.trace_data(1,detected_events(handles.event,2):detected_events(handles.event,3));
        end
        if( ~isempty(handles.voltage_data) )
            voltage = handles.voltage_data(1,detected_events(handles.event,2):detected_events(handles.event,3));
        end
    otherwise
        %The event spans multiple sweeps
        %Get first chunck

        if ( handles.filtered == 1 )
            %add on a bit to the front so it can filter
            if ( detected_events(handles.event,2) > round(1000/handles.SampleInt) + 1 )
                data = handles.trace_data(1,(detected_events(handles.event,2)-round(1000/handles.SampleInt)):end);
            else
                data = [handles.trace_data(1,detected_event(handles.event,2))*ones(1, round(1000/handles.SampleInt)); handles.trace_data(1,detected_events(handles.event,2):end)];
            end
        else
            data = handles.trace_data(1,detected_events(handles.event,2):end);
        end
        if( ~isempty(handles.voltage_data))
            %                 voltage = Trace_Data(detected_events(handles.event,2):end,2);
            voltage = handles.voltage_data(1,detected_events(handles.event,2):end);
        end
        %Get Rest of chunk until end

        for i = 2:length(handles.cursweep)-1
            data = [data handles.trace_data(i,:)];
            if( ~isempty(handles.voltage_data) )
                voltage = [voltage handles.voltage_data(i,:)];
            end
        end
        %Finish it off

        data = [data handles.trace_data(end,1:detected_events(handles.event,3))];

        if( ~isempty(handles.voltage_data) )
            voltage = [voltage handles.voltage_data(end,1:detected_events(handles.event,3))];
        end

end

if ( handles.filtered == 1 )
    data = filter(handles.Hb, data);
    data = data((round(1000/handles.SampleInt)+1):end)';
else
    data = data';
end

if ( handles.get_raw_data == 1 )
    voltage = voltage';
    handles.SampleInt = handles.cur_time_tick*1e6;
    time = 0:(handles.SampleInt*1e-6):(length(data)-1)*handles.SampleInt*1e-6;
    guidata(hObject,handles);
    return;
end
    
%load the frequency in the text box if advanced options isn't set
set(handles.edt_sample_freq, 'String', handles.time_tick_options{handles.matIndex});

%add compensation exponential if coefs present
if( handles.expCoefs_computed{handles.matIndex}(detected_events(handles.event,7)+1) == 1 && get(handles.chkbox_trans_remove_enable, 'Value') == 1)

    set(handles.txt_have_coeff, 'BackgroundColor', [0 1 0]);
    %removing transient
    %trim off amp saturation
%     data = data(80:end);
    time = (0:(handles.cur_time_tick):(length(data)-1)*handles.cur_time_tick)';       
    
    %rebuild exponential
    if( time(end-handles.saturatedSampleTime{handles.matIndex}(detected_events(handles.event,7)+1)) < handles.maxTransientLength{handles.matIndex}(detected_events(handles.event,7)+1) )
        endSample = length(time(1:end-handles.saturatedSampleTime{handles.matIndex}(detected_events(handles.event,7)+1)));
    else
        endSample = floor( handles.maxTransientLength{handles.matIndex}(detected_events(handles.event,7)+1)/(handles.cur_time_tick) );
    end

    %compute exponetial signal
    %
    % W. Dunbar - Changes
    %   Feb. 23, 2009. Changed "handles.expCoefs(event,:)" to
    %   "handles.expCoefs(1,:)", since only coefs for first event
    %   fitted used in the compensation 
    expFit = exponentialFit(handles.expCoefs{handles.matIndex}(handles.detected_events{handles.matIndex}(handles.event,7)+1,:), time(1:endSample));

    %add to data to remove transient.
    %
    % W. Dunbar - Changes
    %   Feb. 23, 2009. Removed IF conditions by commenting them out. Also
    %   changed sign from " ... + (expFit - expFit(end))" to " ... -
    %   (expFit - expFit(end))" since fit is to probing transient (based on
    %   transientComp2.m file).
    %if( abs( expFit(1) ) < 5 )
        data(handles.saturatedSampleTime{handles.matIndex}(detected_events(handles.event,7)+1):endSample+handles.saturatedSampleTime{handles.matIndex}(detected_events(handles.event,7)+1)-1) ...
            = data(handles.saturatedSampleTime{handles.matIndex}(detected_events(handles.event,7)+1):endSample+handles.saturatedSampleTime{handles.matIndex}(detected_events(handles.event,7)+1)-1) - (expFit - expFit(end));
%     else
%         data(handles.saturatedSampleTime:endSample+handles.saturatedSampleTime-1) ...
%             = data(handles.saturatedSampleTime:endSample+handles.saturatedSampleTime-1) + expFit;
%     end
    
    voltage = voltage';

else
    %set(handles.txt_have_coeff, 'BackgroundColor', [1 0 0]);
    voltage = voltage';
    time = 0:(handles.cur_time_tick):(length(data)-1)*handles.cur_time_tick;
end

%Trim off the beginning steady state for ramping events
if ( handles.ramping_file{handles.matIndex} == 1 && tempEventType == 5)
    try
        first2ms_max = max(voltage(1:round(2e-3/handles.cur_time_tick),1));
    catch
        first2ms_max = max(voltage,1);
    end
    start_ramp_index = find(voltage > first2ms_max, 1);
    if( ~isempty(start_ramp_index) )
        data = data(start_ramp_index:end);
        voltage = voltage(start_ramp_index:end);
        time = 0:(handles.cur_time_tick):(length(data)-1)*handles.cur_time_tick;
    end
end

if ( get(handles.chkbox_adv_options, 'Value') == 0 )
    set(handles.edt_sample_freq, 'Value', 1)
    handles.SampleInt = handles.cur_time_tick*1e6;
else
    handles.cur_time_tick = handles.cur_time_tick*get(handles.edt_sample_freq,'Value');
    data = data(1:get(handles.edt_sample_freq,'Value'):end);
    voltage = voltage(1:get(handles.edt_sample_freq,'Value'):end);
    time = time(1:get(handles.edt_sample_freq, 'Value'):end);
    handles.SampleInt = handles.cur_time_tick*1e6;
end

if (length(data) <= beginning_points_to_trim+end_points_to_trim+1)
    if ( beginning_points_to_trim == 0 )
        end_points_to_trim = length(data)-3;
    elseif( end_points_to_trim == 0)
        beginning_points_to_trim = length(data)-3;
    else
        beg_trim_ratio = beginning_points_to_trim/(end_points_to_trim+beginning_points_to_trim);
        end_trim_ratio = end_points_to_trim/(end_points_to_trim+beginning_points_to_trim);
        amount_to_trim = length(data)-3;%+end_points_to_trim+beginning_points_to_trim
        beginning_points_to_trim = floor(beg_trim_ratio*amount_to_trim);
        end_points_to_trim = floor(end_trim_ratio*amount_to_trim);
    end
end

data = data(1+beginning_points_to_trim:end-end_points_to_trim);
voltage = voltage(1+beginning_points_to_trim:end-end_points_to_trim);
time = time(1:end-(end_points_to_trim+beginning_points_to_trim));
guidata(hObject, handles)
