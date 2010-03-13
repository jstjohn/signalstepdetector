function [processedSweep, processedSteps] = ChisquaredFit_Custom(hObject, handles, plot_yn)
%this function adresses a two-dim array 'data'in a specific segment
%and determines the best step-fit there
%data: data
%i1: start ind
%i2: end ind
%As far as I can tell this will just be the whole data segment
term_reached = [];
term_amp = [];
term_dur = [];
term_iqr =[];
wv_amp = [];
wv_dur = [];
wv_iqr = [];
pre_term_amp = [];
pre_term_dur = [];
pre_term_iqr = [];
processedSteps = [];
term_time = 0;
multi_steps_amp = {};
multi_steps_dur = {};

% Get the data for the event
[trace_data, time, voltage] = get_event(handles, hObject);


handles = guidata(hObject);

%preform a double check step. This may actually be a detrement to
%performance?
if(get(handles.edt_check_output,'Value') == get(handles.edt_check_output,'Max'))
    double_check = true;
else
    double_check = false;
end

%see if the user wants an attempt at labeling states.
if (get(handles.edt_label_steps,'Value') == get(handles.edt_label_steps,'Max'))
    label_states = true;
else
    label_states = false;
end


window_size= str2num(get(handles.edtWinLen,'String'));
var = str2num(get(handles.edtVar,'String'));
F = @(xmean,ymean,xlen,ylen) ((xmean - ymean)^2)/(((1/xlen)+(1/ylen))*var);
%erfc((x^2)

%now to calculate an appropriate cutoff given our window size, (but only for the "double" cases)
%we know both window lengths, and we need to find out the appropriate
%amplitude difference cutoff value such that p < 0.05 to conlcude that the
%shift in mean amplitudes was not by chance.

%Solving for this in Wolfram-Alpha solve for M, erfc((M^2)/((2*V)/L)) =
%0.05. I get M = +/- (1.66487 * sqrt(var)/sqrt(L)) Where M is the
%difference in means, L is the window length, and V is the variance.
pval = str2num(get(handles.edtPVal,'String'));
query = strcat('erfc(sqrt(((M^2)/((2*',num2str(var),')/',num2str(window_size),'))/2))=',num2str(pval));
%minStepSize = 3.64277 * sqrt(var)/sqrt(window_size);
%inStepSize = solve();
%disp('minStepSize '+minStepSize)
result=solve(query);
%the result will be in the form number, -number, and type 'sym'. First we
%take the absolute value of the number, and then we cast it to a double so
%that we can use in as a window mean difference cutoff value.
minStepSize=abs(double(result(1)));
try
    exp_description = char(get(handles.edt_exp_descrip,'String'));
catch
    exp_description = '';
end
try
    exp_voltage = char(get(handles.edt_exp_voltage,'String'));
catch
    exp_voltage = '';
end
try
    exp_station = char(get(handles.edt_exp_stn,'String'));
catch
    exp_station = '';
end

if handles.curindex == 1
    % Open the file and append these results to the end
    fid = fopen(strcat(get(handles.edt_save_filename,'String'),'.txt'),'w');
    fprintf(fid, strcat('#Experiment_Description:',exp_description, '\n'));
    fprintf(fid, strcat('#Experiment_Voltage:',exp_voltage, '\n'));
    fprintf(fid, strcat('#Experiment_Station:',exp_station,'\n'));
    fprintf(fid, '#p_value_cutoff:%e\n',pval);
    fprintf(fid, '#window_size:%d\n',window_size);
    fprintf(fid, '#minimum_step_size:%f\n',minStepSize);
    fprintf(fid, '#baseline_variance:%f\n',var);
    fprintf(fid, strcat('#file_name:',handles.filename,'\n\n'));
    if label_states
        fprintf(fid, 'Start_Time\tMean_Amplitude\tStep_Time\tState_Label_By_Order_Seen_In_Sweep\n');
    else
        fprintf(fid, 'Start_Time\tMean_Amplitude\tStep_Time\n');
    end
    fclose(fid);
    
end


alpha = handles.alpha;
%filt = exp_filt(trace_data,1,length(trace_data),alpha);
filt = trace_data;
filt(:,2) = time;
% Get the parameters of the step checking from waveformviewer
minStepLength = str2num(get(handles.edtWinLen, 'String'))*1e-3;
%minStepSize = str2num(get(handles.edtMinStepSize,'String'));
%eject_working = get(handles.chkbox_primer_eject_working, 'Value');
eject_working = 0;


if 0
    if( handles.detected_events{handles.matIndex}(handles.event,8) == 2 && length(trace_data)*handles.cur_time_tick > .24)
        processedSweep = mean(trace_data(end-round(10e-3/handles.cur_time_tick):end))
    else
        processedSweep = [];
        disp(['nope:' num2str(length(trace_data)*handles.cur_time_tick)])
    end
    processedSteps = [];
    return
end
% If the file is ramping and we want to treat it as such
if (handles.ramping_file{handles.matIndex} == 1 && strcmp(handles.data_type, 'ramping'))
    
    %     if( length(trace_data) > length(handles.combined_ramps)  )
    %         difference =trace_data(1:length(handles.combined_ramps))./handles.combined_ramps;
    %         time = time(1:length(handles.combined_ramps));
    %     else
    %         difference = trace_data./handles.combined_ramps(1:length(trace_data));
    %     end
    
    %trying to find where the voltage ramp stops
    try
        last10ms_min = min(voltage(end-round(5e-3/handles.cur_time_tick):end-round(.1e-3/handles.cur_time_tick)));
        kink_point = find(voltage(1:end-round(5e-3/handles.cur_time_tick)) < last10ms_min,1,'last');
    catch
        %         rethrow(lasterror)
        last10ms_min = mean(voltage);
        kink_point = find(voltage < last10ms_min, 1, 'last');
    end
    
    % Get the conductance, will run step detection on this, thanks bill!
    difference = trace_data./voltage;
    
    %     trace_data = filter(handles.b,handles.a,difference);
    trace_data = difference;
    filt = trace_data;
else
    after_kink = -1;
end

if( handles.primer_file{handles.matIndex} == 1 )
    if( handles.detected_events{handles.matIndex}(handles.event,8) == 11 || handles.detected_events{handles.matIndex}(handles.event,8) == 10 )
        %We're just going to use our threshold not going to bother with step
        %detection. It is a annealed primer if two things happen. 1) There was
        %an eject event, or 2) the current falls below our threshold at some
        %point after our trimming line. we'll keep track of a few things,
        %whether there was an association event, was the dna in pore before
        %probing? if the primer dissocated before we ejected it how long was it
        %around. we'll alos keep track of how long the detection took for the
        %the fsm to recgonize it as a dna
        
        %ProcessedSweep will be different then usual
        %[-1-DNA wasn't in pore @ end of hold but saw annealing, 0-DNA wasn't in pore at end of hold, 1-DNA was in pore, but no
        %primer anneal event, 2-DNA in pore, with primer anneal, primer was
        %ejected, 3-DNA in pore, w/ anneal, dissocated on own;
        %
        %Total time DNA was held before probe;
        %
        %Detection time-Won't be able to tell if probe and capture voltage
        %are same;
        %
        %Dissocation Time - If dissociated before eject, how long]
        
        hold_ind = handles.detected_events_voltages{handles.matIndex}(handles.event,1);
        probe_ind = handles.detected_events_voltages{handles.matIndex}(handles.event,2);
        eject_ind = handles.detected_events_voltages{handles.matIndex}(handles.event,3);
        probe_trans_trim = str2num(get(handles.edt_primer_transient_trim,'String'));
        dna_hold_current = str2num(get(handles.edt_primer_hold_current, 'String'));
        probe_threshold_current = str2num(get(handles.edt_primer_open_channel_probe, 'String'));
        
        if( eject_ind == 0 )
            %There wasn't an eject event so were going to look at the rest of
            %the event that was past the probing voltage
            primer_ejected = 0
            eject_ind = length(trace_data);
        else
            primer_ejected = 1
        end
        
        %we'll look at the last 10 samples and if any are above the hold
        %voltage + 10pA then we'll call it a miss held
        if(~isempty(find(trace_data(probe_ind-10:probe_ind) > dna_hold_current+10)))
            %DNA wasn't in pore when change was made
            dna_in_pore = 0
        else
            dna_in_pore = 1
        end
        
        %If it wasn't an primer eject event then we can look to see if the
        %current went below the current threshsold for the holding current
        dissociate_ind = find(trace_data(probe_ind + probe_trans_trim:eject_ind) < probe_threshold_current,1);
        disp(trace_data(eject_ind))
        dissociate_ind_end = trace_data(eject_ind) < probe_threshold_current
        if( primer_ejected == 0 &&  eject_working == 1 )
            if( isempty(dissociate_ind) )
                %Could find dissociation, this is bad since it wasn't eject
                %event
                disp('PRIMER: Couldn''t find dissociation when it wasn''t an eject event')
                dissociate_ind = length(trace_data(probe_ind + probe_trans_trim:eject_ind))-1;
                primer_anneal = 1
            elseif( dissociate_ind == 1)
                %It was already dissociated or never associated when it was
                %probed
                primer_anneal = 0
            else
                %Found where it dissociated.
                primer_anneal = 1
            end
        elseif( primer_ejected == 1 && eject_working == 1)
            primer_anneal = 1
        else
            %Eject is meaningless, will eject no matter what
            if( isempty(dissociate_ind) )
                primer_anneal = 1;
            elseif( dissociate_ind == 1)
                %It was already dissociated or never associated when it was
                %probed
                primer_anneal = 0
            else
                %Found where it dissociated.
                primer_anneal = 1
            end
        end
        
        %Now we'll look at all of the possiblities for our event and
        %compile the data
        dissociation_time = 0;
        
        if( eject_working == 1)
            if( dna_in_pore == 0 && (primer_ejected == 1 || primer_anneal == 1) )
                %Something has gone horribly wrong, basically a false positive
                detection_state = -1;
            elseif( dna_in_pore == 0 )
                %DNA wasn't in pore and didn't see annealing
                detection_state = 0;
            elseif( dna_in_pore == 1 )
                if( primer_ejected == 0 && primer_anneal == 0 )
                    %No observed association
                    detection_state = 1;
                elseif( primer_ejected == 1 )
                    %Primer was ejected, so it was associated
                    detection_state = 2;
                    %Calculate Time To Dissociaiton
                    dissociation_time = eject_ind - probe_ind;
                elseif( primer_ejected == 0 && primer_anneal == 1 )
                    %Association with some found dissociation event
                    detection_state = 3;
                    %Calculate Time To Dissociaiton
                    dissociation_time = dissociate_ind + probe_trans_trim;
                end
            end
        else
            if( dna_in_pore == 0 && primer_anneal == 1)
                %False detect
                detection_state = -1;
            elseif( dna_in_pore == 0 && primer_anneal == 0)
                %No annealing
                detection_state = 0;
            elseif( dna_in_pore == 1 && primer_anneal == 0)
                detection_state = 1;
            elseif( dna_in_pore == 1 && primer_anneal == 1)
                if( dissociate_ind_end == 0 )
                    detection_state = 2;
                    dissociation_time = eject_ind - probe_ind;
                else
                    detection_state = 3;
                    dissociation_time = dissociate_ind + probe_trans_trim;
                end
            end
        end
        
        dissociation_time = dissociation_time*handles.cur_time_tick;
        total_dna_hold_time = probe_ind*handles.cur_time_tick;
        dna_detection_time = hold_ind*handles.cur_time_tick;
        
        if strcmpi(plot_yn,'yes')
            if( detection_state == 3 )
                axes(handles.sweep_plot);
                line([time(probe_ind + probe_trans_trim + dissociate_ind) time(probe_ind + probe_trans_trim + dissociate_ind)], ylim(), 'color', 'c', 'linestyle', '-');
            end
        end
        processedSweep = [detection_state, total_dna_hold_time, dna_detection_time, dissociation_time];
        
    else
        processedSweep = [];
    end
    
    return
    
end

% % For fishing files
% if (handles.fishing_file{handles.matIndex} == 1)
%
%     % Most likely will trim some points from the beginning of the event to
%     % help in analysis
%     try
%         beginning_points_to_trim = round(str2num(get(handles.edt_beginning_trim, 'String')));
%     catch
%         beginning_points_to_trim = 0;
%         disp(lasterr)
%     end
%     wv_dur = time(end) + beginning_points_to_trim * handles.SampleInt*1e-6;
% else
%     wv_dur = time(end);
% end

% Get some general stats on the event
wv_amp = mean(trace_data);
wv_iqr = iqr(trace_data);
n = length(filt);

% We will only look for steps that have a minimum length, meaning we don't
% have to start looking for steps until that minimum step length has pass
% and we can stop 'early'

%start_point = round(minStepLength/handles.cur_time_tick);
%end_point = n - start_point -1;
start_point = 1;
end_point = n;

% This will be a matrix of the data points where there are steps, we will
% be looking in between each of the plautes generated by the steps for
% multiple steps, so to beginning we have a 'step' at the very first point
% and one at the end
%
% The first column will be which data point the step is at the second
% column is the chisquared value of the step and the 3rd column is the
% length ranked value of the step
step_matrix = [start_point 0.0 0];

%Chisq=zeros(length(filt)-3,1); %vector of zeros of the length of the data
%this will be the fit of each of the platueas
% Decide if we should analyze the event for steps or not
ignore_eject = get(handles.chkbox_ignore_eject, 'Value');
max_event_length = str2num(get(handles.edt_max_event_length, 'String'));
number_of_steps = 1; %initialize to one step seen.




% if ( ( (get(handles.chkbox_analyze_probing_events,'Value') == 1 && (handles.detected_events{handles.matIndex}(handles.event,8) == 3)) || get(handles.chkbox_analyze_probing_events,'Value') == 0) &&...
%         ( ignore_eject == 0 || (ignore_eject == 1 && handles.detected_events{handles.matIndex}(handles.event,8) == 2 && max_event_length > wv_dur) || (ignore_eject == 1 && ~strcmp(handles.data_type, 'fishing')  ) )&& ...
%         ( ignore_eject == 0 || (ignore_eject == 1 && max_event_length > wv_dur && strcmp(handles.data_type,'binary_ternary') ) || (ignore_eject == 1 && ~strcmp(handles.data_type, 'binary_ternary') )  ) )
if (end_point -start_point > window_size*2+2)
    %Will store the results for the optimal steps within each of the
    %plateaus, since the step matrix contains the start and end points
    %for each of the platues it will have the number of plateaus + 1
    %points, so to store the optimal steps for each of the plateaus we
    %will only need size(step_matrix) - 1
    %disp(['Step Num: ' num2str(i)])
    event_data = filt(start_point:end_point,1);
    
    window1_sum = sum(event_data(start_point:(start_point+window_size)));
    i = start_point+window_size;
    window2_sum = sum(event_data((i):(i+window_size)));
    sumMinStepSize = minStepSize * window_size;
    while i < (end_point-window_size-window_size-2)
        i = i + 1;
        %simulate advancing the window forward one position...
        %window_mean = window_mean + ((event_data(i-window_size-1)-event_data(i))/window_size);
        window1_sum = window1_sum - (event_data(i-window_size))+ (event_data(i));
        window2_sum = window2_sum - (event_data(i+1)) + (event_data(i+window_size+1));
        %Handle the end case!
        %
        if i >= end_point-window_size-window_size-3
            step_matrix(end+1,1) = end_point; %last step ends at the end point...
            number_of_steps = number_of_steps + 1;
            break
        end
        diff = window1_sum - window2_sum;
        absDiff = abs(diff);
        if( sumMinStepSize <= absDiff)
            %level_amp = mean(event_data(step_matrix(end,1):(i-(window_size+1))));
            %check if still holds true after an update
            lastDiff = 0;
            ldiff = diff;
            %hill climb to the maximum window difference.
            while absDiff > lastDiff && sign(ldiff) == sign(diff) && i < end_point - window_size - 1;
                lastDiff = absDiff;
                ldiff = diff;
                i = i+1;
                window1_sum = window1_sum - (event_data(i-window_size))+ (event_data(i));
                window2_sum = window2_sum - (event_data(i+1)) + (event_data(i+window_size+1));
                diff = window1_sum - window2_sum;
                absDiff = abs(diff);
            end
            %mark this point as a jump.
            number_of_steps = number_of_steps + 1;
            step_matrix(end+1,1) = i;
            %step_matrix(end,2) = sign(diff);
            newdiff = diff;
            %now burn through some space until we are at the end of
            %a step.
            stopVal = end_point - window_size - 1;
            absStop = sumMinStepSize/2;
            while sign(newdiff) == sign(diff) && abs(newdiff) >= absStop && i < stopVal
                i = i+1;
                window1_sum = window1_sum - (event_data(i-window_size))+ (event_data(i));
                window2_sum = window2_sum - (event_data(i+1)) + (event_data(i+window_size+1));
                newdiff = window1_sum - window2_sum;
            end
            
        end
    end
    if step_matrix(end,1) ~= end_point
        step_matrix(end+1,1) = end_point;
    end
    
    
    
else
    stpSize = 0;
    leftStepLength = 0;
    rightStepLength = 0;
    leftStepAmp = wv_amp;
    rightStepAmp = wv_amp;
    h=1;
end
%%%%%%%%%
% JSJ: This algorithm requires post processing of discovered
% Potential breakpoints. Here is where that occures.
if number_of_steps > 2 %need more sophisticated step combining algorithm JSJ
    if double_check
        i = 1;
        numLeft = size(step_matrix,1)-1;
        while i < numLeft
            i = i + 1;
            h = step_matrix(i,1);

            leftStepAmp = mean(filt(step_matrix(i-1,1):step_matrix(i,1)));
            rightStepAmp = mean(filt(step_matrix(i,1):step_matrix(i+1,1)));

            leftLength = step_matrix(i,1) - step_matrix(i-1,1);
            rightLength = step_matrix(i+1,1) - step_matrix(i,1);
            %F = @(xmean,ymean,xlen,ylen) ((xmean -
            %ymean)^2)/(((1/xlen)+(1/ylen))*var);
            p = erfc(sqrt(F(leftStepAmp,rightStepAmp,leftLength,rightLength)/2));

            if( p >= pval)
                %stepsToDel = [stepsToDel;i];
                numLeft = numLeft - 1;
                number_of_steps = number_of_steps - 1;
                step_matrix(i,:)=[];
                i = i - 1;
                %step_matrix(i,:)=[];
            end

        end
    end
    if label_states
        coord = struct('start',{-1},'end',{-1});
        state = struct('label',{-1},'mean',{-1},'length',{-1},'coords',{[coord]});
        %% Get first state
        stateCounter = 1;
        state1 = struct(state);
        %state1 = struct();
        stepAmp = mean(filt(step_matrix(1,1):step_matrix(2,1)));
        coordtmp = struct(coord);
        coordtmp.start = step_matrix(1,1);
        coordtmp.end = step_matrix(2,1);
        state1.coords = [ coordtmp ];
        state1.mean = stepAmp;
        state1.label = stateCounter;
        step_matrix(1,2) = state1.label;
        state1.length = step_matrix(2,1) - step_matrix(1,1);
        states = [state1];
        
        %% Get rest of states
        i = 1;
        while i < size(step_matrix,1)-1
            i = i+1;
            %get minimum difference between the current state amp, and the
            %seen state amps
            minIndex = 1;
            stateAmp = states(1).mean;
            currAmp = mean(filt(step_matrix(i,1):step_matrix(i+1,1)));
            minDiff = abs(stateAmp - currAmp);
            for j=2:stateCounter
                stateAmp = states(j).mean;
                thisDiff = abs(stateAmp - currAmp);
                if thisDiff < minDiff
                    minDiff = thisDiff;
                    minIndex = j;
                end
            end
            
            %now see if the minimum difference is significant
            p = erfc(sqrt(F(states(minIndex).mean,currAmp,states(minIndex).length, step_matrix(i+1,1) - step_matrix(i,1))/2));
            if p >= pval %same state
                coordtmp = struct(coord);
                %coordtmp = struct();
                coordtmp.start = step_matrix(i,1);
                coordtmp.end = step_matrix(i+1,1);
                states(minIndex).coords(end+1) =  coordtmp;
                tmpSum = states(minIndex).mean * states(minIndex).length;
                tmpSum = tmpSum + sum(filt(step_matrix(i,1):step_matrix(i+1,1)));
                states(minIndex).length = states(minIndex).length + step_matrix(i+1,1) - step_matrix(i,1);
                states(minIndex).mean = tmpSum / states(minIndex).length;
                step_matrix(i,2) = states(minIndex).label;
                
            else
                %make a new state!
                stateCounter = stateCounter + 1;
                nstate = struct(state);
                %nstate = struct();
                nstate.mean = currAmp;
                nstate.length = step_matrix(i+1,1) - step_matrix(i,1);
                nstate.label = stateCounter;
                step_matrix(i,2) = nstate.label;
                nCoords = struct(coord);
                %nCoords = struct();
                nCoords.start = step_matrix(i,1);
                nCoords.end = step_matrix(i+1,1);
                nstate.coords = [nCoords];
                states(end+1) = nstate;
                
                
            end
            
        end
    end
    
    %% rest
    
    
    %uncommont this to print the pa time values
    %             if handles.curindex == 1
    %                 format long;
    %                 disp(filt);
    %             end
    % Open the file and append these results to the end
    fid = fopen(strcat(get(handles.edt_save_filename,'String'),'.txt'),'a');
    if label_states
        for i = 2:size(step_matrix,1)
            fprintf(fid, '%20.15f\t%20.15f\t%20.15f\t%d\n', time(step_matrix(i-1,1)),mean(filt(step_matrix(i-1,1):step_matrix(i,1))),time(step_matrix(i,1))-time(step_matrix(i-1,1)), step_matrix(i-1,2));
        end
        fprintf(fid, '%20.15f\t%d\t%d\t%d\n',time(step_matrix(end,1)),-10000, -10000, -10000);%value to signify end of step matrix.
    else
        for i = 2:size(step_matrix,1)
            fprintf(fid, '%20.15f\t%20.15f\t%20.15f\n', time(step_matrix(i-1,1)),mean(filt(step_matrix(i-1,1):step_matrix(i,1))),time(step_matrix(i,1))-time(step_matrix(i-1,1)));
        end
        fprintf(fid, '%20.15f\t%d\t%d\n',time(step_matrix(end,1)),-10000, -10000);%value to signify end of step matrix.
    end
    fclose(fid);
    
    
    
    for i=2:size(step_matrix,1)
        %Construct our multiple step info for when we have lots of steps
        processedSteps(i-1,:) = [mean(filt(step_matrix(i-1,1):step_matrix(i,1))) time(step_matrix(i,1))-time(step_matrix(i-1,1))];
        %         multi_steps_amp(i-1) = mean(filt(step_matrix(i-1,1):step_matrix(i,1)));
        %         multi_steps_dur(i-1) = time(step_matrix(i,1))-time(step_matrix(i-1,1));
    end
    
else
    processedSteps = [];
end
% End the post processing section
%%%%%%%%%%%%





%     set(handles.sweep_plot,'YLim',yaxis);
%     ydata = [fin_avg fin_avg];
%     xdata = get(handles.sweep_plot,'XLim');
%     line(xdata,ydata,'color','g');

%Now that we have a step, lets see if it fit our criteria

%Now we will collapse any steps that are smaller than our min step size


%if we have found a step or steps, then we will mark it down and calculate
%stats
r = 0;
l = 0;
h=1;
stpSize = 0;

term_reached = -1;
pre_term_amp = 0;
pre_term_dur = 0;
pre_term_iqr = 0;
term_amp = 0;
term_dur = 0;
term_iqr =0;
term_time = 0;
after_kink = 0;
dissociation_voltage = 0;


% disp('Number of Steps found')
% disp(size(step_matrix, 1)-2)
% disp(time(step_matrix(:,1)).*1000)
% disp('Multiple Steps Amp')
% disp(multi_steps_amp)
% disp('Multiple Steps Dur');
% disp(multi_steps_dur);


% This is what goes into processed sweeps (handles.sweep_data)
processedSweep = [term_reached, term_amp, term_dur, term_iqr, ...
    wv_amp, wv_dur, wv_iqr, pre_term_amp, pre_term_dur, ...
    pre_term_iqr, term_time, after_kink, dissociation_voltage];
guidata(hObject,handles);
