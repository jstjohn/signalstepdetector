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


%See if we need to ignore the terminal step.
if get(handles.edt_ignore_term, 'Value') == get(handles.edt_ignore_term,'Max')
    contents = get(handles.term_algo,'String');
    term_algo = contents{get(handles.term_algo,'Value')};
    if strcmp(term_algo,'Chisquared')
        ignore_term = true;
        my_ignore_term = false;
    else
        ignore_term = false;
        my_ignore_term = true;
%         set(handles.edt_term_amp,'Visible','off');        
%      set(handles.edt_term_p,'Visible','off');
        terminal_amp = str2double(get(handles.edt_term_amp,'String'));
        terminal_p = str2double(get(handles.edt_term_p,'String'));
    end
else
    my_ignore_term = false;
    ignore_term = false;
end



window_size= str2num(get(handles.edtWinLen,'String'));
var = str2num(get(handles.edtVar,'String'));
F = @(xmean,ymean,xlen,ylen) ((xmean - ymean)^2)/(((1/xlen)+(1/ylen))*var);
%erfc((x^2)

%now to calculate an appropriate cutoff given our window size, (but only for the "double" cases)
%we know both window lengths, and we need to find out the appropriate
%amplitude difference cutoff value such that p < 0.05 to conlcude that the
%shift in mean amplitudes was not by chance.

pval = str2num(get(handles.edtPVal,'String'));


if isempty(get(handles.txt_min_step_size,'String'))
    %Solving for this in Wolfram-Alpha solve for M, erfc((M^2)/((2*V)/L)) =
    %0.05. I get M = +/- (1.66487 * sqrt(var)/sqrt(L)) Where M is the
    %difference in means, L is the window length, and V is the variance.
    query = strcat('erfc(sqrt(((M^2)/((2*',num2str(var),')/',num2str(window_size),'))/2))=',num2str(pval));
    %minStepSize = 3.64277 * sqrt(var)/sqrt(window_size);
    %inStepSize = solve();
    %disp('minStepSize '+minStepSize)
    result=solve(query);
    %the result will be in the form number, -number, and type 'sym'. First we
    %take the absolute value of the number, and then we cast it to a double so
    %that we can use in as a window mean difference cutoff value.
    minStepSize=abs(double(result(1)));
    set(handles.txt_min_step_size,'String',num2str(minStepSize));
else
    minStepSize = str2double(get(handles.txt_min_step_size,'String'));
    
end

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
    fprintf(fid,'#ignore_terminal_step:');
    if ignore_term || my_ignore_term %ignore the terminal step?
        fprintf(fid, 'Yes\n');
        fprintf(fid, '#terminal_removal_algo:');
        fprintf(fid, strcat(term_algo,'\n'));
        if my_ignore_term
            fprintf(fid,'#terminal_removal_p_value:%e\n',terminal_p);
            fprintf(fid,'#terminal_removal_amplitude:%f\n',terminal_amp);
        end
        
    else
        fprintf(fid,'No\n');
    end
    fprintf(fid, strcat('#file_name:',handles.filename,'\n\n'));
    if label_states
        fprintf(fid, 'Start_Time\tEnd_Time\tStep_Time\tMean_Amplitude\tState_Label_By_Order_Seen_In_Sweep\n');
    else
        fprintf(fid, 'Start_Time\tEnd_Time\tStep_Time\tMean_Amplitude\n');
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


term_step_matrix = [];
if (ignore_term)
    term_step_matrix = termStep(handles,trace_data,time,voltage);
    end_point = term_step_matrix(2,1); %start point of term step
end
if ((ignore_term && (size(term_step_matrix,1) > 2)) || ~ignore_term )
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
        while i < (end_point-(window_size)-1)
            i = i + 1;
            %simulate advancing the window forward one position...
            %window_mean = window_mean + ((event_data(i-window_size-1)-event_data(i))/window_size);
            window1_sum = window1_sum - (event_data(i-window_size))+ (event_data(i));
            window2_sum = window2_sum - (event_data(i+1)) + (event_data(i+window_size+1));
            %Handle the end case!
            %
            if i >= end_point-window_size
                step_matrix(end,2) = end_point;
                break
            end
            diff = window1_sum - window2_sum;
            absDiff = abs(diff);
            if( sumMinStepSize <= absDiff)
                
                %maximize the slope
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
                
                
                %Now move the right window forward until the difference stops
                %growing
                lastDiff = 0;
                backPoint = i;
                while absDiff > lastDiff && i < (end_point - (window_size) - 1)
                    lastDiff = absDiff;
                    i = i + 1;
                    window2_sum = window2_sum - (event_data(i+1)) + (event_data(i+window_size+1));
                    absDiff = abs(window1_sum - window2_sum);
                end %i points to the end point of the transition
                
                
                %Then move the left window backwards until the difference stops
                %growing even further
                lastDiff = 0;
                while absDiff > lastDiff && backPoint > start_point + (window_size) + 1 && backPoint > step_matrix(end,1)+1
                    lastDiff = absDiff;
                    window1_sum = window1_sum - event_data(backPoint) + event_data(backPoint-window_size-1);
                    backPoint = backPoint - 1;
                    absDiff = abs(window1_sum - window2_sum);
                end %backPoint points to the start of the transition
                %(and end of last state)
                
                step_matrix(end,2) = backPoint;
                step_matrix(end+1,1) = i;
                
                
                window1_sum = sum(filt(i-window_size:i));
                diff = window1_sum - window2_sum;
                %now burn through some space until we are at the end of
                %a step.
                stopVal = end_point - window_size - 2;
                absStop = sumMinStepSize;
                newDiff = diff;
                while sign(newDiff) == sign(diff) && abs(newDiff) >= absStop && i < stopVal
                    i = i+1;
                    window1_sum = window1_sum - (event_data(i-window_size))+ (event_data(i));
                    window2_sum = window2_sum - (event_data(i+1)) + (event_data(i+window_size+1));
                    newDiff = window1_sum - window2_sum;
                end
                
                
                
                %             %level_amp = mean(event_data(step_matrix(end,1):(i-(window_size+1))));
                %             %check if still holds true after an update
                %             lastDiff = 0;
                %             ldiff = diff;
                %             %hill climb to the maximum window difference.
                %             while absDiff > lastDiff && sign(ldiff) == sign(diff) && i < end_point - window_size - 1;
                %                 lastDiff = absDiff;
                %                 ldiff = diff;
                %                 i = i+1;
                %                 window1_sum = window1_sum - (event_data(i-window_size))+ (event_data(i));
                %                 window2_sum = window2_sum - (event_data(i+1)) + (event_data(i+window_size+1));
                %                 diff = window1_sum - window2_sum;
                %                 absDiff = abs(diff);
                %             end
                %             %mark this point as a jump.
                %             number_of_steps = number_of_steps + 1;
                %             step_matrix(end+1,1) = i;
                %             %step_matrix(end,2) = sign(diff);
                %             newdiff = diff;
                %             %now burn through some space until we are at the end of
                %             %a step.
                %             stopVal = end_point - window_size - 1;
                %             absStop = sumMinStepSize/2;
                %             while sign(newdiff) == sign(diff) && abs(newdiff) >= absStop && i < stopVal
                %                 i = i+1;
                %                 window1_sum = window1_sum - (event_data(i-window_size))+ (event_data(i));
                %                 window2_sum = window2_sum - (event_data(i+1)) + (event_data(i+window_size+1));
                %                 newdiff = window1_sum - window2_sum;
                %             end
                
            end %end if we have seen a jump
        end
        
        
    else
        stpSize = 0;
        leftStepLength = 0;
        rightStepLength = 0;
        leftStepAmp = wv_amp;
        rightStepAmp = wv_amp;
        h=1;
    end
    
    if step_matrix(end,2) == 0
        step_matrix(end,2) = end_point;
    end
    
    if ignore_term && size(term_step_matrix,1) >= 3
        step_matrix(end+1,1) = term_step_matrix(2,1);
        step_matrix(end,2) = term_step_matrix(3,1);
    end
    %%%%%%%%%
    % JSJ: This algorithm requires post processing of discovered
    % Potential breakpoints. Here is where that occures.
    if size(step_matrix,1) > 0 %need more sophisticated step combining algorithm JSJ
        if double_check && size(step_matrix,1) > 1
            i = 1;
            numLeft = size(step_matrix,1);
            while i < numLeft
                i = i + 1;
                h = step_matrix(i,1);
                
                leftStepAmp = mean(filt(step_matrix(i-1,1):step_matrix(i-1,2)));
                rightStepAmp = mean(filt(step_matrix(i,1):step_matrix(i,2)));
                
                leftLength = step_matrix(i-1,2) - step_matrix(i-1,1);
                rightLength = step_matrix(i,2) - step_matrix(i,1);
                %F = @(xmean,ymean,xlen,ylen) ((xmean -
                %ymean)^2)/(((1/xlen)+(1/ylen))*var);
                p = erfc(sqrt(F(leftStepAmp,rightStepAmp,leftLength,rightLength)/2));
                
                if( p >= pval)
                    %stepsToDel = [stepsToDel;i];
                    numLeft = numLeft - 1;
                    step_matrix(i-1,2) = step_matrix(i,2);
                    step_matrix(i,:)=[];
                    i = i - 1; %don't move forward one i, just stay here...
                    %step_matrix(i,:)=[];
                end
                
            end
        end
        
        %% ignore terminal steps by this simple algorithm if selected
        if my_ignore_term
            i = 0;
            numLeft = size(step_matrix,1);
            while i < numLeft
                i = i + 1;
                
                leftStepAmp = mean(filt(step_matrix(i,1):step_matrix(i,2)));
                rightStepAmp = terminal_amp;
                
                leftLength = step_matrix(i,2) - step_matrix(i,1);
                rightLength = step_matrix(end,2) - step_matrix(i,2); %assume the rest is baseline
                %F = @(xmean,ymean,xlen,ylen) ((xmean -
                %ymean)^2)/(((1/xlen)+(1/ylen))*var);
                p = erfc(sqrt(F(leftStepAmp,rightStepAmp,leftLength,rightLength)/2));
                
                if( p >= terminal_p || leftStepAmp < terminal_amp)
                    %stepsToDel = [stepsToDel;i];
                    
                    step_matrix(i:end,:)=[];
                    break; %everything else is terminal step.
                    %step_matrix(i,:)=[];
                end
                
            end
        end
        
        
        %    number_of_steps = size(step_matrix,1);
        if label_states && size(step_matrix,1) > 0
            coord = struct('start',{-1},'end',{-1});
            state = struct('label',{-1},'mean',{-1},'length',{-1},'coords',{[coord]});
            %% Get first state
            stateCounter = 1;
            state1 = struct(state);
            %state1 = struct();
            stepAmp = mean(filt(step_matrix(1,1):step_matrix(1,2)));
            coordtmp = struct(coord);
            coordtmp.start = step_matrix(1,1);
            coordtmp.end = step_matrix(1,2);
            state1.coords = [ coordtmp ];
            state1.mean = stepAmp;
            state1.label = stateCounter;
            step_matrix(1,3) = state1.label;
            state1.length = step_matrix(1,2) - step_matrix(1,1);
            states = [state1];
            
            %% Get rest of states
            i = 1;
            while i < size(step_matrix,1)
                i = i+1;
                %get minimum difference between the current state amp, and the
                %seen state amps
                minIndex = 1;
                stateAmp = states(1).mean;
                currAmp = mean(filt(step_matrix(i,1):step_matrix(i,2)));
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
                p = erfc(sqrt(F(states(minIndex).mean,currAmp,states(minIndex).length, step_matrix(i,2) - step_matrix(i,1))/2));
                if p >= pval %same state
                    coordtmp = struct(coord);
                    %coordtmp = struct();
                    coordtmp.start = step_matrix(i,1);
                    coordtmp.end = step_matrix(i,2);
                    states(minIndex).coords(end+1) =  coordtmp;
                    tmpSum = states(minIndex).mean * states(minIndex).length;
                    tmpSum = tmpSum + sum(filt(step_matrix(i,1):step_matrix(i,2)));
                    states(minIndex).length = states(minIndex).length + step_matrix(i,2) - step_matrix(i,1);
                    states(minIndex).mean = tmpSum / states(minIndex).length;
                    step_matrix(i,3) = states(minIndex).label;
                    
                else
                    %make a new state!
                    stateCounter = stateCounter + 1;
                    nstate = struct(state);
                    %nstate = struct();
                    nstate.mean = currAmp;
                    nstate.length = step_matrix(i,2) - step_matrix(i,1);
                    nstate.label = stateCounter;
                    step_matrix(i,3) = nstate.label;
                    nCoords = struct(coord);
                    %nCoords = struct();
                    nCoords.start = step_matrix(i,1);
                    nCoords.end = step_matrix(i,2);
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
        if size(step_matrix,1) > 0
            fid = fopen(strcat(get(handles.edt_save_filename,'String'),'.txt'),'a');
            if ignore_term
                if label_states
                    for i = 1:size(step_matrix,1)-1 %ignore the terminal step
                        fprintf(fid, '%20.15f\t%20.15f\t%20.15f\t%20.15f\t%d\n', time(step_matrix(i,1)),time(step_matrix(i,2)),time(step_matrix(i,2))-time(step_matrix(i,1)),mean(filt(step_matrix(i,1):step_matrix(i,2))), step_matrix(i,3));
                    end
                    fprintf(fid, '%20.15f\t%20.15f\t%d\t%d\t%d\n',time(end_point),time(end_point),0, -10000, -10000);%value to signify end of step matrix.
                else
                    for i = 1:size(step_matrix,1)-1 %ignore the terminal step
                        fprintf(fid, '%20.15f\t%20.15f\t%20.15f\t%20.15f\n', time(step_matrix(i,1)),time(step_matrix(i,2)),time(step_matrix(i,2))-time(step_matrix(i,1)),mean(filt(step_matrix(i,1):step_matrix(i,2))));
                    end
                    fprintf(fid, '%20.15f\t%20.15f\t%d\t%d\n',time(end_point),time(end_point),0, -10000);%value to signify end of step matrix.
                end
            else
            end
            if label_states
                for i = 1:size(step_matrix,1)
                    fprintf(fid, '%20.15f\t%20.15f\t%20.15f\t%20.15f\t%d\n', time(step_matrix(i,1)),time(step_matrix(i,2)),time(step_matrix(i,2))-time(step_matrix(i,1)),mean(filt(step_matrix(i,1):step_matrix(i,2))), step_matrix(i,3));
                end
                fprintf(fid, '%20.15f\t%20.15f\t%d\t%d\t%d\n',time(end_point),time(end_point),0, -10000, -10000);%value to signify end of step matrix.
            else
                for i = 1:size(step_matrix,1)
                    fprintf(fid, '%20.15f\t%20.15f\t%20.15f\t%20.15f\n', time(step_matrix(i,1)),time(step_matrix(i,2)),time(step_matrix(i,2))-time(step_matrix(i,1)),mean(filt(step_matrix(i,1):step_matrix(i,2))));
                end
                fprintf(fid, '%20.15f\t%20.15f\t%d\t%d\n',time(end_point),time(end_point),0, -10000);%value to signify end of step matrix.
            end
            fclose(fid);
        end
        
        
        step = 1;
        for i=1:size(step_matrix,1)
            %Construct our multiple step info for when we have lots of steps
            processedSteps(step,:) = [mean(filt(step_matrix(i,1):step_matrix(i,2))) time(step_matrix(i,2))-time(step_matrix(i,1))];
            
            step = step + 1;
            
            %Now lets add in a pseudoState for the transition
            if i < size(step_matrix,1)
                processedSteps(step,:) = [-100  time(step_matrix(i+1,1))-time(step_matrix(i,2))];
                step = step+1;
            end
            %         multi_steps_amp(i-1) = mean(filt(step_matrix(i-1,1):step_matrix(i,1)));
            %         multi_steps_dur(i-1) = time(step_matrix(i,1))-time(step_matrix(i-1,1));
        end
        
    else
        processedSteps = [];
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
