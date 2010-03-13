%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----Checks for a terminal step
% handles        handles of GUI
% index          index in trace_data of sweep we want to look at
% plot_yn        'yes' or 'no' for whether we want to display the check or not
% term_reached   -1 never cross terminal threshold, 0 cross term threshold but
%                ruled out as a real term step, 1 term step
% function [term_reached, term_amp, term_dur, term_iqr, ... 
%     wv_amp, wv_dur, wv_iqr, pre_term_amp, pre_term_dur, ...
%     pre_term_iqr] = analyize_term(hObject,handles,index,plot_yn, force_term)
function processedSweep = analyze_term(hObject,handles,index,plot_yn, force_term)
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
term_time = 0;
epsilon = 0.4; %used to be 0.8 at 18 june 2008 16:23
amp_thresh = 0;
% epsilon = 0.2;

deriv_window = 40; %make it even plz
deriv_alpha = .9995; %filtering on derivative data

if ~strcmpi(plot_yn,'yes') && ~strcmpi(plot_yn,'no')
    disp('must enter yes or no for plot option');
    return;
end

[trace_data, time] = get_event(handles, hObject, index);

handles = guidata(hObject);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPERIMENTAL CHISQUARED STEP FINDER
% 
% data = double(trace_data(:,1));
% 
% % X0 = [double(mean(trace_data(:,1))); double(mean(trace_data(:,1))); double(time(round(length(trace_data(:,1))/2)))];
% % X0 = [double(mean(trace_data(:,1))); double(mean(trace_data(:,1))); double(round(length(trace_data(:,1))/2))];
% % newoptions=optimset('display','iter');
% % [X,FVAL,EXITFLAG,OUTPUT] = fminsearch(@(x) chisquaredfun(x,data, double(time)),X0, newoptions);
% 
% i = round(length(data)/2);
% n = length(data);
% 
% A = [ones(i,1) zeros(i,1); zeros((n-i),1) ones((n-i),1)];
% 
% V = pinv(A)*data

% if strcmpi(plot_yn,'yes')
%     hold on;
%     plot([0 X(3)], [X(1) X(1)],'-k')
%     plot([X(3) time(end)], [X(2) X(2)], '-k');
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stop_point =  length(trace_data);% - round((0.2/1000)/handles.time_tick);

if( handles.event_trans == 0 )
    start_point = 1;
%     start_point = round((0.6/1000)/handles.time_tick);
    end_point = length(trace_data);
    baseline_offset = 1; %in num of data points
else
%    start_point = ceil((handles.event_trans/1000)/handles.time_tick);
    start_point = handles.detected_events(index,9); % + handles.transAddSamples;
    if(start_point == 0)
	    start_point = 1;
    end

    end_point = length(trace_data); %-ceil((handles.event_trans/1000)/handles.time_tick)+1;
    %baseline_offset = 1+start_point; %in num of data points
    baseline_offset = handles.transAddSamples+start_point; %in num of data points
end
% 
% disp(start_point)
% disp(end_point)
% disp(handles.time_tick)
% disp(floor((handles.event_trans/1000)/handles.time_tick))
if( force_term == 1 )
    term_reached = 1;
    wv_amp = sum(trace_data(start_point:end_point))/length(trace_data(start_point:end_point));
    wv_dur = (end_point)*handles.time_tick;
    wv_iqr = iqr(trace_data);
    midpoint = round(handles.termStep/handles.time_tick);
    term_time = time(midpoint);
    term_amp = sum(trace_data(midpoint:end_point))/length(midpoint:end_point);
    term_dur = (end_point-midpoint)*handles.time_tick;
    term_iqr = iqr(trace_data(midpoint:end_point));
    pre_term_dur = (midpoint)*handles.time_tick;
    pre_term_amp = sum(trace_data(start_point:midpoint))/length(start_point:midpoint);
    pre_term_iqr = iqr(trace_data(start_point:midpoint));
    if strcmpi(plot_yn,'yes')
        set(handles.edt_thresh, 'String', [num2str(wv_amp) ' pA']);
        set(handles.edt_termstep, 'String', 'Yes');
        set(handles.edt_termamp,'String', [num2str(term_amp) ' pA']);
        set(handles.edt_termdur,'String', [num2str(term_dur*1000) ' ms']);
        set(handles.edt_pretermamp,'String', [num2str(pre_term_amp) ' pA']);
        set(handles.edt_pretermdur,'String', [num2str(pre_term_dur*1000) ' ms']);
    end
    processedSweep = [term_reached, term_amp, term_dur, term_iqr, ... 
    wv_amp, wv_dur, wv_iqr, pre_term_amp, pre_term_dur, ...
    pre_term_iqr, term_time];
    return
end

alpha = handles.alpha;

if ( (time(end) <= 1e-3) || (start_point >= end_point) || (time(end_point) - time(start_point) <= 1e-3))
    wv_amp = sum(trace_data(start_point:end_point))/length(trace_data(start_point:end_point));
    wv_dur = (end_point)*handles.time_tick;
    wv_iqr = iqr(trace_data);
    term_reached = -1;
    pre_term_amp = 0;
    pre_term_dur = 0;
    pre_term_iqr = 0;
    term_amp = 0;
    term_dur = 0;
    term_iqr =0;
    term_time = 0;
    processedSweep = [term_reached, term_amp, term_dur, term_iqr, ... 
    wv_amp, wv_dur, wv_iqr, pre_term_amp, pre_term_dur, ...
    pre_term_iqr, term_time];
    return
end
%exponential filtering

% ind = 1;
% filt = zeros(length(1:end_point-2),1);
% filt(ind) = handles.trace_data(1,index+1);
% 
% for i=2:end_point-2
%     filt(ind+1) = handles.trace_data(i,index+1)*alpha + (1-alpha)*filt(ind);
%     ind = ind+1;
% end

%We'll cut off the last two points to help w/ our detection
filt = exp_filt(trace_data,1,length(trace_data),alpha);
%Get the basline from whatever we want the baseline window to be

temp_data = guidata(handles.settings);
t_points = length(trace_data);

switch get(temp_data.popmnu_baseline_type,'Value')
    case 1
        handles.baseline_window = str2num(get(temp_data.edt_baseline_window,'String'));
    case 2
        handles.baseline_window = ((str2num(get(temp_data.edt_baseline_window,'String'))*10^-3)/handles.time_tick)/t_points;
end

if( (ceil(t_points*handles.baseline_window)+baseline_offset) > length(trace_data) )
    avg = sum(trace_data(baseline_offset:end));
    fin_avg = avg/size(trace_data(baseline_offset:end),1);
else
    avg = sum(trace_data(baseline_offset:(ceil(t_points*handles.baseline_window)+baseline_offset)));
    fin_avg = avg/size(trace_data(baseline_offset:(ceil(t_points*handles.baseline_window)+baseline_offset)),1);
end

% disp(fin_avg)
% if fin_avg > amp_thresh
%     disp('yes')
% end
wv_amp = sum(trace_data(start_point:end_point))/length(trace_data(start_point:end_point));
wv_dur = end_point*handles.time_tick;
wv_iqr = iqr(trace_data(start_point:end_point));
deriv = zeros(size(filt));

if( handles.rem_trans )
    for i=deriv_window+1:length(filt)
        deriv(i-deriv_window) = (filt(i) - filt(i-deriv_window))/deriv_window;
    end

    filt_deriv = exp_filt(deriv,1,length(deriv),deriv_alpha);

    %we're going to build our breakpoint vector

    if( time(end-2) > 0.3 )
        b = [0:0.005:0.02 0.03:0.01:0.05 0.10:0.05:0.3];
        b = [b 0.4:0.1:(floor((time(end-2)-0.3)/0.1)*0.1+0.3)];
        if(time(end-2)-b(end) > 0.3/2)
            b(end+1) = time(end-2);
        end
    elseif( time(end-2) > 0.05 )
        b = [0:0.005:0.02 0.03:0.01:0.05];
        b = [b 0.1:0.05:(floor((time(end-2)-0.05)/0.05)*0.05+0.05)];
        if(time(end-2)-b(end) > 0.05/2)
            b(end+1) = time(end-2);
        end
    elseif( time(end-2) > 0.02 )
        b = 0:0.005:0.02;
        b = [b 0.03:0.01:(floor((time(end-2)-0.02)/0.01)*0.01+0.02)];
        if(time(end-2)-b(end) > 0.02/2)
            b(end+1) = time(end-2);
        end
    else
        b = 0:0.005:floor(time(end-2)/0.005)*0.005;
        b(end+1) = time(end-2);
    end

%     disp(b);
    %         tic
    fitfnc = spline(b,filt(start_point:end_point)'/spline(b,eye(length(b)),time(start_point:end_point)'));
    fitvals = fnval(time(start_point:end_point),fitfnc);
    derfnc = fnder(fitfnc);
    dervals = fnval(time(start_point:end_point),derfnc);

    if( time(end_point) > 0.02)
        %             disp(min(fitvals(floor(0.02/handles.time_tick):end)))
        [close_zero, zero_ind] = min(abs(dervals(floor(0.02/handles.time_tick):end)));
        flatend_filt = filt(start_point:end_point)-(fitvals-fitvals(zero_ind+floor(0.02/handles.time_tick)-1))';
        %             flatend_filt = filt-(fitvals-min(fitvals(floor(0.02/handles.time_tick):end)))';
        inside_iqr = sort(filt(start_point:end_point)-(fitvals-min(fitvals(floor(0.02/handles.time_tick):end)))');
    else
        flatend_filt = filt(start_point:end_point)-(fitvals-min(fitvals))';
        inside_iqr = sort(filt(start_point:end_point)-(fitvals-min(fitvals))');
    end

    inside_iqr = inside_iqr(floor(length(inside_iqr)*0.25):floor(length(inside_iqr)*.75));

    %         fin_avg = median(inside_iqr);
    avg = sum(flatend_filt(baseline_offset:(ceil(t_points*handles.baseline_window)+baseline_offset)));
    fin_avg = avg/size(flatend_filt(baseline_offset:(ceil(t_points*handles.baseline_window)+baseline_offset)),1);

%             figure(3);plot(time(1:end-2),filt,'c',time(1:end-2),dervals/100,'r',time(1:end-2),fitvals,'k'); grid on;%axis([0 time(end-2) -15 50]);grid on;
%     %         toc
%             figure(2);
%             plot(time(1:end-2), flatend_filt, 'r');
%             grid on
    axes(handles.sweep_plot);
else
    fitvals = 0;
end
    
if strcmpi(plot_yn,'yes')
    axes(handles.sweep_plot);
    hold on;
    yaxis = get(handles.sweep_plot,'YLim');
    plot(time(start_point:end_point), filt(start_point:end_point),'-r')
    if(fitvals ~= 0)
        plot(time(start_point:end_point),fitvals,'c');
    end
    set(handles.sweep_plot,'YLim',yaxis);
    ydata = [fin_avg fin_avg];   
    xdata = get(handles.sweep_plot,'XLim');
    line(xdata,ydata,'color','g');
    plot(time(baseline_offset),trace_data(baseline_offset),'co');
    plot(time(start_point),trace_data(start_point),'gs');
    plot(time(stop_point), trace_data(stop_point),'gs');

    if( (ceil(t_points*handles.baseline_window)+baseline_offset) > length(trace_data) )
        plot(time(end),trace_data(end),'co');
    else
        plot(time(ceil(t_points*handles.baseline_window)+baseline_offset), trace_data(ceil(t_points*handles.baseline_window)+baseline_offset),'co');
    end
    hold off;
    set(handles.edt_thresh,'String',[num2str(wv_amp) ' pA']);
end
num_7ms_points = ceil(7e-3/handles.time_tick);
num_4ms_points = ceil(4e-3/handles.time_tick);
num_2ms_points = ceil(2e-3/handles.time_tick);
num_p2ms_points = ceil(2e-4/handles.time_tick);

trigger = fin_avg + handles.thresh_amp;     %this will be our terminal trigger
                                            %we'll check it at what we set
if strcmpi(plot_yn,'yes')
    axes(handles.sweep_plot);
    hold on;
    ydata = [trigger trigger];
    xdata = get(handles.sweep_plot,'XLim');
%     xdata(1) = trace_data(1);
    line(xdata,ydata,'color','g','linestyle','--');
    hold off;
end

%First pass at determining terminal step, historesis thresholding, we'll
%look back 2 ms from the end, if we don't find a terminal step threshold
%within .2 ms from the end of the trace we'll say couldn't find one
%---------------------
%Second Iteration - We look at the whole waveform w/ the threshold we've
%setup. But now we have a 'flip-flop' bit that decides whether we have a
%thershold or not, we don't break out of our search if we think its a false
%detect. If we have a positive response we'll work backwards from the
%previous 'detect' and average the waveform to the end, if the average is within a
%certain amount of the thershold we'll assume that where the terminal step
%really started and adjust our number accordingly

term_reached = -1;
term_index = [];
aboveThresh = find(filt(start_point:end_point) > fin_avg);
if ~isempty(aboveThresh)
%     firstPoint = find((filt((aboveThresh(1)+start_point-1):end_point) < fin_avg), 1);
      start_point = aboveThresh(1)+start_point-1;
%     if ~isempty(firstPoint)
%         start_point = firstPoint + aboveThresh(1)+start_point-1;
%     else
% 	start_point = end_point;
%     end
end
belowThresh = find(filt(start_point:stop_point) < trigger);
aboveThresh = find(filt(start_point:stop_point) > fin_avg);

% disp(belowThresh(1:5))
%since starting at start_point, re-normalized indices
if fin_avg > amp_thresh
    belowThresh = belowThresh + (start_point-1);
    aboveThresh = aboveThresh + (start_point-1);
else
    belowThresh = [];
    aboveThresh = [];
end

if ~isempty(belowThresh) && ~isempty(aboveThresh)
    aboveThreshInd = 1;
    belowThreshInd = 1;
    threshCrosses = zeros(5000,1);
    threshCrosses(1) = belowThresh(1);
    threshInd = 2;

    %Generate a vector where the first point is the first crossing of the
    %term step threshold, the second point is the next crossing of the
    %baseline, and so forth. Alternating crosses of each threshold starting
    %w/ the terminal step threshold. So, if the vector is odd then it has
    %crossed the terminal step threshold and not come back.
    while(1)

        aboveThreshIndnew = find(aboveThresh(aboveThreshInd:end) > belowThresh(belowThreshInd),1);
        aboveThreshInd = aboveThreshInd + (aboveThreshIndnew-1);
        if(isempty(aboveThreshIndnew))
            break;
        end
        threshCrosses(threshInd) = aboveThresh(aboveThreshInd);

        belowThreshIndnew = find(belowThresh(belowThreshInd:end) > aboveThresh(aboveThreshInd),1);
        belowThreshInd = belowThreshInd + (belowThreshIndnew-1);
        if(isempty(belowThreshIndnew))
            break;
        end
        threshCrosses(threshInd+1) = belowThresh(belowThreshInd);

        threshInd = threshInd+2;
    end

    threshInd = find(threshCrosses,1,'last');
    threshCrosses(threshInd+1:end) = [];

    %if the first point is 'crossed' the term thresh, then ignore it and
    %the subsequent baseline crossage, it assumes that there must be a
    %baseline cross point since to get that baseline past the first point
    %there must be some point that crosses the baseline

%     disp(threshCrosses)
    
    if (mod(length(threshCrosses),2) == 0)
        term_reached = 0;        
    else
        term_reached = 1;
    end

    threshCross_terms = 1:2:length(threshCrosses);
    threshCross_baseline = 2:2:length(threshCrosses);
    term_index = threshCrosses(threshCross_terms);

    if strcmpi(plot_yn,'yes')
        if(term_reached == 0)
            set(handles.edt_termstep, 'String', 'No');
        else
            set(handles.edt_termstep, 'String', 'Yes');
        end
        hold on;
        plot(time(threshCrosses(threshCross_terms)),filt(threshCrosses(threshCross_terms)),'ks');
        plot(time(threshCrosses(threshCross_baseline)),filt(threshCrosses(threshCross_baseline)), 'ko');
        hold off;
    end

else

    %if isempty(aboveThresh) && ~isempty(belowThresh)
    %    term_reached = 1;
    %    term_index = belowThresh(1);
    %end
    
    if term_reached ~= 1
        if strcmpi(plot_yn,'yes')
            hold on;
            plot(time(end_point-num_p2ms_points),filt(end_point-num_p2ms_points),'go');
            hold off;
            set(handles.edt_termstep, 'String', 'No');
        end
        %            disp('Ran outta time');
        %            disp(filt(end_point-num_p2ms_points-2));
    else
        if strcmpi(plot_yn,'yes')
            set(handles.edt_termstep, 'String', 'Yes');
        end
    end
end

% for i=end_point-num_2ms_points:end_point

% for i=start_point:end_point-num_p2ms_points
%    if term_reached ~= 1 && i >= end_point-num_p2ms_points
%        if strcmpi(plot_yn,'yes')
%            hold on;
%            plot(time(end_point-num_p2ms_points),filt(end_point-num_p2ms_points),'go');
%            hold off;
%            set(handles.edt_termstep, 'String', 'No');
%        end
% %        disp('Ran outta time');
% %        disp(filt(end_point-num_p2ms_points-2));
%        break;
%    end
%    if filt(i) < trigger && term_reached <= 0
% %         disp('Found Term Step');
%         if strcmpi(plot_yn,'yes')
%            hold on;
%            plot(time(i),filt(i),'ks')
%            hold off;
%            set(handles.edt_termstep, 'String', 'Yes');
%         end
%        term_index = [term_index; i];
%        term_reached = 1;
%    end
%    if term_reached == 1 && filt(i) > fin_avg
% %         disp('Found Term Step but appears to be false');
%         if strcmpi(plot_yn,'yes')
%            hold on;
%            plot(time(i),filt(i),'ko')
%            hold off;
%            set(handles.edt_termstep, 'String', 'False Detect');
%         end
%        %term_index = i;
%        term_reached = 0;
%    end
% 
% end



%Now we'll average what appeared to be the terminal step to see its
%amplitude
if term_reached ~= -1
    term_amp = sum(trace_data(term_index(end):end_point-1))/length(term_index(end):end_point-1);
    term_dur = (end_point-term_index(end))*handles.time_tick;
    term_iqr = iqr(trace_data(term_index(end):end_point-1));
    pre_term_dur = (term_index(end))*handles.time_tick;
    pre_term_amp = sum(trace_data(start_point:term_index(end)))/length(start_point:term_index(end));
    pre_term_iqr = iqr(trace_data(start_point:term_index(end)));

%     disp(['Terminal Step Amp: ' num2str(term_amp) '   Terminal Step Duration: ' num2str(term_dur)]);
else
    term_amp = 0;
    term_dur = 0;
    term_iqr = 0;
    pre_term_dur = 0;
    pre_term_amp = 0;
    pre_term_iqr = 0;
end

%For the false detects, we'll see the average of the amplitude of the
%wrongly detected terminal step, if the average of it is less then our
%threshold amplitude then most likely it just had a spike that went back up
%to our baseline.
%We'll check each of the previous ones also since we might have a correct
%terminal step before a false detect

% disp(term_amp)
% disp(fin_avg +(handles.thresh_amp*(1/2)))
if term_reached ~= -1
    %Can't start with a terminal step for the event
    startind = 1;

    for thresh_ind = startind:length(term_index)
        term_amp = sum(trace_data(term_index(thresh_ind):end_point-1))/length(term_index(thresh_ind):end_point-1);
        
        if term_amp <= fin_avg +(handles.thresh_amp*(1-epsilon))
            term_reached = 1;
            set(handles.edt_termstep, 'String', 'Yes');
            term_index = term_index(thresh_ind);
            break;
        end
    end
end

if strcmp(get(handles.edt_termstep,'String'),'False Detect')
    term_reached = -2;
end

%We'll now adjust where we start counting the terminal step, since right
%now it waits until the threshold is crossed and although this works it
%give more 'time' to the terminal step, we'll also recalculate the amp and
%duration. We'll take the midway point between when the filtered signal
%crossed the level 0 thershold and the trigger threshold

if term_reached == 1
    for i=term_index:-1:1
        if filt(i) >= (fin_avg+(handles.thresh_amp.*0.5))
            %found level '.5' crossing
                midpoint = i;
%                 midpoint = ceil((term_index(end)-i)/2)+i;
                if(start_point >= midpoint)
                    pre_term_amp = 0;
                else
                    term_amp = sum(trace_data(midpoint:end_point-1))/length(midpoint:end_point-1);
                    term_dur = (end_point-midpoint)*handles.time_tick;
                    term_iqr = iqr(trace_data(midpoint:end_point-1));
                    pre_term_dur = (midpoint)*handles.time_tick;
                    pre_term_amp = sum(trace_data(start_point:midpoint))/length(start_point:midpoint);
                    pre_term_iqr = iqr(trace_data(start_point:midpoint));
                    term_time = time(midpoint);
                end
                if strcmpi(plot_yn,'yes')
                   hold on;
                   ydata = get(handles.sweep_plot,'YLim');
                   xdata = [time(midpoint) time(midpoint)];
                   handles.termStepLine = line(xdata, ydata,'color','r','linestyle','-');
                   handles.termStep = time(midpoint);
%                    plot(trace_data(midpoint,1),filt(midpoint),'ro');
                   hold off;
                end
                break;
        end
    end
end



%we'll make the amplitudes in the form of negative pA from the baseline

% wv_amp = wv_amp - handles.level0;
% if term_amp ~= 0 
%     term_amp = term_amp - handles.level0;
%     pre_term_amp = pre_term_amp - handles.level0;
% end

if strcmpi(plot_yn,'yes')
    set(handles.edt_termamp,'String', [num2str(term_amp) ' pA']);
    set(handles.edt_termdur,'String', [num2str(term_dur*1000) ' ms']);
    set(handles.edt_pretermamp,'String', [num2str(pre_term_amp) ' pA']);
    set(handles.edt_pretermdur,'String', [num2str(pre_term_dur*1000) ' ms']);
end

hold off

processedSweep = [term_reached, term_amp, term_dur, term_iqr, ... 
wv_amp, wv_dur, wv_iqr, pre_term_amp, pre_term_dur, ...
pre_term_iqr, term_time];

guidata(hObject,handles);
