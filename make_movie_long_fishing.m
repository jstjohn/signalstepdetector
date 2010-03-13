clear all;

generate_movie = 1;

filename = 'C:\Users\BLGyarfas\Documents\MATLAB\datafiles\09519011_Filt.abf';
% filename = 'C:\Documents and Settings\Administrator\My Documents\MATLAB\datafiles\09519011_Filt.abf';
% filename = 'Z:\pClamp\Station BN Data\09519011_Filt.abf';
video_output_name = 'DNA_Synthesis4.avi';

% sweeps_to_concat = [35:39];
sweeps_to_concat = [41:45];
% sweeps_to_concat = [20:23];
animation_data_raw = [];

%If the start time for a sweep is the beginning, just put 0
% start_time_for_each_sweep = [1.5 0 0 0 0];
start_time_for_each_sweep = [4.4 0 0 0 0]; 
% start_time_for_each_sweep = [0.2 0 0 0]; 

%if end, just put inf
% end_time_for_each_sweep = [inf inf inf inf 2.756];
end_time_for_each_sweep = [inf inf inf inf 3.3];
% end_time_for_each_sweep = [inf inf inf 3.2];

if (length(start_time_for_each_sweep) ~= length(sweeps_to_concat) || ...
        length(end_time_for_each_sweep) ~= length(sweeps_to_concat) )
    disp('start times or end times do match the number of sweeps to concat')
    return
end

for i=1:length(sweeps_to_concat)
    [data, sampleIntOrg] = import_abf(filename, sweeps_to_concat(i));

    if start_time_for_each_sweep(i) == 0
        start_time = 1;
    else
        start_time = round(start_time_for_each_sweep(i)/(sampleIntOrg*1e-6));
    end

    if end_time_for_each_sweep(i) == inf
        end_time = length(data);
    else
        end_time = round(end_time_for_each_sweep(i)/(sampleIntOrg*1e-6));
    end

    animation_data_raw = [animation_data_raw; data(start_time:end_time,:)];

end

down_sample = 1;
down_sample_factor = 5;

if down_sample == 1
    animation_data = animation_data_raw(1:down_sample_factor:end,:);
    time = 0:(sampleIntOrg*down_sample_factor)*1e-6:((sampleIntOrg*down_sample_factor*1e-6*(length(animation_data)-1)));
    sampleInt = sampleIntOrg*down_sample_factor;
else
    animation_data = animation_data_raw;
    time = 0:sampleIntOrg*1e-6:(sampleIntOrg*1e-6*(length(animation_data)-1));
    sampleInt = sampleIntOrg;
end

h1 = figure('Position',[10 50 720 480]);
set(h1, 'DoubleBuffer', 'on', 'Color', [1 1 1]);

if generate_movie == 1
    mov = avifile(video_output_name,'Compression', 'None', 'fps', 30, 'quality', 100);
end

%trying for 30 fps so plot every 33.3 ms

frame_per_sec = 30;
points_to_plot = round((1/15)/(sampleInt*1e-6));

seconds_to_plot = time(end);

frame_counter = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FIRST PART - Real Time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
acurrent = subplot(2,1,1);
axis([0 22 -20 65]);
ylabel('pA')
set(acurrent, 'XTickLabel', '', 'NextPlot', 'replacechildren');
set(acurrent, 'XLimMode', 'manual', 'YLimMode', 'manual');
tempInsetCurrent = get(acurrent, 'TightInset');

avoltage = subplot(2,1,2);
axis([0 22 -50 200]);
xlabel('time (s)')
ylabel('mV');
tempInsetVoltage = get(avoltage, 'TightInset');
set(avoltage, 'XLimMode', 'manual', 'YLimMode', 'manual', 'NextPlot', 'replacechildren');

shared_width = 1-(max([tempInsetCurrent(1) tempInsetVoltage(1)]) ...
    + max([tempInsetCurrent(3) tempInsetVoltage(3) ]));

shared_left = max([tempInsetCurrent(1) tempInsetVoltage(1)]);

set(acurrent,'Position', [shared_left 0.3+tempInsetCurrent(2) shared_width 0.7-(tempInsetCurrent(4)+tempInsetCurrent(2))])
set(avoltage,'Position', [shared_left tempInsetVoltage(2) shared_width 0.3-(tempInsetVoltage(4)+tempInsetVoltage(2))])

axes(acurrent)
hplotc = plot(acurrent, time(1:(points_to_plot+1)), animation_data(1:(points_to_plot+1), 1), '-k');
%         text('Units', 'normalized', 'Position', [0.9 0.9 0], 'String', '1x','FontSize',24, 'Color', [1 0 0])
htext = text('Units', 'pixels', 'Position', [600 300 0], 'String', '0.5x','FontSize',24, 'Color', [1 0 0]);

axes(avoltage)
hplotv = plot(avoltage, time(1:(points_to_plot+1)), animation_data(1:(points_to_plot+1), 2), '-k');

for i=(1+points_to_plot):points_to_plot:(round(seconds_to_plot/(sampleInt*1e-6)-points_to_plot))

    delete(hplotc)
    delete(htext)
    delete(hplotv)

    axes(acurrent)
    hplotc = plot(acurrent, time(1:i), animation_data(1:i, 1), '-k');
    %         text('Units', 'normalized', 'Position', [0.9 0.9 0], 'String', '1x','FontSize',24, 'Color', [1 0 0])
    htext = text('Units', 'pixels', 'Position', [600 300 0], 'String', '0.5x','FontSize',24, 'Color', [1 0 0]);

    axes(avoltage)
    hplotv = plot(avoltage, time(1:i), animation_data(1:i, 2), '-k');

    if generate_movie == 1
        F = getframe(gcf);
        mov = addframe(mov, F);
    else
        pause(0.001)
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FADE OUT OLD TRACE AND CHANGE TEXT VALUE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotchange = zeros(150,3);
plotchange(:,1) = cos(linspace(pi,2*pi,150))./2+.5;
plotchange(:,2) = cos(linspace(pi,2*pi,150))./2+.5;
plotchange(:,3) = cos(linspace(pi,2*pi,150))./2+.5;

textchange = zeros(152,3);
textchange(:,1) = ones(152,1);
textchange(:,2) = cos(linspace(pi,3*pi,152))./2+.5;
textchange(:,3) = cos(linspace(pi,3*pi,152))./2+.5;


for i=1:length(plotchange)
    set(hplotc, 'Color', plotchange(i,:))
    set(hplotv, 'Color', plotchange(i,:))

    set(htext, 'Color', textchange(i,:));
    if i > length(plotchange)/2
        set(htext, 'String', 'x1')
    end
    if generate_movie == 1
        F = getframe(gcf);
        mov = addframe(mov, F);
    else
        pause(0.01)
    end
end

delete(hplotc)
delete(hplotv)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%START UP TRACE AND ZOOM IN AS ITS PLAYIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

down_sample = 0;
down_sample_factor = 5;

if down_sample == 1
    animation_data = animation_data_raw(1:down_sample_factor:end,:);
    time = 0:(sampleIntOrg*down_sample_factor)*1e-6:((sampleIntOrg*down_sample_factor*1e-6*(length(animation_data_raw)-1)));
    sampleInt = sampleIntOrg*down_sample_factor;
else
    animation_data = animation_data_raw;
    time = 0:sampleIntOrg*1e-6:(sampleIntOrg*1e-6*(length(animation_data)-1));
    sampleInt = sampleIntOrg;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ATTEMPT AT GETTING AMPLITUDE VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lower_threshold = -15;
upper_threshold = 20;

inds_below = find(animation_data(:,1) < lower_threshold);
inds_above = find(animation_data(:,1) > upper_threshold);

i = 1;
inds_at_change = [];

while i < length(inds_below)
%     fprintf('\n\nBEFORE: inds_below(i): %d   inds_above(next_val): %d    i: %d\n', inds_below(i), inds_above(next_val), i);
    inds_at_change = [inds_at_change; inds_below(i)];
    next_val = find(inds_above > inds_below(i), 1);
    if isempty(next_val)
        break
    else
        i = find(inds_below > inds_above(next_val), 1);
        if isempty(i)
            break
        end
    end
%     fprintf('AFTER:inds_below(i): %d   inds_above(next_val): %d    i: %d\n', inds_below(i), inds_above(next_val), i);
end

inds_at_change = inds_at_change - 50;

%now take the mean of 10ms before each of the points

inds_at_change = inds_at_change(2:end-1);

mean_at_inds_at_change = [];

for i=1:length(inds_at_change)
    mean_at_inds_at_change(i) = mean(animation_data((inds_at_change(i)-(round(0.01/(sampleInt*1e-6)))):inds_at_change(i),1));
end

%REMOVE OUTLIERS, FOR NOW ANYTHING BELOW 6pA

figure(4);subplot(2,1,1); plot(mean_at_inds_at_change,'-k.')
inds_to_remove = find(mean_at_inds_at_change < 6 | mean_at_inds_at_change > 7.5);

disp('---------NUMBER OF POINTS REMOVED-------------')
disp(num2str(length(inds_to_remove)))
disp('----------------------------------------------')
%UNTIL CODE CAN BE CHANGE TO TRULY REMOVE OUTLIERS FROM STEP ANALYSIS WILL
%WILL FUDGE IT BY AVERAGE THE POINT BEFORE AND AFTER THE OUTLIER AND
%REPLACE THE OUTLIER WITH THE RESULTING AMPLITUDE

mean_at_inds_at_change(inds_to_remove) = mean([mean_at_inds_at_change(inds_to_remove-1);mean_at_inds_at_change(inds_to_remove+1)]);

subplot(2,1,2); plot(mean_at_inds_at_change,'-k.')
start_cur_axis = axis(acurrent);
start_volt_axis = axis(avoltage);

end_cur_axis = [0 3.5 5.5 8];
end_volt_axis = [0 3.5 -50 100];

%Take 3 seconds to zoom in

%CALCULATE ZOOM

seconds_of_init_zoom = 3;
trans_cur_axis = [];
trans_volt_axis = [];
for i=1:4
    trans_cur_axis(:,i) = start_cur_axis(i)-((tanh(linspace(-3,3,seconds_of_init_zoom*30)).*(1/2) + 0.5).*(start_cur_axis(i)-end_cur_axis(i)));
    trans_volt_axis(:,i) = start_volt_axis(i)-((tanh(linspace(-3,3,seconds_of_init_zoom*30)).*(1/2) + 0.5).*(start_volt_axis(i)-end_volt_axis(i)));
end

trans_cur_axis = [start_cur_axis; trans_cur_axis; end_cur_axis];
trans_volt_axis = [start_volt_axis; trans_volt_axis; end_volt_axis];

x=pi:0.05:2*pi;
y = (1:length(x)).*(1/30);
z = linspace(0,1,63);
offset = (((1-z).*(cos(x)+1)) + (z.*y));

%2s through start adding in the axis that we move us forward
%all were doing is going to be shifting the center of the axis as fast as
%real-time is (1/30) of a second

% offset = linspace(0,1/30,2*30);

%now we'll just extend the offset to the beginning of time

% offset = [zeros(1,length(trans_cur_axis)-length(offset)) offset];

%now we have to adjust trans_cur_axis & trans_volt_axis, just extending them

% trans_cur_axis(:,1:2) = trans_cur_axis(:,1:2) + [offset' offset'];
% trans_volt_axis(:,1:2) = trans_volt_axis(:,1:2) + [offset' offset'];

%now we need to shift the x axis for the rest of the time, another 7
%seconds

total_time = 22;

trans_cur_axis = [trans_cur_axis; repmat(trans_cur_axis(end,:),total_time*30,1)];
trans_volt_axis = [trans_volt_axis; repmat(trans_volt_axis(end,:),total_time*30,1)];

% trans_cur_axis(93:end,1:2) = trans_cur_axis(93:end,1:2) + [offset(end).*(1:length(trans_cur_axis(93:end,:)))' offset(end).*(1:length(trans_cur_axis(93:end,:)))'];
% trans_volt_axis(93:end,1:2) = trans_volt_axis(93:end,1:2) + [offset(end).*(1:length(trans_volt_axis(93:end,:)))' offset(end).*(1:length(trans_volt_axis(93:end,:)))'];

%This should change unless we change our zoom rate

incorp_trans_beg = 80;
incorp_trans_end = incorp_trans_beg + 62;

trans_cur_axis(incorp_trans_beg:incorp_trans_end,1:2) = trans_cur_axis(incorp_trans_beg:incorp_trans_end,1:2) + [offset' offset'];
trans_volt_axis(incorp_trans_beg:incorp_trans_end,1:2) = trans_volt_axis(incorp_trans_beg:incorp_trans_end,1:2) + [offset' offset'];

trans_cur_axis((incorp_trans_end+1):end,:) = repmat(trans_cur_axis(incorp_trans_end,:), length(trans_cur_axis((incorp_trans_end+1):end,:)),1);
trans_volt_axis((incorp_trans_end+1):end,:) = repmat(trans_volt_axis(incorp_trans_end,:), length(trans_cur_axis((incorp_trans_end+1):end,:)),1);

trans_cur_axis((incorp_trans_end+1):end,1:2) = trans_cur_axis((incorp_trans_end+1):end,1:2) + [(1/30).*(1:length(trans_cur_axis((incorp_trans_end+1):end,:)))' (1/30).*(1:length(trans_cur_axis((incorp_trans_end+1):end,:)))'];
trans_volt_axis((incorp_trans_end+1):end,1:2) = trans_volt_axis((incorp_trans_end+1):end,1:2) + [(1/30).*(1:length(trans_volt_axis((incorp_trans_end+1):end,:)))' (1/30).*(1:length(trans_volt_axis((incorp_trans_end+1):end,:)))'];

%Will add a small zoom out at 15 seconds to 18 seconds, only to the y axis
%though
seconds_of_init_zoom = 3;
start_time_of_zoom = 6;
start_cur_axis = trans_cur_axis(start_time_of_zoom*30,3:4);
end_cur_axis = [5.5 13];

points_for_zoom = length(tanh(linspace(-3,3,seconds_of_init_zoom*30)))-1;

for i=1:2
    trans_cur_axis((start_time_of_zoom*30):(start_time_of_zoom*30+points_for_zoom),i+2) = start_cur_axis(i)-((tanh(linspace(-3,3,seconds_of_init_zoom*30)).*(1/2) + 0.5).*(start_cur_axis(i)-end_cur_axis(i)));
    trans_volt_axis((start_time_of_zoom*30):(start_time_of_zoom*30+points_for_zoom),i+2) = start_volt_axis(i)-((tanh(linspace(-3,3,seconds_of_init_zoom*30)).*(1/2) + 0.5).*(start_volt_axis(i)-end_volt_axis(i)));
end

disp(size([ones(length(trans_cur_axis((15*30+points_for_zoom+1):end,1)),1)*trans_cur_axis((15*30+points_for_zoom),3) ones(length(trans_cur_axis((15*30+points_for_zoom+1):end,1)),1)*trans_cur_axis((15*30+points_for_zoom),4)]))
disp(size(trans_cur_axis((15*30+points_for_zoom+1):end,3:4)))

% trans_cur_axis((15*30+points_for_zoom+1):end,3:4) =
% [ones(length(trans_cur_axis((15*30+points_for_zoom+1):end,1)),1)*trans_cur_axis((15*30+points_for_zoom),3); ones(length(trans_cur_axis((15*30+points_for_zoom+1):end,1)),1)*trans_cur_axis((15*30+points_for_zoom),4)]; 
trans_cur_axis((start_time_of_zoom*30+points_for_zoom+1):end,3) = ones(length(trans_cur_axis((start_time_of_zoom*30+points_for_zoom+1):end,1)),1)*trans_cur_axis((start_time_of_zoom*30+points_for_zoom),3);
trans_cur_axis((start_time_of_zoom*30+points_for_zoom+1):end,4) = ones(length(trans_cur_axis((start_time_of_zoom*30+points_for_zoom+1):end,1)),1)*trans_cur_axis((start_time_of_zoom*30+points_for_zoom),4);
points_to_plot = round((1/30)/(sampleInt*1e-6));

% Custom multicolored line
% number_of_steps = 7;
number_of_steps = 6;
minimum_step_length = 9;
inds_of_steps = ChisquaredFit_Custom_Multi(mean_at_inds_at_change,number_of_steps, minimum_step_length);

inds_of_steps_for_plot = inds_at_change(inds_of_steps);

% z = .05;
% 
% datay = mean_at_inds_at_change;
% datax = time(inds_at_change); %1:(length(test_data));
% 
% vert_inds = [1; 1; sort([[2:length(datax)]'; [2:length(datax)]'; [2:length(datax)]'; [2:length(datax)]'])];
% 
% vert_inds(end-1:end,:) = [];
% % 
% % taspectratio = get(gca, 'DataAspectRatio');
% % dxdyratio = (taxis(2)-taxis(1))/(taxis(4)-taxis(3));
% 
% % zadd = z.*ones(1,length(datax)-1);
% % zadd = zadd + z*dxdyratio*abs(sin(atan((datay(2:end)-datay(1:(end-1)))./(datax(2:end)-datax(1:(end-1))))));
% % disp(atand(((datay(2:end)-datay(1:(end-1))).*taspectratio(2))./((datax(2:end)-datax(1:(end-1))).*taspectratio(1))))\
% % angofpoints = atand(((datay(2:end)-datay(1:(end-1))))./((datax(2:end)-datax(1:(end-1)))./taspectratio(1)))
% dy = (z./2).*cos(atan(((datay(2:end)-datay(1:(end-1))))./((datax(2:end)-datax(1:(end-1))))));
% dx = (z./2).*sin(atan(((datay(2:end)-datay(1:(end-1))))./((datax(2:end)-datax(1:(end-1))))));
% 
% dx = reshape(repmat(dx,4,1),length(dx).*4,1);
% dy = reshape(repmat(dy,4,1),length(dy).*4,1);
%     
% tempmult = ones(length(dx),1);
% tempmult1 = tempmult;
% tempmult2 = tempmult;
% tempmult1(2:2:end) = tempmult(2:2:end).*-1;
% tempmult2(1:2:(end-1)) = tempmult(1:2:(end-1)).*-1;
% dx = dx.*tempmult1;
% dy = dy.*tempmult2;
% 
% vertices = [datax(vert_inds)'+dx datay(vert_inds)'+dy];
% % vertices = [datax(1:(end-1))+dx datay(1:(end-1))-dy; datax(1:(end-1))-dx datay(1:(end-1))+dy; datax(2:end)+dx datay(2:end)-dy; datax(2:end)-dx datay(2:end)+dy];
% faces = repmat([1 2 3; 2 3 4], length(datay)-1,1) + repmat(sort([[0:4:4*(length(datay)-2)]'; [0:4:4*(length(datay)-2)]']),1,3);
% 
% cmaptemp = colormap('lines');
% 
% cmapt=zeros(length(mean_at_inds_at_change)-1,3);
% 
% for i=1:length(inds_of_steps)-1
%     cmapt(inds_of_steps(i):(inds_of_steps(i+1)-1),:) = repmat(cmaptemp(i,:),length(inds_of_steps(i):(inds_of_steps(i+1)-1)),1);
% end
% 
% cmap = zeros((length(mean_at_inds_at_change)-1)*2,3);
% 
% cmap(1:2:end-1) = cmapt;
% cmap(2:2:end) = cmapt;

% seconds_to_plot = length(trans_cur_axis)/30;

axes(acurrent)
hplotc = plot(acurrent, time(1:1), animation_data(1:1, 1), '-k');% hold on;
hplotm1 = hplotc;

%         text('Units', 'normalized', 'Position', [0.9 0.9 0], 'String', '1x','FontSize',24, 'Color', [1 0 0])
% htext = text('Units', 'pixels', 'Position', [600 300 0], 'String', '1x','FontSize',24, 'Color', [1 0 0]);

axes(avoltage)
hplotv = plot(avoltage, time(1:1), animation_data(1:1, 2), '-k');

% j = 1:points_to_plot:length(trans_cur_axis)*points_to_plot;
j = 1:points_to_plot:length(animation_data);

number_of_steps_to_plot = 0;

for i=1:length(j)
    
    delete(hplotc)
    delete(hplotv)
    if ishandle(hplotm1)
        for k=1:number_of_steps_to_plot
            eval(['delete(hplotm' num2str(k) ')'])
        end
    end
    
    axes(acurrent)
    
    beginning_time_point = round(trans_cur_axis(i,1)/(sampleInt*1e-6));
    
    if beginning_time_point <= 0
        beginning_time_point = 1;
    end
    hplotc = plot(acurrent, time(beginning_time_point:j(i)), animation_data(beginning_time_point:j(i), 1), 'LineStyle', '-', 'Color', [0.5 0.5 0.5]);
    %         text('Units', 'normalized', 'Position', [0.9 0.9 0], 'String', '1x','FontSize',24, 'Color', [1 0 0])
    htext = text('Units', 'pixels', 'Position', [600 300 0], 'String', '1x','FontSize',24, 'Color', [1 0 0]);

    if (j(i) > inds_at_change(1))
        hold on
        %         tfaces = faces(find(inds_at_change >= beginning_time_point,1)*2:find(inds_at_change < j(i),1,'last')*2,:);
        %         tcmap = cmap(find(inds_at_change >= beginning_time_point,1)*2:find(inds_at_change < j(i),1,'last')*2,:);
        %         phandle =
        %         patch('Faces',tfaces,'Vertices',vertices,'FaceVertexCData',tcmap,'FaceColor','flat', 'EdgeColor', 'none', 'FaceAlpha', 1);


        %         inds_of_steps
        %         inds_at_change
        %range of plot (relative to the inds_at_change vector)
        range_of_plot = inds_at_change(find(inds_at_change >= beginning_time_point,1):find(inds_at_change < j(i),1,'last'));
        range_of_plot_inds = find(inds_at_change >= beginning_time_point,1):find(inds_at_change < j(i),1,'last');

        if ~isempty(range_of_plot)
            starting_step_num = find(inds_of_steps_for_plot <= range_of_plot(1),1, 'Last')
            ending_step_num = find(inds_of_steps_for_plot <= range_of_plot(end),1,'Last')
            colormap_to_use = colormap('lines');
            colormap_to_use(7,:) = [];
            if ~isempty(starting_step_num)
                number_of_steps_to_plot = (ending_step_num-starting_step_num)+1

                target_start_point = range_of_plot(1);%inds_at_change(find(inds_at_change >= beginning_time_point,1));
                target_end_point = range_of_plot(end);%inds_at_change(find(inds_at_change < j(i),1,'last'));


                for k=1:number_of_steps_to_plot

                    if k==1
                        %the first line is the only line, start to end
                        eval(['hplotm' num2str(k) ' = plot(acurrent, time(range_of_plot), mean_at_inds_at_change(range_of_plot_inds), ''LineStyle'', ''-'', ''Color'',colormap_to_use(starting_step_num+(k-1),:), ''LineWidth'', 3)']);
                    else
                        eval(['hplotm' num2str(k) ' = plot(acurrent, time(inds_at_change(inds_of_steps(starting_step_num+(k-1)):range_of_plot_inds(end))), mean_at_inds_at_change(inds_of_steps(starting_step_num+(k-1)):range_of_plot_inds(end)), ''LineStyle'', ''-'', ''Color'',colormap_to_use(starting_step_num+(k-1),:), ''LineWidth'', 3)']);
                    end
                end
            end
        end
    end

    axes(avoltage)
    hplotv = plot(avoltage, time(beginning_time_point:j(i)), animation_data(beginning_time_point:j(i), 2), 'LineStyle', '-', 'Color', [0.5 0.5 0.5]);
    
    axis(acurrent, trans_cur_axis(i,:));
    axis(avoltage, trans_volt_axis(i,:));

    if generate_movie == 1
        F = getframe(gcf);
        mov = addframe(mov, F);
    else
        pause(0.01)
    end
end

%%%%%%%%%%%%%%%%%%
%Pause for 5 seconds
%%%%%%%%%%%%%%%%%%%%

amount_of_pause = 5;
for m=1:amount_of_pause*30
    if generate_movie == 1
        F = getframe(gcf);
        mov = addframe(mov, F);
    else
        pause(0.01)
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FADE OUT OLD TRACE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotchange = zeros(75,3);
%     textchange(:,1) = cos(linspace(0,2*pi,152))./2+.5;
plotchange(:,1) = cos(linspace((3*pi)/2,2*pi,75))./2+.5;
plotchange(:,2) = cos(linspace((3*pi)/2,2*pi,75))./2+.5;
plotchange(:,3) = cos(linspace((3*pi)/2,2*pi,75))./2+.5;

if ishandle(hplotm1)
    for k=1:number_of_steps_to_plot
        eval(['plotchange' num2str(k) '=zeros(75,3);'])
        eval(['plotchange' num2str(k) '(1,:)= get(hplotm' num2str(k) ',''Color'');'])
        eval(['plotchange' num2str(k) '(:,1) = linspace(plotchange' num2str(k) '(1,1), 1, 75);'])
        eval(['plotchange' num2str(k) '(:,2) = linspace(plotchange' num2str(k) '(1,2), 1, 75);'])
        eval(['plotchange' num2str(k) '(:,3) = linspace(plotchange' num2str(k) '(1,3), 1, 75);'])
    end
end

for i=1:length(plotchange)
    set(hplotc, 'Color', plotchange(i,:))
    set(hplotv, 'Color', plotchange(i,:))

    if ishandle(hplotm1)
        for k=1:number_of_steps_to_plot
            eval(['set(hplotm' num2str(k) ', ''Color'', plotchange(i,:))'])
        end
    end
    
    if generate_movie == 1
        F = getframe(gcf);
        mov = addframe(mov, F);
    else
        pause(0.01)
    end
end

delete(hplotc)
delete(hplotv)

mov = close(mov);