clear all;

generate_movie = 1;

% filename = 'Z:\pClamp\Station BN Data\09519009_Filt.abf';
filename = 'C:\Users\BLGyarfas\Documents\MATLAB\datafiles\09519009_Filt.abf';
video_output_name = 'FastVid3Step_With_Zoom2_Raw.avi';
animation_data = [];

start_time = 3.5;
end_time = 3.3;

[data, sampleInt] = import_abf(filename, 32);

animation_data = [animation_data; data(round(start_time/(sampleInt*1e-6)):end,:)];

[data, sampleInt] = import_abf(filename, 33);

animation_data = [animation_data; data(1:round(end_time/(sampleInt*1e-6)),:)];

time = 0:sampleInt*1e-6:(sampleInt*1e-6*(length(animation_data)-1));

h1 = figure('Position',[200 50 720 480]);
set(h1, 'DoubleBuffer', 'on', 'Color', [1 1 1]);

if generate_movie == 1
    mov = avifile(video_output_name,'Compression', 'none', 'fps', 30, 'quality', 100);
end

%trying for 30 fps so plot every 33.3 ms

frame_per_sec = 30;
points_to_plot = round((1/30)/(sampleInt*1e-6));

seconds_to_plot = time(end);

frame_counter = 1;
with_voltage = 1;

if with_voltage == 0
    %Without Voltage
    for i=1:points_to_plot:round(seconds_to_plot/(sampleInt*1e-6))
        h = plot(time(1:i), animation_data(1:i, 1));
        axis([0 5 -20 30]);
        xlabel('time (s)')
        ylabel('pA')
        tempInset = get(gca, 'TightInset');
        tempPos = get(gca, 'Position');
        tempOuter = get(gca, 'OuterPosition');
        set(gca, 'Position', [tempInset(1) tempInset(2) 1-tempInset(3)-tempInset(1) 1-tempInset(4)-tempInset(2)])
        if generate_movie == 1
            F = getframe(gcf);
            mov = addframe(mov, F);
        end
    end
else
    %with voltage
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %FIRST PART - Real Time
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    acurrent = subplot(2,1,1);
    axis([0 5 -20 65]);
    ylabel('pA')
    set(acurrent, 'XTickLabel', '', 'NextPlot', 'replacechildren');
    set(acurrent, 'XLimMode', 'manual', 'YLimMode', 'manual');
    tempInsetCurrent = get(acurrent, 'TightInset');

    avoltage = subplot(2,1,2);
    axis([0 5 -50 200]);
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
    hplotc = plot(acurrent, time(1:i), animation_data(1:i, 1), '-k');
    %         text('Units', 'normalized', 'Position', [0.9 0.9 0], 'String', '1x','FontSize',24, 'Color', [1 0 0])
    htext = text('Units', 'pixels', 'Position', [600 300 0], 'String', '1x','FontSize',24, 'Color', [1 0 0]);

    axes(avoltage)
    hplotv = plot(avoltage, time(1:i), animation_data(1:i, 2), '-k');

    for i=1:points_to_plot:round(seconds_to_plot/(sampleInt*1e-6))

        delete(hplotc)
        delete(htext)
        delete(hplotv)

        axes(acurrent)
        hplotc = plot(acurrent, time(1:i), animation_data(1:i, 1), '-k');
        %         text('Units', 'normalized', 'Position', [0.9 0.9 0], 'String', '1x','FontSize',24, 'Color', [1 0 0])
        htext = text('Units', 'pixels', 'Position', [600 300 0], 'String', '1x','FontSize',24, 'Color', [1 0 0]);

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
    %%%%% HIGHLIGHT ZOOM AREA - FIRST ZOOM
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    hold on;
    %Highlight Zoom, blink the line for 4 seconds
    x = linspace(0,5*pi, 150); % gives 120 values
    y = cos(x)/2 + 0.5;

    axes(acurrent)

    linehandle1 = line([2.75 2.75],[-20 65], [0 0], 'LineStyle', '--', 'Color', [y(1) y(1) y(1)], 'LineWidth', 2);
    linehandle2 = line([4.25 4.25],[-20 65], [0 0], 'LineStyle', '--', 'Color', [y(1) y(1) y(1)], 'LineWidth', 2);

    for i=2:length(y)

        set(linehandle1, 'Color', [y(i) y(i) y(i)]);
        set(linehandle2, 'Color', [y(i) y(i) y(i)]);

        if generate_movie == 1
            F = getframe(gcf);
            mov = addframe(mov, F);
        else
            pause(0.01)
        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ZOOM IN - FIRST PART, CHANGE THE TIME TIME INDICIATOR FOR X1 TO X10
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    start_cur_axis = axis(acurrent);
    start_volt_axis = axis(avoltage);

    end_cur_axis = [2.75 4.25 -10 20];
    end_volt_axis = [2.75 4.25 -50 100];

    %Take 4 seconds to zoom in

    %CALCULATE ZOOM

    for i=1:4
        trans_cur_axis(:,i) = start_cur_axis(i)-((tanh(linspace(-3,3,150)).*(1/2) + 0.5).*(start_cur_axis(i)-end_cur_axis(i)));
        trans_volt_axis(:,i) = start_volt_axis(i)-((tanh(linspace(-3,3,150)).*(1/2) + 0.5).*(start_volt_axis(i)-end_volt_axis(i)));
    end

    trans_cur_axis = [start_cur_axis; trans_cur_axis; end_cur_axis];
    trans_volt_axis = [start_volt_axis; trans_volt_axis; end_volt_axis];

    %CACLUATE TEXT CHANGE

    textchange = zeros(152,3);
    %     textchange(:,1) = cos(linspace(0,2*pi,152))./2+.5;
    textchange(:,1) = ones(152,1);
    textchange(:,2) = cos(linspace(pi,3*pi,152))./2+.5;
    textchange(:,3) = cos(linspace(pi,3*pi,152))./2+.5;

    for i=1:length(trans_cur_axis)
        axis(acurrent, trans_cur_axis(i,:));
        axis(avoltage, trans_volt_axis(i,:));

        set(htext, 'Color', textchange(i,:));
        if i > length(trans_cur_axis)/2
            set(htext, 'String', 'x10')
        end

        if generate_movie == 1
            F = getframe(gcf);
            mov = addframe(mov, F);
        else
            pause(0.01)
        end
    end

    delete(linehandle1)
    delete(linehandle2)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %FADE OUT OLD TRACE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    plotchange = zeros(150,3);
    %     textchange(:,1) = cos(linspace(0,2*pi,152))./2+.5;
    plotchange(:,1) = cos(linspace(pi,2*pi,150))./2+.5;
    plotchange(:,2) = cos(linspace(pi,2*pi,150))./2+.5;
    plotchange(:,3) = cos(linspace(pi,2*pi,150))./2+.5;

    for i=1:length(plotchange)
        set(hplotc, 'Color', plotchange(i,:))
        set(hplotv, 'Color', plotchange(i,:))

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
    %REPLAY AT 10X
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    slow_down_factor = 10;
    points_to_plot = round((1/30)/(sampleInt*1e-6*slow_down_factor));

    start_time = 2.75;
    end_time = 4.25;

    axes(acurrent)
    hplotc = plot(acurrent, time(1:i), animation_data(1:i, 1), '-k');
    %         text('Units', 'normalized', 'Position', [0.9 0.9 0], 'String', '1x','FontSize',24, 'Color', [1 0 0])
    htext = text('Units', 'pixels', 'Position', [600 300 0], 'String', '10x','FontSize',24, 'Color', [1 0 0]);

    axes(avoltage)
    hplotv = plot(avoltage, time(1:i), animation_data(1:i, 2), '-k');

    for i=round(start_time/(sampleInt*1e-6)):points_to_plot:round(end_time/(sampleInt*1e-6))

        delete(hplotc)
        delete(htext)
        delete(hplotv)

        axes(acurrent)
        hplotc = plot(acurrent, time(1:i), animation_data(1:i, 1), '-k');
        %         text('Units', 'normalized', 'Position', [0.9 0.9 0], 'String', '1x','FontSize',24, 'Color', [1 0 0])
        htext = text('Units', 'pixels', 'Position', [600 300 0], 'String', '10x','FontSize',24, 'Color', [1 0 0]);

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
    %%ZOOM HIGHLIGHT AGAIN - 2ND TIME
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x = linspace(0,5*pi, 150); % gives 120 values
    y = cos(x)/2 + 0.5;

    axes(acurrent)

    linehandle1 = line([3.92 3.92],[-20 65], [0 0], 'LineStyle', '--', 'Color', [y(1) y(1) y(1)], 'LineWidth', 2);
    linehandle2 = line([4.05 4.05],[-20 65], [0 0], 'LineStyle', '--', 'Color', [y(1) y(1) y(1)], 'LineWidth', 2);

    for i=2:length(y)

        set(linehandle1, 'Color', [y(i) y(i) y(i)]);
        set(linehandle2, 'Color', [y(i) y(i) y(i)]);

        if generate_movie == 1
            F = getframe(gcf);
            mov = addframe(mov, F);
        else
            pause(0.01)
        end

    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% ZOOM IN - 2ND TIME
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    start_cur_axis = axis(acurrent);
    start_volt_axis = axis(avoltage);

    end_cur_axis = [3.92 4.05 -2.5 17.5];
    end_volt_axis = [3.92 4.05 -50 100];

    %Take 4 seconds to zoom in

    %awesome zoom

    trans_cur_axis = [];
    trans_volt_axis = [];
    for i=1:4
        trans_cur_axis(:,i) = start_cur_axis(i)-((tanh(linspace(-3,3,150)).*(1/2) + 0.5).*(start_cur_axis(i)-end_cur_axis(i)));
        trans_volt_axis(:,i) = start_volt_axis(i)-((tanh(linspace(-3,3,150)).*(1/2) + 0.5).*(start_volt_axis(i)-end_volt_axis(i)));
    end

    trans_cur_axis = [start_cur_axis; trans_cur_axis; end_cur_axis];
    trans_volt_axis = [start_volt_axis; trans_volt_axis; end_volt_axis];

    %text change from x10 to x30

    textchange = zeros(152,3);
    %     textchange(:,1) = cos(linspace(0,2*pi,152))./2+.5;
    textchange(:,1) = ones(152,1);
    textchange(:,2) = cos(linspace(pi,3*pi,152))./2+.5;
    textchange(:,3) = cos(linspace(pi,3*pi,152))./2+.5;

    for i=1:length(trans_cur_axis)
        axis(acurrent, trans_cur_axis(i,:));
        axis(avoltage, trans_volt_axis(i,:));

        set(htext, 'Color', textchange(i,:));
        if i > length(trans_cur_axis)/2
            set(htext, 'String', 'x100')
        end

        if generate_movie == 1
            F = getframe(gcf);
            mov = addframe(mov, F);
        else
            pause(0.01)
        end
    end

    delete(linehandle1)
    delete(linehandle2)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %FADE OUT OLD TRACE - 2ND TIME
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    plotchange = zeros(150,3);
    %     textchange(:,1) = cos(linspace(0,2*pi,152))./2+.5;
    plotchange(:,1) = cos(linspace(pi,2*pi,150))./2+.5;
    plotchange(:,2) = cos(linspace(pi,2*pi,150))./2+.5;
    plotchange(:,3) = cos(linspace(pi,2*pi,150))./2+.5;

    for i=1:length(plotchange)
        set(hplotc, 'Color', plotchange(i,:))
        set(hplotv, 'Color', plotchange(i,:))

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
    %REPLAY AT 100X - 3RD REPLAY
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    slow_down_factor = 100;
    points_to_plot = round((1/30)/(sampleInt*1e-6*slow_down_factor));

    start_time = 3.92;
    end_time = 4.05;

    axes(acurrent)
    hplotc = plot(acurrent, time(1:i), animation_data(1:i, 1), '-k');
    %         text('Units', 'normalized', 'Position', [0.9 0.9 0], 'String', '1x','FontSize',24, 'Color', [1 0 0])
    htext = text('Units', 'pixels', 'Position', [600 300 0], 'String', '100x','FontSize',24, 'Color', [1 0 0]);

    axes(avoltage)
    hplotv = plot(avoltage, time(1:i), animation_data(1:i, 2), '-k');

    for i=round(start_time/(sampleInt*1e-6)):points_to_plot:round(end_time/(sampleInt*1e-6))

        delete(hplotc)
        delete(htext)
        delete(hplotv)

        axes(acurrent)
        hplotc = plot(acurrent, time(1:i), animation_data(1:i, 1), '-k');
        %         text('Units', 'normalized', 'Position', [0.9 0.9 0], 'String', '1x','FontSize',24, 'Color', [1 0 0])
        htext = text('Units', 'pixels', 'Position', [600 300 0], 'String', '100x','FontSize',24, 'Color', [1 0 0]);

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
    %FADE OUT OLD TRACE - 3RD TIME
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    plotchange = zeros(150,3);
    %     textchange(:,1) = cos(linspace(0,2*pi,152))./2+.5;
    plotchange(:,1) = cos(linspace(pi,2*pi,150))./2+.5;
    plotchange(:,2) = cos(linspace(pi,2*pi,150))./2+.5;
    plotchange(:,3) = cos(linspace(pi,2*pi,150))./2+.5;

    for i=1:length(plotchange)
        set(hplotc, 'Color', plotchange(i,:))
        set(hplotv, 'Color', plotchange(i,:))

        if generate_movie == 1
            F = getframe(gcf);
            mov = addframe(mov, F);
        else
            pause(0.01)
        end
    end

    delete(hplotc)
    delete(hplotv)

end

if generate_movie == 1
    mov = close(mov);
end