%BILL-makes sure everything is in the same directory and this should work.
%I just pulled out the path part of the ABF loading so it should search the
%current directory

clear all
close all

%set some parameters
maxTransientLength = 7e-3; %take full length of event's worth of following transient up to this long
saturatedSampleTime = 75;  %number of samples from end of event to when amp is out of saturation
showPlots = 0;
maxIterations = 750;

firstRun = 1;


%initial guess for fit
x0 = [50; -1e3; 50; -1e2; 50; -10; 0];                

%define abf file data variables outside loop for persistance between cycles
traceData = [];
sampleInt = 0;
numEventsTot = 0;

% [ABFfilename, ABFpathname] = uigetfile('*.abf','Open ABF directory', [pwd '/'], 'MultiSelect', 'on');
% ABFpathname
ABFpathname = 'U:\pCLAMP\Station X Data\';

%make sure you have the ABF file listed on next line
[filename, pathname] = uigetfile('*.mat','Open Trace Data', [pwd '/'], 'MultiSelect', 'on');

if isequal(filename,0) || isequal(pathname,0)
    %do nothing, user pressed cancel
    
    poreEventData.defaultMATPathname = [pwd '\'];
else
    %update file info
    poreEventData.defaultMATPathname = pathname;
    %handles.MATpathname = pathname;

    %check if multiple MAT files selected
    if( iscell( filename ) )
        numMATfiles = length(filename);
        %only one file selected
    else
        numMATfiles = 1;
    end

    for currentMATindex = 1:numMATfiles
        if( numMATfiles > 1 )
            load( [pathname filename{currentMATindex}] );
        else
            load( [pathname filename] );
        end
        
        %pre-allocate memory for exponential coefs
        expCoefs = zeros(length(poreEventData.detectedEvents(:,1)), length(x0));

        %loop through events detected in mat file
        for currentEvent = 1:length(poreEventData.detectedEvents(:,1))

            %only perform compensation on fishing events
            if( poreEventData.detectedEvents(currentEvent,8) == 3 &&  ...
                    poreEventData.detectedEvents(currentEvent,4) == poreEventData.detectedEvents(currentEvent,5))

                %load current sweep and current file
                %     fprintf('prev. sweep: %d, current sweep: %d\n', poreEventData.currentSweep, poreEventData.detectedEvents(currentEvent,1) );
                %     loadFile = ( poreEventData.currentSweep ~=
                %     poreEventData.detectedEvents(currentEvent,1) );
                poreEventData.currentSweep = poreEventData.detectedEvents(currentEvent,1);
                poreEventData.currentFileIndex = poreEventData.detectedEvents(currentEvent,10);

                %load data file
                if( 1 ) %loadFile )
                    if( iscell(poreEventData.filename) )
                        [traceData, sampleInt, numEventsTot] = import_abf( ...
                            [ABFpathname poreEventData.filename{poreEventData.currentFileIndex}], ...
                            poreEventData.currentSweep );
                        %             [traceData, sampleInt, numEventsTot] = import_abf( ...
                        %                 [poreEventData.filename{poreEventData.currentFileIndex}], ...
                        %                 poreEventData.currentSweep );
                    else
                        [traceData, sampleInt, numEventsTot] = import_abf( [poreEventData.pathname poreEventData.filename], poreEventData.currentSweep );
                        fprintf('should never go here\n');
                    end
                end

                traceData = double(traceData);

                %convert sample interval to microseconds
                sampleInt = double(sampleInt/1e6);

                eventData = traceData( ...
                    poreEventData.detectedEvents(currentEvent,2)+saturatedSampleTime:poreEventData.detectedEvents(currentEvent,3), 1 );

                %create time vector
                time = 0:sampleInt:(length(eventData)-1)*sampleInt;
                
                %estimate transient length for first event
                if( isfield(poreEventData, 'maxTransientLength') )
                    maxTransientLength = poreEventData.maxTransientLength
                    saturatedSampleTime = poreEventData.saturatedSampleTime
                else
                    if( firstRun == 1 && time(end) > 35e-3 )
                        maxTransientLength = 1e-3;

                        while( (maxTransientLength < (time(end)-5e-3)) && (var( eventData( (floor(maxTransientLength/sampleInt)):(floor(maxTransientLength/sampleInt)+floor(5e-3/sampleInt)) ) ) > var( eventData( (end-floor(5e-3/sampleInt)):end ) )) )

                            maxTransientLength = maxTransientLength + 1e-3;
                        end

                        if( maxTransientLength > 7e-3 )
                            maxTransientLength = maxTransientLength*2.5
                        else
                            maxTransientLength = 7e-3;
                        end

                        firstRun = 0;
                    end
                end


                %take the post-event transient (length of event up to maxTransientLength)
                transientEndIndex = (2*poreEventData.detectedEvents(currentEvent,3))-poreEventData.detectedEvents(currentEvent,2);
                if( (transientEndIndex-(poreEventData.detectedEvents(currentEvent,3)+1)) > (maxTransientLength/sampleInt) )
                    transientEndIndex = (maxTransientLength/sampleInt) + poreEventData.detectedEvents(currentEvent,3) + 1;
                end


                %get transient from after fishing event
                transStartIndex = (poreEventData.detectedEvents(currentEvent,3)+saturatedSampleTime+1+7);
                if( length(traceData) > transientEndIndex+1 )
                    transientData = traceData( transStartIndex:(transientEndIndex+1), 1 );
                else
                    transientData = traceData( transStartIndex:end, 1 );
                    fprintf('event near end of sweep: %d\n',transStartIndex)
                end

                %Make sure it runs long enough to get a good fit 'display', 'iter',
                options = optimset('maxFunEvals', 1e16, 'maxIter', maxIterations);

                %Check out documentation, lsqcurvefit is just an application of lsqnonlin
                transientData = double(transientData);
                
                if( isfield(poreEventData, 'expCoefs') )
                    expFit = poreEventData.expCoefs(currentEvent, :);
                else
                    [expFit, resnorm, RESIDUAL, EXITFLAG] = lsqcurvefit( ...
                        @(x,t) exponentialFit(x,time(1:length(transientData))), x0, double( time(1:length(transientData)) ), transientData', [], [], options);
%                     disp(EXITFLAG)

                    %save fit coefs
                    expCoefs(currentEvent,:) = expFit;
                end
                
                %update inital guess to previous fit
                x0 = expFit;
                

                %add inverse transient to compensate event
                fitData = exponentialFit( expFit,time(1:length(transientData)) );
                
%                 if( abs( expFit(1) ) < 5 )
                    compensatedEvent = [( fitData - fitData(end) )'+eventData(1:length(transientData)); eventData(length(transientData)+1:end)];
%                 else
%                     compensatedEvent = [( fitData )'+eventData(1:length(transientData)); eventData(length(transientData)+1:end)];
%                 end

                %calc mean
                meanAmplitude = mean( compensatedEvent);

                %%%%%%%%%%%%%%%%%%%%
                %plot results
                if( showPlots == 1 )
                    expFit
                    plot( time , eventData, 'r'); hold on, grid on
                    % plot( [length(eventData):(length(eventData)+length(transientData)-1)], transientData, 'b');
                    plot( time(1:length(transientData)), transientData, 'b');
                    %plot(b, trans_decimate, 'go', time(1:length(transientData)), transSpline, 'm')
                    plot( time(1:length(transientData)), exponentialFit(expFit,time(1:length(transientData))), 'r');

                    plot( time, compensatedEvent, 'c');
                    % axis([0 ((length(eventData)+length(transientData)-1)) -100 100]);
                    axis([0 time(end) -100 100]);

                    plot([0, length(eventData)-1], [meanAmplitude, meanAmplitude], 'g--');

                    fprintf('eventData[%.4f:%.4f]s, transientData[%.4f,%.4f]s\n', ...
                        poreEventData.detectedEvents(currentEvent,2)*sampleInt, ...
                        poreEventData.detectedEvents(currentEvent,3)*sampleInt, ...
                        (poreEventData.detectedEvents(currentEvent,3)+1)*sampleInt, ...
                        transientEndIndex*sampleInt );
                    pause
                    hold off
                end
                %%%%%%%%%%%%%%%%%%%%%%

                %save amplitude back to poreEventData
                poreEventData.detectedEvents_ms(currentEvent, 5) = meanAmplitude;

            end

            fprintf('current event: %d/%d\n', currentEvent, length(poreEventData.detectedEvents(:,1)) );
        end

        
        %add parameters to mat file
        poreEventData.maxTransientLength = maxTransientLength;
        poreEventData.saturatedSampleTime = saturatedSampleTime;
        poreEventData.expCoefs = expCoefs;
        
        %save mat file to disk        
        fprintf('saving data...\n');
        if( numMATfiles > 1 )
            save_filename = [pathname filename{currentMATindex}(1:end-4) '_spline.mat']
        else
            save_filename = [pathname filename(1:end-4) '_spline.mat']
        end
        
        save( save_filename, 'poreEventData' );

    end
end
