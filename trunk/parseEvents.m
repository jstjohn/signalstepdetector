% Parse events takes raw data in sweep form and returns
% start and stop indices
% eventTiming = [startIndex0, endIndex0; startIndex1, endIndex1;...]
% eventTiming(:,6) is the type of event, strange meaning for most, for the
% primer annealing version, the type of event is 10 for no eject, 11 for
% eject
% sweepData - [n,2] vector containing current trace data
% samplePeriod - sample period in sec
% baseline - baseline current for data
% stdDev - standard deviation for baseline current
% BIG ASSUMPTION, WERE GOING TO BE 'RUNNING' THE PORE IN POSITIVE VOLTAGE
% MODE, SO IF WERE IN NEGATIVE VOLTAGE,WERE EJECTING A SAMPLE
function [eventTiming, cutoffEventSweep, cutoffEventStartTime, currentDNAnumber, beginHighOut, cutoffVChangeStart, cutoffEventFoundRamp, cutoffEventFoundPrimer, eventVoltageTimingPrimer] = parseEvents( sweepData, samplePeriod, baseLine, stdDev, minEventTime, fishingBaseline, currentSweep, continuedEventSweep, continuedEventStartTime, currentDNAnumber, passedBeginHigh, passedVChangeStart, passedEventFoundRamp, passedEventVoltageTimingPrimer, rmsThreshold, holdingVoltage, probingVoltage, ejectVoltage )
disp(['-------------- NEW SWEEP: ' num2str(currentSweep) '----------------'])
disp(['-------------- NEW SWEEP: ' num2str(currentSweep) '----------------'])
disp(['-------------- NEW SWEEP: ' num2str(currentSweep) '----------------'])

%get size of sweepData
[numSamples, numSignals] = size(sweepData);

%flag to mark beginning for event by voltage change
%calculate beginning during transition cutoff
%only valid during fishing mode
eventStartAtVoltageChange = 1;

cutoffVChangeStart = 0;
cutoffEventFoundRamp = [0,0];
cutoffEventFoundPrimer = [nan,nan,nan];

%min voltage change to count for event cutoff on voltage change
voltageChangeEpsilon = 50;

%negative fishing voltage
fishingVoltage = -20;

%don't count events that occur with primer anneal voltage (50mV)
ignorePrimerAnneal = 1;

%min duration of detected events
if( minEventTime <= 10e-3 )
    minEventSamples = 2;
else
    minEventSamples = round( (minEventTime./1000) ./ samplePeriod );
    if(minEventSamples > numSamples)
        minEventSamples = numSamples;
    end
end

%is there a cutoff event from previous sweep?
if( continuedEventSweep ~= 0 )
    disp(['Cutoff Present! continuedEventStartTime: ' num2str(continuedEventStartTime)])
end

if( (continuedEventSweep ~= 0) && (continuedEventStartTime ~= 0) )
    cutoffEventPresent = 1;
else
    cutoffEventPresent = 0;
end

%init return variables
cutoffEventSweep = 0;
cutoffEventStartTime = 0;
beginHighOut = 0;

%event needs to be over for 0.2ms for end to be detected
minOutEventSamples = round( (0.3/1000) / samplePeriod );

%define thresholds
threshold = baseLine - rmsThreshold*stdDev
baselineThreshold = 45;
thresholdFishing = 0; %fishingBaseline;% + 7*stdDev;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if non fishing experiment
if( fishingBaseline == -250 )

    %array that thresholds the sampled data
    %1 if sample < threshold & sample > 5
    inEventArray = ( (sweepData(:,1) < threshold) & (sweepData(:,1) > -5) );
    outEventArray = ~inEventArray;

    fishingData = 0;
    rampingData = 0;
    primerData = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ramp case
elseif( fishingBaseline == -500 )
        %array that thresholds the sampled data
    %1 if sample < threshold & sample > 5
    inEventArray = ( (sweepData(:,1) < threshold) & (sweepData(:,1) > -20) );
    outEventArray = ~inEventArray;

    fishingData = 0;
    rampingData = 1;
    primerData = 0;
    %use voltage delta to find where voltage switches from capture to
    %dropped voltage, and for eject after
    
    %first find changes of >= 50mV

    voltageChangeIndices = find( abs(sweepData(1:end-2,2)-sweepData(2:end-1,2)) >= 50);
    voltageChangeIndices = voltageChangeIndices+2; % line up voltage index
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %primer annealing case
elseif( fishingBaseline == -750 )
    %array that thresholds the sampled data
    %1 if sample < threshold & sample > 5
%     inEventArray = ( (sweepData(:,1) < threshold)) & (sweepData(:,2) < (holdingVoltage+2.5) & ( sweepData(:,2) > probingVoltage - 2.5) );

    %We'll find the start/end time for those events that actually have
    %voltage changes in them
    inEventArray = ( (sweepData(:,1) < threshold) );
    outEventArray = ~inEventArray;

    fishingData = 0;
    rampingData = 0;
    primerData = 1;
    %use voltage delta to find where voltage switches from capture to
    %dropped voltage, and for eject after
    
    %first find changes of >= 50mV

    voltageChangeIndices = find( abs(sweepData(1:end-2,2)-sweepData(2:end-1,2)) >= 30);
    voltageChangeIndices = voltageChangeIndices+2; % line up voltage index
    
else
    %set flag for fishing data
    fishingData = 1;
    rampingData = 0;
    primerData = 0;
    %find last occurence of primer anneal voltage in sweep
    primerAnnealVoltage = 50;

    %array that thresholds the sampled data
    %1 if sample < threshold & sample > -5
    if( numSignals == 2 && ignorePrimerAnneal == 1 )
        %doesn't include primer anneal voltage in event detection
        inEventArray = ( ((sweepData(:,1) > thresholdFishing) & (sweepData(:,1) < threshold)) & ((sweepData(:,2) > (primerAnnealVoltage+2.5)) | (sweepData(:,2) < (primerAnnealVoltage-2.5))) & ( sweepData(:,2) > fishingVoltage + 2.5) );
        outEventArray = ~inEventArray;
    else
        inEventArray = ( (sweepData(:,1) > thresholdFishing) & (sweepData(:,1) < threshold) );
        outEventArray = ~inEventArray;
    end

    %use voltage delta to find where voltage switches to and from 50mV
    %first find changes of >= 35mV

    voltageChangeIndices = find( abs(sweepData(1:end-2,2)-sweepData(2:end-1,2)) >= 35);
    voltageChangeIndices = voltageChangeIndices+2; % line up voltage index

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%event detection kernel
%an event is detected if sequence of consecutive 1's is found in
%inEventArray of length minEventSamples
inEventKernel = ones(minEventSamples,1);
outEventKernel = ones(minOutEventSamples,1);

%filter data to find event start/stop indices
inEventArrayFiltered = filter(inEventKernel, 1, inEventArray);
outEventArrayFiltered = filter(outEventKernel, 1, outEventArray);

%reverse search for min event return time?

%find where events begin and end
eventTimingStart = find( inEventArrayFiltered == (minEventSamples-1) );
eventTimingEnd = find( outEventArrayFiltered == (minOutEventSamples-1) );

%loop through found event transitions and pair start/stop times
i=1;        %eventTimingStart index
j=1;        %eventTimingEnd index
count=0;    %keep track of number of events found
cutoffEventFound = 0;   %flag is set if start but not end is found

%pre-allocate eventTiming matrix
eventTiming = zeros( max([length(eventTimingStart),length(eventTimingEnd)]), 4 );

while( (i <= length(eventTimingStart)) && (j <= length(eventTimingEnd)) )
    %find event end past event start
    j = find( eventTimingEnd > eventTimingStart(i), 1 );
    if(isempty(j))
        fprintf('Break out due to no more end events\n');
        cutoffEventFound = 1;
        break;
    end   %break out of loop if no more event ends

    %extract start/stop times from find result
    eventTiming(count+1,:) = [ eventTimingStart(i)-(minEventSamples-2), eventTimingEnd(j)-(minOutEventSamples-1), currentSweep, currentSweep ];
    count = count + 1;

    %find next event start
    i = find( eventTimingStart > eventTimingEnd(j), 1 );
    if(isempty(i))
        %         fprintf('Break out due to no more start events\n');
        break;
    end   %break out of loop if no more event starts
end

disp('AFTER FIRST LOOP, eventTiming')
temp = eventTiming;
temp(:,1:2) = temp(:,1:2).*samplePeriod;
disp(temp)
disp('-------')
if( cutoffEventPresent == 1 && isempty(eventTimingEnd) )
    %event continues through sweep, pass cutoff event info through
    disp('Thinks that it is a passthrough event')
    cutoffEventSweep = continuedEventSweep;
    cutoffEventStartTime = continuedEventStartTime;
end

%check and see if unfinished event is at end of eventTimingStart
if( ~isempty(eventTimingStart) )   %there is an event start
    if( cutoffEventFound ) %and it doesn't have a matching end

        disp('Getting the info together')
        %pass info for next round of parseEvents
        cutoffEventSweep = currentSweep;
        cutoffEventStartTime = eventTimingStart(i)-(minEventSamples-2);

        %remove start transition for cutoff event
        maxTransitionLength = round( (0.1/1000) / samplePeriod );
        if (length(sweepData) < cutoffEventStartTime+maxTransitionLength)
            maxTransitionLength = length(sweepData) - cutoffEventStartTime;
        end
        %for each event, loop back until a sharp change in current is spotted
        if( ~rampingData )
            for j = 0:maxTransitionLength
                currentDelta = abs( sweepData(cutoffEventStartTime+j,1) - sweepData(cutoffEventStartTime+maxTransitionLength, 1) );
                if( currentDelta <= stdDev )
                    %update cutoff event start time
                    cutoffEventStartTime = cutoffEventStartTime+j;
                    break
                end
            end
        end
        
        if( primerData == 1)
            %Need to pass through the voltage inds if there are any
            vChangeHold = find(sweepData(cutoffEventStartTime:end,2) > holdingVoltage - 2 & sweepData(cutoffEventStartTime:end,2) < holdingVoltage + 2,1);
            vChangeProbe = find(sweepData(cutoffEventStartTime:end,2) > probingVoltage - 2 & sweepData(cutoffEventStartTime:end,2) < probingVoltage + 2,1);
            vChangeEject = find(sweepData(cutoffEventStartTime:end,2) > ejectVoltage - 2 & sweepData(cutoffEventStartTime:end,2) < ejectVoltage + 2,1);

            if( ~isempty(vChangeHold) ), cutoffEventFoundPrimer(1) = vChangeHold; end
            if( ~isempty(vChangeProbe) ), cutoffEventFoundPrimer(2) = vChangeProbe; end
            if( ~isempty(vChangeEject) ), cutoffEventFoundPrimer(3) = vChangeEject; end

        end

    end
end

%trim off un-used space due to pre-allocating
if( isequal(count, 0) )
    eventTiming = [];
else    
    eventTiming = [eventTiming(1:count,:), zeros(count,4)];
end

% Initialize our eventVoltageTimingPrimer

eventVoltageTimingPrimer = zeros(count, 3);

% disp('AFTER TRIMMING TRANSIENT, eventTiming')
% temp = eventTiming;
% temp(:,1:2) = temp(:,1:2).*samplePeriod;
% disp(temp)
% disp('-------')
%separate start/end events and adjust indices
if( count > 0 )
    %trim off events that happened before beginning and end of sweep
    %event started before sweep
    if( (eventTiming(1,1) == 1) )
        if(cutoffEventPresent == 0) %remove event, no previous info saved
            eventTiming = eventTiming(2:end,:);
            count = count - 1;
        else                        %add event start info from previous sweep
            eventTiming(1,1) = continuedEventStartTime;
            eventTiming(1,3) = continuedEventSweep;
        end
    elseif( rampingData == 1  && cutoffEventPresent == 1)
        disp('Filling in gap from prev event')
        endRampEvent = find(sweepData(1:eventTiming(1,1),2) < -20, 1);
        if( endRampEvent ~= 0 )
            eventTiming = [[continuedEventStartTime endRampEvent continuedEventSweep currentSweep 0 5 0 0];eventTiming];
            count = count + 1;
        end
    end
    if( primerData == 1 && cutoffEventPresent == 1)
        disp('Passing the previous primer info from the last sweep')
        eventVoltageTimingPrimer(1,:) = passedEventVoltageTimingPrimer
    end
    %event continues into next sweep
    if( ~isempty(eventTiming) ) %make sure didn't remove only event
        if( (eventTiming(end,2) == length(sweepData(:,1))) )
            disp('passing info again?')
            %return cutoff event for processing in next sweep
            cutoffEventSweep = currentSweep;
            cutoffEventStartTime = eventTiming(end,1);

            %remove event from eventTiming
            eventTiming = eventTiming(1:end-1,:);
            count = count - 1;

        end
    end
end

% disp('AFTER ADJUSTING INDICIES, eventTiming')
% temp = eventTiming;
% temp(:,1:2) = temp(:,1:2).*samplePeriod;
% disp(temp)
% disp(['cutoffEventFound: ' num2str(cutoffEventFound)])
% disp('-------')
%keep track of what events have voltage changes
if( ~isempty(eventTiming) )
    hasVoltageChange = zeros(length(eventTiming(:,1)));
end

%limit the distance back the voltage change can be
%only want to skip back to beginning of capacitive transient
vChangeStartLimit = 5000;

%initialize vChangeLoggingColumn (capacitive transient length)
% and event type
if( count ~= 0 )
    disp('hey blah')
    disp(['count: ' num2str(count) ' size(eventTiming)' num2str(size(eventTiming))])
    eventTiming = [eventTiming, zeros(count,3) ];
end

%if still events left do voltage event end check
%for cutoff event
if( (count > 0 || cutoffEventFound) && numSignals == 2 )

    %find where voltage changes and end events when it does
%     voltageChangeIndices = filter( [1,-1],1, sweepData(:,2) );
    vChangeWindow = 50;
    voltageChangeIndices = filter( [ones(1, vChangeWindow) -1*ones(1, vChangeWindow)],1, sweepData(:,2) );
    
%     figure(15); plot(0:samplePeriod:samplePeriod*(length(sweepData)-1), voltageChangeIndices, '-b.');

    %Finish up voltage cutoffs for events that spanned sweeps
    if( cutoffEventFound && (fishingData || rampingData || primerData))
        
        vChangeStart = find( abs( voltageChangeIndices(cutoffEventStartTime:-1:2) ) >= 5, 1 );
        if( vChangeStart ~= 0 )
            %get where event first goes above thresholdFishing
            if ( rampingData )
                fishingStart = find( sweepData (cutoffEventStartTime:end, 1 ) > -10, 1);
            else
                fishingStart = find( sweepData( cutoffEventStartTime:end, 1 ) > thresholdFishing, 1 );
            end
            
            if( ~isempty( fishingStart )  && (cutoffEventStartTime + fishingStart < length( sweepData(:,1) ) ) )
                disp('fixing cutoffEventStartTime')
                cutoffEventStartTime = cutoffEventStartTime + fishingStart
            end

            if( vChangeStart < cutoffEventStartTime  && vChangeStart <= vChangeStartLimit )
                %keep track of vChange start for amplitude calc purposes
                cutoffVChangeStart = vChangeStart;
            end
        end
    end
end

%keep track of what events have voltage changes
if( ~isempty(eventTiming) )
    hasVoltageChange = zeros(length(eventTiming(:,1)));
end
try
disp('eventTiming after changes, but before bulk of detection for primer')
temp = eventTiming;
temp(:,1:2) = temp(:,1:2).*samplePeriod;
disp(temp)
disp('------------')
catch
end
%for other events (non-cutoff)
if( count > 0 && numSignals == 2 )

    vChange = 0;
    vChangeStart = 0;
    %loop through found events and see if voltage changes during event
    for i = 1:length( eventTiming(:,1) )
        %search for voltage change from beginning of sweep if cutoff event
        if( eventTiming(i,3) ~= eventTiming(i,4) )
            vChange = find( abs( voltageChangeIndices(2:eventTiming(i,2)+30) ) >= voltageChangeEpsilon, 1 );

            %update new end
            if( vChange ~= 0 && ~primerData )
                eventTiming(i,2) = vChange; %don't need +1 because we start at 2 above
                hasVoltageChange(i) = 1;
                vChange = 0;
                
                %if voltage change at end of event;
                eventTiming(i,6) = 2;

            end

            if( fishingData == 1)
                %if cutoff event, bring in passed vChangeStart
                vChangeStart = passedVChangeStart;

                if( vChangeStart ~= 0 )
                    if( vChangeStart < eventTiming(i,1)  && vChangeStart <= vChangeStartLimit )
                        %keep track of vChange start for amplitude calc purposes
                        eventTiming(i,7) = vChangeStart;
                        eventTiming(i,1) = eventTiming(i,1) - vChangeStart;
                        
                        %if event type is 'initial fishing event',
                        %eventTiming(i,6) == 2
                        %make a 'regular fishing event' since there is a
                        %voltage change
                        %at the beginning (end fishing event wouldn't have
                        %that)
                        if( eventTiming(i,6) == 2 )
                            eventTiming(i,6) = 3;
                        elseif( eventTiming(i,6) == 0 )
                            eventTiming(i,6) = 4;
                        end                        
                    end
                end
            elseif( rampingData == 1)
                vChangeSecond = find(sweepData(eventTiming(i,2):end,2) < -20,1)
                if( vChangeSecond ~= 0 )
                    fprintf('Looking for change to neg\n length(sweepData(eventTiming(i,2):end)): %d\t eventTiming(i,2): %d\n', length(sweepData(eventTiming(i,2):end)), eventTiming(i,2));
                    eventTiming(i,2) = eventTiming(i,2) + vChangeSecond - 1;
                    hasVoltageChange(i) = 1;
                    eventTiming(i,6) = 5;
                else
                    disp('oops')
                end
            elseif( primerData == 1)
                %First Event @ start of sweep is a cutoff, see if the rest
                %of the volage changes are in this part, we'll check the
                %pass through to make sure that it didn't find the last of
                %the things before
                disp(['Event #' num2str(i) ' spans sweeps?'])
                vChangeHold = find(sweepData(1:eventTiming(i,2),2) > holdingVoltage - 2 & sweepData(1:eventTiming(i,2),2) < holdingVoltage + 2,1)
                vChangeProbe = find(sweepData(1:eventTiming(i,2),2) > probingVoltage - 2 & sweepData(1:eventTiming(i,2),2) < probingVoltage + 2,1)
                vChangeEject = find(sweepData(1:eventTiming(i,2),2) > ejectVoltage - 2 & sweepData(1:eventTiming(i,2),2) < ejectVoltage + 2,1)
                eventStartLength = length(sweepData(:,1)) - eventTiming(i,1);
                
                switch nnz(isnan(eventVoltageTimingPrimer(i,:)))
                    case 0
                        disp('case 0')
                        %No nans, do nothing except make the passed vector
                        %the correct one
                        %Which should of happened earlier
                        %eventVoltageTimingPrimer = passedEventVoltageTimingPrimer;
                    case 1
                        disp('case 1')
                        %Should only be the last one
                        if( isnan(eventVoltageTimingPrimer(i,end)) )
                            %see if there is an eject
                            if( ~isempty(vChangeEject) )
                                %Add it to the end
                                eventVoltageTimingPrimer(i,:) = [eventVoltageTimingPrimer(i,1:2), vChangeEject+eventStartLength];
                                eventTiming(i,6) = 11;
                            else
                                %if there isn't it must not be there
                                eventVoltageTimingPrimer(i,:) = [eventVoltageTimingPrimer(i,1:2), 0];
                                eventTiming(i,6) = 10;
                            end
                        else
                            disp('PRIMER ERROR: only 1 nan and it is not the last one')
                            eventVoltageTimingPrimer(i,:) = [0,0,0];
                        end
                    case 2
                        disp('case 2')
                        %Should be the last 2
                        if( nnz(isnan(eventVoltageTimingPrimer(i,end-1:end))) == 2)
                            %Check for probing && Eject
                            if( ~isempty(vChangeProbe) && ~isempty(vChangeEject) )
                                %Found the probe and eject
                                eventVoltageTimingPrimer(i,:) = [eventVoltageTimingPrimer(i,1), vChangeProbe + eventStartLength, vChangeEject + eventStartLength];
                                eventTiming(i,6) = 11;
                            elseif( ~isempty(vChangeProbe) && isempty(vChangeEject) )
                                %only the probe, so no eject
                                eventVoltageTimingPrimer(i,:) = [eventVoltageTimingPrimer(i,1), vChangeProbe + eventStartLength, 0];
                                eventTiming(i,6) = 10;
                            else
                                %something has gone terribly wrong
                                disp('PRIMER ERROR: 2 nans and the probe is not found')
                                eventVoltageTimingPrimer(i,:) = [0,0,0];
                            end
                        else
                            disp('PRIMER ERROR: 2 nans and it is not the last two')
                            eventVoltageTimingPrimer(i,:) = [0,0,0];
                        end
                    case 3
                        disp('case 3')
                        %Totally empty, it might be a non primer event or
                        %it might just have start it the previous sweep
                        if( ~isempty(vChangeHold) && ~isempty(vChangeProbe) && isempty(vChangeEject) )
                            eventTiming(i,6) = 10;
                            eventVoltageTimingPrimer(i,:) = [vChangeHold + eventStartLength, vChangeProbe + eventStartLength, 0];
                            fprintf('vChangeHold/Probe = %.3f %.3f [%.3f,%.3f]\n', vChangeHold, vChangeProbe, eventTiming(i,1).*samplePeriod, eventTiming(i,2).*samplePeriod);
                        elseif( ~isempty(vChangeHold) && ~isempty(vChangeProbe) && ~isempty(vChangeEject) )
                            eventTiming(i,6) = 11;
                            eventVoltageTimingPrimer(i,:) = [vChangeHold + eventStartLength, vChangeProbe + eventStartLength, vChangeEject + eventStartLength];
                            fprintf('vChangeHold/Probe = %.3f %.3f %.3f [%.3f,%.3f]\n', vChangeHold, vChangeProbe, vChangeEject, eventTiming(i,1).*samplePeriod, eventTiming(i,2).*samplePeriod);
                        elseif( isempty(vChangeHold) && isempty(vChangeProbe) && isempty(vChangeEject) )
                            fprintf('vChangeHold/Probe = nan [%.3f, %.3f]\n', eventTiming(i,1).*samplePeriod, eventTiming(i,2).*samplePeriod);
                        else
                            disp('PRIMER ERROR: passed event timing was empty but voltage don''t match up')
                            eventVoltageTimingPrimer(i,:) = [0,0,0];
                        end
       
                end
                    %Couldn't find at least one of the voltages
            end

            %if same sweep events
        else

            %If we have captured and reduced for the ramp
            if ( rampingData && sweepData(eventTiming(i,1)+30,2) < 100 )
                %for rampingData we are looking for the change to negative
                vChangeSecond = find(sweepData(eventTiming(i,2):end,2) < -20,1)
                if( vChangeSecond ~= 0 )
                    fprintf('Looking for change to neg\n length(sweepData(eventTiming(i,2):end)): %d\t eventTiming(i,2): %d\n', length(sweepData(eventTiming(i,2):end)), eventTiming(i,2));
                    eventTiming(i,2) = eventTiming(i,2) + vChangeSecond - 1;
                    hasVoltageChange(i) = 1;
                    eventTiming(i,6) = 5;
                else
                    %if the eject is in the next sweep
                    %return cutoff event for processing in next sweep
                    cutoffEventFoundRamp = [1,i];
                    cutoffEventSweep = currentSweep
                    cutoffEventStartTime = eventTiming(end,1)
                end
%                 figure(16);
%                 plot(eventTiming(i,1)*samplePeriod:samplePeriod:(length(vChangeSecond)+eventTiming(i,1)-1)*samplePeriod,vChangeSecond)
%                 figure(12);
            elseif ( primerData )
                
                vChangeHold = find(sweepData(eventTiming(i,1):eventTiming(i,2),2) > holdingVoltage - 2 & sweepData(eventTiming(i,1):eventTiming(i,2),2) < holdingVoltage + 2,1);
                vChangeProbe = find(sweepData(eventTiming(i,1):eventTiming(i,2),2) > probingVoltage - 2 & sweepData(eventTiming(i,1):eventTiming(i,2),2) < probingVoltage + 2,1);
                vChangeEject = find(sweepData(eventTiming(i,1):eventTiming(i,2),2) > ejectVoltage - 2 & sweepData(eventTiming(i,1):eventTiming(i,2),2) < ejectVoltage + 2,1);
                %We have an event that has both the probing and holding
                %voltage present, this is an event we want
                if( ~isempty(vChangeHold) && ~isempty(vChangeProbe) && isempty(vChangeEject) )
                    eventTiming(i,6) = 10;
                    eventVoltageTimingPrimer(i,:) = [vChangeHold, vChangeProbe, 0];
                    fprintf('vChangeHold/Probe = %.3f %.3f [%.3f,%.3f]\n', vChangeHold, vChangeProbe, eventTiming(i,1).*samplePeriod, eventTiming(i,2).*samplePeriod);
                elseif( ~isempty(vChangeHold) && ~isempty(vChangeProbe) && ~isempty(vChangeEject) )
                    eventTiming(i,6) = 11;
                    eventVoltageTimingPrimer(i,:) = [vChangeHold, vChangeProbe, vChangeEject];
                    fprintf('vChangeHold/Probe = %.3f %.3f %.3f [%.3f,%.3f]\n', vChangeHold, vChangeProbe, vChangeEject, eventTiming(i,1).*samplePeriod, eventTiming(i,2).*samplePeriod);
                else
                    fprintf('vChangeHold/Probe = nan [%.3f, %.3f]\n', eventTiming(i,1).*samplePeriod, eventTiming(i,2).*samplePeriod);
                end
                
            else
            
                vChange = find( abs( voltageChangeIndices(eventTiming(i,1):eventTiming(i,2)+30) ) >= voltageChangeEpsilon, 1 );

                if( vChange ~= 0)
                    fprintf('vChange = %.3f [%.3f,%.3f]\n', vChange, eventTiming(i,1).*samplePeriod, eventTiming(i,2).*samplePeriod);
                else
                    fprintf('vChange = nan [%.3f, %.3f]\n', eventTiming(i,1).*samplePeriod, eventTiming(i,2).*samplePeriod);
                end
                
                %update new end
                if( vChange ~= 0 )
                    eventTiming(i,2) = eventTiming(i,1) + vChange - 2;
                    %fprintf('new event time: [%.3f,%.3f]\n', eventTiming(i,1).*samplePeriod, eventTiming(i,2).*samplePeriod);

                    hasVoltageChange(i) = 1;
                    vChange = 0;

                    %if voltage change at end of event;
                    eventTiming(i,6) = 2;


                end
            end
            %search for voltage change at beginning of event if fishing data
            if( fishingData == 1 || rampingData == 1)
                vChangeStart = find( abs( voltageChangeIndices(eventTiming(i,1):-1:2) ) >= 50, 1 );
                
                if( vChangeStart ~= 0 )
                    %get where event first goes above thresholdFishing
                    if ( rampingData )
                        fishingStart = find( sweepData( eventTiming(i,1):eventTiming(i,2), 1 ) > -10, 1 );
                    else
                        fishingStart = find( sweepData( eventTiming(i,1):eventTiming(i,2), 1 ) > thresholdFishing, 1 );
                    end

                    if( ~isempty( fishingStart ) && ( (eventTiming(i,1) + fishingStart) < eventTiming(i,2)) )
                        eventTiming(i,1) = eventTiming(i,1) + fishingStart; %move "start" of event back to above threshold
                    end
                    
                    if( vChangeStart < eventTiming(i,1)  && vChangeStart <= vChangeStartLimit )
                        %keep track of vChange start for amplitude calc purposes
                        eventTiming(i,7) = vChangeStart;
                        eventTiming(i,1) = eventTiming(i,1) - vChangeStart;
                        
                        %if event type is 'initial fishing event',
                        %eventTiming(i,6) == 2
                        %make a 'regular fishing event' since there is a
                        %voltage change
                        %at the beginning (end fishing event wouldn't have
                        %that)
                        if( eventTiming(i,6) == 2 )
                            eventTiming(i,6) = 3;
                        elseif( eventTiming(i,6) == 0 )
                            eventTiming(i,6) = 4;
                        end

                    end
                end
            end
        end

    end
end

if( cutoffEventFoundRamp(1) == 1)
    cutoffEventSweep = currentSweep;
    cutoffEventStartTime = eventTiming(cutoffEventFoundRamp(2),1);

    %remove event from eventTiming
    eventTiming(cutoffEventFoundRamp(2):end,:) = [];
    count = count - 1;

end

%PUT CHECK TO SEE IF THE START OF EVENTS ARE IN THE RANGE OF PERVIOUS EVENTS
%WERE PUSHING FORWARD THE LIMIT QUITE A BIT

if( rampingData == 1)
    del_inds = [];
    for i=2:length(eventTiming(:,1))
        if( eventTiming(i,3) == eventTiming(i,4)  && hasVoltageChange(i,1) == 0)
            eventToCompare = find(hasVoltageChange(1:i-1,1) == 1,1,'last');
            if( ~isempty(eventToCompare) )
                if( eventTiming(i,1) > eventTiming(eventToCompare,1) && eventTiming(i,2) < eventTiming(eventToCompare,2) )
                    del_inds = [del_inds;i];
                end
            end
        end
    end

    eventTiming(del_inds,:) = [];
    count = count - length(del_inds);
end

if( count > 0 )
    %find where event transitions start/end and chop off transitions
    %event needs to be over for 0.2ms for end to be detected
    maxTransitionLength = round( (0.1/1000) / samplePeriod );


    for i = 1:length(eventTiming(:,1))

        %don't remove start transition for cutoff events
        %that is done in call to parseEvents when event start
        if( eventTiming(i,3) == eventTiming(i,4) && ~rampingData && ~primerData )

%             if( eventTiming(i,6) ~= 3 && eventTiming(i,6) ~= 4) %want event to start with voltage change
                %for each event, loop back until a sharp change in current is spotted
                for j = 0:maxTransitionLength
                    currentDelta = abs( sweepData(eventTiming(i,1)+j,1) - sweepData(eventTiming(i,1)+maxTransitionLength, 1) );
                    if( currentDelta <= stdDev )
                        eventTiming(i,1) = eventTiming(i,1)+j;
                        break
                    end
                end
%             end
        end
        %for each event, loop back until a sharp change in current is spotted
        %unless voltage change has already been used to mark end of event
%         hasVoltage = 1;
%         if( numSignals < 2 )
%             hasVoltage = 0;
%         elseif( isequal( hasVoltageChange(i), 1 ) )
%             hasVoltage = 1;
%         end
% 
%         if( hasVoltage )
%             for j = 0:maxTransitionLength
%                 currentDelta = abs( sweepData(eventTiming(i,2)-j,1) - sweepData(eventTiming(i,2)-maxTransitionLength, 1) );
%                 if( currentDelta <= stdDev )
%                     eventTiming(i,2) = eventTiming(i,2)-j;
%                     break
%                 end
%             end
%         end
% 

        %tag events to a particular DNA molecule if fishing using fishing start
        %events
        if( eventTiming(i,6) == 2 )
            fprintf('\nnew DNA %d in sweep %d at %fms (%d)', currentDNAnumber, currentSweep, eventTiming(i,1).*samplePeriod, eventTiming(i,1) );            
            currentDNAnumber = currentDNAnumber + 1;
            eventTiming(i,5) = currentDNAnumber;            
        elseif( eventTiming(i,6) == 4 )
            fprintf(' and ended in sweep %d at %fms (%d)', currentSweep, eventTiming(i,1).*samplePeriod, eventTiming(i,1) );            
            eventTiming(i,5) = currentDNAnumber;            
        elseif( eventTiming(i,6) == 3 )
            eventTiming(i,5) = currentDNAnumber;
        else
            eventTiming(i,5) = 0;            
        end
    
    end

end


if( count > 0 )

    %find events that are still the minium length and update eventTiming
    %if events start and end in different sweeps, add a sweep length for
    %the number of sweeps the event lasts for
    indicesToKeep = find( ( ((eventTiming(:,4)-eventTiming(:,3)).*numSamples)+(eventTiming(:,2)-eventTiming(:,1)) ) >= minEventSamples );
    eventTiming = eventTiming( indicesToKeep, :);
    eventVoltageTimingPrimer = eventVoltageTimingPrimer(indicesToKeep,:);
    %find min current for each event
    %store min amplitude for each event
    if( ~isempty( eventTiming ) )
        for i = 1:length(eventTiming(:,1))
            %fprintf('eventnumber %d: sweepData[%d,%d] %d\n', i, eventTiming(i,1), eventTiming(i,2), length(sweepData(:,1) ) );
            if( eventTiming(i,2) > eventTiming(i,1) )
                eventTiming(i,8) = min( sweepData( eventTiming(i,1):eventTiming(i,2), 1 ) );
            end
        end
    end
end


try
disp('eventVoltageTimingPrimer: ')
disp(eventVoltageTimingPrimer)
disp('eventTiming(:,6) == ')
disp(eventTiming(:,6))
disp('cutoffEventFoundPrimer: ')
disp(cutoffEventFoundPrimer)
catch
end

disp('-------------- END SWEEP ----------------')
disp('-------------- END SWEEP ----------------')
disp('-------------- END SWEEP ----------------')



