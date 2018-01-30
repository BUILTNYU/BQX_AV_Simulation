function [ output ] = BQX_Project1FINAL_v1 ( fleetSize )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EZ10 Simulation for the rush hour
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all

test = 0;
contSimulation = 0;

% Inputs
lambda = xlsread('StreetcarData', 'Data1'); %arrival rates (passenger per hour)
stopsPosition = xlsread('StreetcarData', 'Data2');  %miles
numberOfStops = size(stopsPosition,2);
shuttleAverageSpeed = 10.6; %miles/hour
stopTime = 15/3600;  %hour
shuttleCapacity = 12;
simulationTime = 2; %hour
numberShuttles = fleetSize; %total number of shuttles

if mod(numberOfStops,2) == 0
    for i=1:numberOfStops/2
        stopsPosition(numberOfStops + i) = (stopsPosition(2*i) + stopsPosition(2*i-1))/2;
    end
else
    for i=1:(numberOfStops-1)/2
        stopsPosition(numberOfStops + i) = (stopsPosition(2*i) + stopsPosition(2*i-1))/2;
    end
    stopsPosition((3*numberOfStops + 1)/2) = (stopsPosition(numberOfStops) + stopsPosition(numberOfStops-1))/2;
end

% Passengers Arrival Time (Poisson Distribution)

arrivalTime = zeros(numberOfStops, numberOfStops); %arrival time of each passenger in each station
auxTime = zeros(numberOfStops, numberOfStops);
numberPassenger = zeros(numberOfStops, numberOfStops); %total number of passengers in each station

for i=1:numberOfStops
    for j=1:numberOfStops
        while auxTime(i,j) < simulationTime
            a = rand();
            b = -log(a)/lambda(i,j);
            auxTime(i,j) = auxTime(i,j) + b; 
            if auxTime(i,j) < simulationTime
                numberPassenger(i,j) = numberPassenger(i,j)+1;                
                contSimulation = contSimulation + 1;
                arrivalTime(i,j,numberPassenger(i,j)) = auxTime(i,j);
            end
        end
    end
end

% Variables

error1 = zeros(numberShuttles,numberOfStops,1+shuttleCapacity);
error2 = zeros(numberShuttles,numberOfStops,1+shuttleCapacity);

passengerServed = zeros(numberOfStops, numberOfStops); %total number of passengers served per station
shuttlePosition = zeros(1,4);
shuttleIdleTime = zeros(1, numberShuttles); %time event 3 will occur
shuttleArrivalTime = zeros(1, numberShuttles); %time event 2 will occur
shuttleLeavingTime = zeros(1, numberShuttles); %time event 4 will occur
shuttleDest1 = zeros(numberShuttles,size(stopsPosition,2)); %binary variable with every destination of each shuttle moving fowards
shuttleDest2 = zeros(numberShuttles,size(stopsPosition,2)); %binary variable with every destination of each shuttle moving backwards
shuttleDestPass1 = zeros(numberShuttles,numberOfStops,1+shuttleCapacity); %number of passenger going to each destination for each shuttle moving fowards
shuttleDestPass2 = zeros(numberShuttles,numberOfStops,1+shuttleCapacity); %number of passenger going to each destination for each shuttle moving backwards
auxTrack1 = zeros(numberShuttles,numberOfStops,1+shuttleCapacity); %auxiliary variable to help track the passengers
auxTrack2 = zeros(numberShuttles,numberOfStops,1+shuttleCapacity); %auxiliary variable to help track the passengers
shuttlePickPass1 = zeros(numberShuttles,numberOfStops,1+shuttleCapacity); %number of passenger that is going to be picked at every destination for each shuttle moving fowards
shuttlePickPass2 = zeros(numberShuttles,numberOfStops,1+shuttleCapacity); %number of passenger that is going to be picked at every destination for each shuttle moving backwards
shuttleNumberPass = int8(zeros(1, numberShuttles)); %[shuttle totalNumberPass]
travelTime = zeros(1,5);
waitingTime = zeros(1,4);
totalTime = zeros(1,5);
for i=1:numberShuttles
    if (mod(i, numberOfStops) == 0)
        shuttlePosition(i,1) = stopsPosition(numberOfStops); %position o the shuttle at a specified moment 
    else
        shuttlePosition(i,1) = stopsPosition(mod(i, numberOfStops)); 
    end
    shuttlePosition(i,2) = 0; %time in which shuttlePosition was accounted for
    shuttlePosition(i,3) = 0; %destination (next stop of the shuttle)
    shuttlePosition(i,4) = 0; %direction ( -1 means shuttle moving backwards; 0 means that the shuttle isn't moving; 1 means shuttle moving forwards)
    shuttleIdleTime(1,i) = 0; %time event 3 will occur
    shuttleArrivalTime(1,i) = 0; %time event 2 will occur
    shuttleLeavingTime(1,i) = 0; %time event 4 will occur
    shuttleDest1(i,:) = 0; 
    shuttleDest2(i,:) = 0;
    shuttleDestPass1(i,:,1) = 0; 
    shuttleDestPass2(i,:,1) = 0;
    shuttlePickPass1(i,:,1) = 0;
    shuttlePickPass2(i,:,1) = 0;    
    auxTrack1(i,:,1) = 0; 
    auxTrack2(i,:,1) = 0;
    shuttleNumberPass(1,i) = 0; %[shuttle totalNumberPass]    
    error1(i,:,1) = 0; 
    error2(i,:,1) = 0; 
end

time = 0;

% Events Simulation
while time < simulationTime
    %Choose event
    aux = simulationTime;
    event = zeros(1,4);
    
    %Passenger Arrival
    for i=1:numberOfStops
        for j=1:numberOfStops
            if (passengerServed(i,j)+1 <= numberPassenger(i,j))
                if (arrivalTime(i,j,passengerServed(i,j)+1) < aux && arrivalTime(i,j,passengerServed(i,j)+1) > 0)
                    aux = arrivalTime(i,j,passengerServed(i,j)+1);
                    event = [1 i j passengerServed(i,j)+1];
                end
            end
        end
    end
    
    for i=1:numberShuttles
        %Passenger pick/drop
        if (shuttleArrivalTime(1,i) < aux && shuttleArrivalTime(1,i) > 0)
            aux = shuttleArrivalTime(1,i);
            event = [2 i 0 0];
        end
        %Shuttle IDLE
        if (shuttleIdleTime(1,i) < aux && shuttleIdleTime(1,i) > 0)
            aux = shuttleIdleTime(1,i);
            event = [3 i 0 0];
        end
        %Shuttle leaving station
        if (shuttleLeavingTime(1,i) < aux && shuttleLeavingTime(1,i) > 0)
            aux = shuttleLeavingTime(1,i);
            event = [4 i 0 0];
        end
    end
        
    time = aux;
    
    % Event 1 (Passenger Arrival)
    if (event(1,1) == 1)
        %Shuttle expected arrival time to the passenger
        shuttleTime = zeros(1,numberShuttles);
        for i=1:numberShuttles
            pos = shuttlePosition(i,1) + shuttlePosition(i,4)*shuttleAverageSpeed*(time - shuttlePosition(i,2));
            shuttleTime(1,i) = 100;
            %Shuttle is full
            if (shuttleNumberPass(1,i) >= shuttleCapacity)
                shuttleTime(1,i) = 100; %random large number
                auxbin = 0;
                %Determine direction
                if (stopsPosition(shuttlePosition(i,3)) > pos)
                    auxbin = 1;
                end
                if (event(1,2) == shuttlePosition(i,3) && shuttleCapacity - shuttleNumberPass(1,i) + shuttleDestPass1(i,shuttlePosition(i,3),1)*auxbin + shuttleDestPass2(i,shuttlePosition(i,3),1)*(1-auxbin) > 0)
                    if (auxbin == 1 && event(1,2) < event(1,3))
                        if (shuttlePosition(i,4) ~= 0)
                            shuttleTime(1,i) = shuttleArrivalTime(1,i) - time;
                        else
                            shuttleTime(1,i) = shuttleLeavingTime(1,i) + abs(stopsPosition(event(1,2)) - pos)/shuttleAverageSpeed - time;
                        end
                    end
                    if (auxbin == 0 && event(1,2) > event(1,3))
                        if (shuttlePosition(i,4) ~= 0)
                            shuttleTime(1,i) = shuttleArrivalTime(1,i) - time;
                        else
                            shuttleTime(1,i) = shuttleLeavingTime(1,i) + abs(stopsPosition(event(1,2)) - pos)/shuttleAverageSpeed - time;
                        end
                    end
                end
            %Passenger wants to move forwards    
            elseif (event(1,2) < event(1,3) && shuttleNumberPass(1,i) ~= 0)
                %Shuttle moving forwards
                if (shuttlePosition(i,4) == 1)
                    %Shuttle moving forwards and shuttle_pos < passenger_pos
                    if (stopsPosition(event(1,2)) > pos)
                        %Determine between which stations the shuttle is located
                        posaux = numberOfStops;
                        for j=2:numberOfStops
                            if (stopsPosition(j-1) <= pos && stopsPosition(j) > pos)
                                posaux = j-1;
                            end
                        end
                        timeaux = 0;
                        for j=posaux+1:event(1,2)
                            if (shuttleDest1(i,j) == 1 && j ~= event(1,2))
                                timeaux = timeaux + stopTime;
                            end
                        end
                        shuttleTime(1,i) = abs(stopsPosition(event(1,2)) - pos)/shuttleAverageSpeed + timeaux;
                    %Shuttle moving forwards and shuttle_pos > passenger_pos
                    else
                        shuttleTime(1,i) = 100; %random large number
                    end
                end
                %Shuttle moving backwards
                if (shuttlePosition(i,4) == -1)
                    %Determine last destination backwards of the shuttle
                    posaux2 = 0;
                    auxbin = 0;
                    while (auxbin == 0 && posaux2 < numberOfStops)
                        posaux2 = posaux2 + 1;
                        if (shuttleDest2(i,posaux2) == 1)
                            auxbin =1;
                        end
                    end
                    %Determine between which stations the shuttle is located
                    posaux = numberOfStops;
                    for j=2:numberOfStops
                        if (stopsPosition(j-1) <= pos && stopsPosition(j) > pos)
                            posaux = j-1;
                        end
                    end
                    timeaux = 0;
                    for j=posaux2:posaux
                        if (shuttleDest2(i,j) == 1)
                            if (j == posaux2)
                                if (posaux2 ~= event(1,2))
                                    timeaux = timeaux + stopTime;
                                end
                            else
                                timeaux = timeaux + stopTime;
                            end
                        end
                    end
                    %Passenger_pos >= last_shuttle_backwards_dest
                    if (posaux2 <= event(1,2) && posaux2 == 1)
                        for j=posaux2:event(1,2)
                            if (shuttleDest1(i,j) == 1 && j ~= posaux2 && j ~= event(1,2))
                                timeaux = timeaux + stopTime;
                            end
                        end
                        shuttleTime(1,i) = abs(stopsPosition(event(1,2)) + pos - 2*stopsPosition(posaux2))/shuttleAverageSpeed + timeaux; 
                    %Passenger_pos < last_shuttle_backwards_dest
                    else
                        shuttleTime(1,i) = 100; 
                    end
                end
                %Shuttle not moving
                if (shuttlePosition(i,4) == 0)
                    %Shuttle not moving and heading forwards with shuttle_pos <= passenger_pos 
                    if(pos <= stopsPosition(event(1,2)) && pos < stopsPosition(shuttlePosition(i,3)))
                        %Determine between which stations the shuttle is located
                        posaux = numberOfStops;
                        for j=2:numberOfStops
                            if (stopsPosition(j-1) <= pos && stopsPosition(j) > pos)
                                posaux = j-1;
                            end
                        end
                        timeaux = 0;
                        for j=posaux:event(1,2)
                            if (shuttleDest1(i,j) == 1 && j ~= event(1,2))
                                timeaux = timeaux + stopTime;
                            end
                        end
                        if (posaux ~= event(1,2))
                            timeaux = timeaux + shuttleLeavingTime(1,i) - time;
                        end
                        shuttleTime(1,i) = abs(stopsPosition(event(1,2)) - pos)/shuttleAverageSpeed + timeaux;
                    %Shuttle not moving and heading forwards with shuttle_pos > passenger_pos 
                    elseif (pos > stopsPosition(event(1,2)) && pos < stopsPosition(shuttlePosition(i,3)))
                        shuttleTime(1,i) = 100; %random large number
                    %Shuttle not moving and heading backwards 
                    elseif (pos > stopsPosition(shuttlePosition(i,3)))     
                        %Determine last destination backwards of the shuttle
                        posaux2 = 0;
                        auxbin = 0;
                        while (auxbin == 0 && posaux2 < numberOfStops)
                            posaux2 = posaux2 + 1;
                            if (shuttleDest2(i,posaux2) == 1)
                                auxbin = 1;
                            end
                        end
                        %Determine in which station the shuttle is located
                        for j=1:numberOfStops
                            if (stopsPosition(j) == pos)
                                posaux = j;
                            end
                        end
                        timeaux = 0;
                        for j=posaux2:posaux
                            if (shuttleDest2(i,j) == 1)
                                if (j == posaux2)
                                    if (posaux2 ~= event(1,2))
                                        timeaux = timeaux + stopTime;
                                    end
                                else
                                    timeaux = timeaux + stopTime;
                                end
                            end
                        end
                        if (posaux ~= event(1,2))
                            timeaux = timeaux + shuttleLeavingTime(1,i) - time;
                        end
                        %Passenger_pos >= last_shuttle_backwards_dest
                        if (posaux2 <= event(1,2) && posaux2 == 1)
                            for j=posaux2:event(1,2)
                                if (shuttleDest1(i,j) == 1 && j ~= posaux2 && j ~= event(1,2))
                                    timeaux = timeaux + stopTime;
                                end
                            end
                            shuttleTime(1,i) = abs(stopsPosition(event(1,2)) + pos - 2*stopsPosition(posaux2))/shuttleAverageSpeed + timeaux; 
                        %Passenger_pos < last_shuttle_backwards_dest
                        else
                            shuttleTime(1,i) = 100; 
                        end
                    end
                end
            %Passenger wants to move backwards
            elseif (event(1,2) > event(1,3) && shuttleNumberPass(1,i) ~= 0)
                %Shuttle moving forwards
                if (shuttlePosition(i,4) == 1)                    
                    %Determine last destination forwards of the shuttle
                    for j=1:numberOfStops
                        if (shuttleDest1(i,j) == 1)
                            posaux2 = j;
                        end
                    end
                    %Determine between which stations the shuttle is located
                    posaux = numberOfStops;
                    for j=2:numberOfStops
                        if (stopsPosition(j-1) <= pos && stopsPosition(j) > pos)
                            posaux = j-1;
                        end
                    end
                    timeaux = 0;
                    for j=posaux+1:posaux2
                        if (shuttleDest1(i,j) == 1)
                            timeaux = timeaux + stopTime;
                        end
                    end
                    %Passenger_pos <= last_shuttle_forwards_dest
                    if (posaux2 >= event(1,2) && posaux2 == numberOfStops)
                        for j=event(1,2):posaux-1
                            if (shuttleDest2(i,j) == 1 && j ~= posaux2 && j ~= event(1,2))
                                timeaux = timeaux + stopTime;
                            end
                        end
                        shuttleTime(1,i) = abs(stopsPosition(event(1,2)) + pos - 2*stopsPosition(posaux2))/shuttleAverageSpeed + timeaux;                        
                    %Passenger_pos > last_shuttle_forwards_dest
                    else
                        shuttleTime(1,i) = 100;
                    end      
                end
                %Shuttle moving backwards
                if (shuttlePosition(i,4) == -1) 
                    %Shuttle moving backwards and shuttle_pos >= passenger_pos
                    if (stopsPosition(event(1,2)) <= pos)
                        %Determine between which stations the shuttle is located
                        posaux = numberOfStops;
                        for j=2:numberOfStops
                            if (stopsPosition(j-1) <= pos && stopsPosition(j) > pos)
                                posaux = j-1;
                            end
                        end
                        timeaux = 0;
                        for j=event(1,2):posaux
                            if (shuttleDest2(i,j) == 1 && j ~= event(1,2))
                                timeaux = timeaux + stopTime;
                            end
                        end   
                        shuttleTime(1,i) = abs(stopsPosition(event(1,2)) - pos)/shuttleAverageSpeed + timeaux;
                    %Shuttle moving backwards and shuttle_pos < passenger_pos
                    else
                        shuttleTime(1,i) = 100; %random large number
                    end
                end
                %Shuttle not moving
                if (shuttlePosition(i,4) == 0)
                    %Shuttle not moving and heading backwards with shuttle_pos >= passenger_pos
                    if (pos >= stopsPosition(event(1,2)) && pos > stopsPosition(shuttlePosition(i,3)))
                        %Determine between which stations the shuttle is located
                        posaux = numberOfStops;
                        for j=2:numberOfStops
                            if (stopsPosition(j-1) <= pos && stopsPosition(j) > pos)
                                posaux = j-1;
                            end
                        end
                        timeaux = 0;
                        for j=event(1,2):posaux
                            if (shuttleDest2(i,j) == 1 && j ~= event(1,2))
                                timeaux = timeaux + stopTime;
                            end
                        end                     
                        if (posaux ~= event(1,2))
                            timeaux = timeaux + shuttleLeavingTime(1,i) - time;
                        end
                        shuttleTime(1,i) = abs(stopsPosition(event(1,2)) - pos)/shuttleAverageSpeed + timeaux;
                    %Shuttle not moving and heading backwards with shuttle_pos < passenger_pos
                    elseif (pos < stopsPosition(event(1,2)) && pos > stopsPosition(shuttlePosition(i,3)))
                            shuttleTime(1,i) = 100; %random large number
                    %Shuttle not moving and heading forwards
                    elseif (pos < stopsPosition(shuttlePosition(i,3)))                                      
                        %Determine last destination forwards of the shuttle
                        posaux2 = 0;
                        for j=1:numberOfStops
                            if (shuttleDest1(i,j) == 1)
                                posaux2 = j;
                            end
                        end
                        %Determine in which station the shuttle is located
                        for j=1:numberOfStops
                            if (stopsPosition(j) == pos)
                                posaux = j;
                            end
                        end
                        timeaux = 0;
                        for j=posaux:posaux2
                            if (shuttleDest1(i,j) == 1 && j ~= posaux)
                                timeaux = timeaux + stopTime;
                            end
                        end    
                        %Passenger_pos <= last_shuttle_forwards_dest
                        if (posaux2 >= event(1,2) && posaux2 == numberOfStops) 
                            for j=event(1,2):posaux2
                                if (shuttleDest2(i,j) == 1 && j ~= posaux2 && j ~= event(1,2))
                                    timeaux = timeaux + stopTime;
                                end
                            end
                            shuttleTime(1,i) = abs(stopsPosition(event(1,2)) + pos - 2*stopsPosition(posaux2))/shuttleAverageSpeed + timeaux;  
                        %Passenger_pos < last_shuttle_forwards_dest
                        else   
                            shuttleTime(1,i) = 100;
                        end
                    end
                end                
            %Shuttle that is serving no passenger
            elseif (shuttleNumberPass(1,i) == 0)
                shuttleTime(1,i) = abs(stopsPosition(event(1,2)) - pos)/shuttleAverageSpeed;
            end
        end
        %Choose the faster shuttle
        fasterShuttle = [0 1000]; %[shuttle time]
        for i=1:numberShuttles
            if (fasterShuttle(2) > shuttleTime(1,i))
                fasterShuttle(1) = i;
                fasterShuttle(2) = shuttleTime(1,i);
            end
        end
        %Update faster shuttle information
        passengerServed(event(1,2),event(1,3)) = passengerServed(event(1,2),event(1,3)) + 1;        
        pos = shuttlePosition(fasterShuttle(1),1) + shuttlePosition(fasterShuttle(1),4)*shuttleAverageSpeed*(time - shuttlePosition(fasterShuttle(1),2));   
        if (event(1,2) < event(1,3))
            shuttleDestPass1(fasterShuttle(1), event(1,3),1) = shuttleDestPass1(fasterShuttle(1), event(1,3),1) + 1; 
            shuttleDestPass1(fasterShuttle(1), event(1,3),1+shuttleDestPass1(fasterShuttle(1), event(1,3),1)) = time;
            auxTrack1(fasterShuttle(1), event(1,3),1+shuttleDestPass1(fasterShuttle(1), event(1,3),1)) = event(1,2);
            shuttleDest1(fasterShuttle(1), event(1,3)) = 1;
        else
            shuttleDestPass2(fasterShuttle(1), event(1,3),1) = shuttleDestPass2(fasterShuttle(1), event(1,3),1) + 1; 
            shuttleDestPass2(fasterShuttle(1), event(1,3),1+shuttleDestPass2(fasterShuttle(1), event(1,3),1)) = time;
            auxTrack2(fasterShuttle(1), event(1,3),1+shuttleDestPass2(fasterShuttle(1), event(1,3),1)) = event(1,2);
            shuttleDest2(fasterShuttle(1), event(1,3)) = 1;
        end
        %Shuttle Moving and serving passengers | No need for a new shuttle
        if (shuttlePosition(fasterShuttle(1),4) ~= 0 && shuttleNumberPass(1,fasterShuttle(1)) ~= 0)
            shuttlePosition(fasterShuttle(1),1) = pos; 
            shuttlePosition(fasterShuttle(1),2) = time;
            shuttleIdleTime(1,fasterShuttle(1)) = 0;
            shuttleLeavingTime(1,fasterShuttle(1)) = 0;
            shuttleNumberPass(1,fasterShuttle(1)) = shuttleNumberPass(1,fasterShuttle(1)) + 1; %[shuttle totalNumberPass] 
            %Passenger moving forwards
            if (event(1,2) < event(1,3))
                %Shuttle moving forwards
                if (shuttlePosition(fasterShuttle(1),4) == 1)
                    %Shuttle moving forwards and shuttle_pos < passenger_pos
                    if(pos < stopsPosition(event(1,2)))
                        shuttleDest1(fasterShuttle(1), event(1,2)) = 1;
                        shuttlePickPass1(fasterShuttle(1), event(1,2),1) = shuttlePickPass1(fasterShuttle(1), event(1,2),1) + 1;
                        shuttlePickPass1(fasterShuttle(1), event(1,2),1+shuttlePickPass1(fasterShuttle(1), event(1,2),1)) = time;
                        error1(fasterShuttle(1), event(1,2),1+shuttlePickPass1(fasterShuttle(1), event(1,2),1)) = 11;
                        shuttleDest1(fasterShuttle(1), event(1,3)) = 1;
                        if (abs(stopsPosition(event(1,2)) - pos) < abs(stopsPosition(shuttlePosition(fasterShuttle(1),3)) - pos))
                            shuttlePosition(fasterShuttle(1),3) = event(1,2); 
                        end
                    %Shuttle moving forwards and shuttle_pos >= passenger_pos (don't happen)
                    else
                        test = 1;
                    end
                end
                %Shuttle moving backwards
                if (shuttlePosition(fasterShuttle(1),4) == -1)
                    %Determine last destination backwards of the shuttle
                    posaux2 = 0;
                    auxbin = 0;
                    while (auxbin == 0 && posaux2 < numberOfStops)
                        posaux2 = posaux2 + 1;
                        if (shuttleDest2(fasterShuttle(1),posaux2) == 1)
                            auxbin =1;
                        end
                    end
                    %Passenger_pos >= last_shuttle_backwards_dest
                    if (posaux2 <= event(1,2) && posaux2 == 1)
                        if (posaux2 ~= event(1,2))
                            shuttleDest1(fasterShuttle(1), event(1,2)) = 1;
                            shuttlePickPass1(fasterShuttle(1), event(1,2),1) = shuttlePickPass1(fasterShuttle(1), event(1,2),1) + 1;
                            shuttlePickPass1(fasterShuttle(1), event(1,2),1+shuttlePickPass1(fasterShuttle(1), event(1,2),1)) = time;                            
                            error1(fasterShuttle(1), event(1,2),1+shuttlePickPass1(fasterShuttle(1), event(1,2),1)) = 12;
                            shuttleDest1(fasterShuttle(1), event(1,3)) = 1;
                        else
                            shuttlePickPass2(fasterShuttle(1), event(1,2),1) = shuttlePickPass2(fasterShuttle(1), event(1,2),1) + 1;
                            shuttlePickPass2(fasterShuttle(1), event(1,2),1+shuttlePickPass2(fasterShuttle(1), event(1,2),1)) = time;                            
                            error2(fasterShuttle(1), event(1,2),1+shuttlePickPass1(fasterShuttle(1), event(1,2),1)) = 29;
                            shuttleDest1(fasterShuttle(1), event(1,3)) = 1;
                        end
                    %Passenger_pos < last_shuttle_backwards_dest (don't happen)
                    else
                        shuttleDest2(fasterShuttle(1), event(1,2)) = 1;
                        shuttlePickPass2(fasterShuttle(1), event(1,2),1) = shuttlePickPass2(fasterShuttle(1), event(1,2),1) + 1;
                        shuttlePickPass2(fasterShuttle(1), event(1,2),1+shuttlePickPass2(fasterShuttle(1), event(1,2),1)) = time;
                        error2(fasterShuttle(1), event(1,2),1+shuttlePickPass2(fasterShuttle(1), event(1,2),1)) = 21;
                        shuttleDest1(fasterShuttle(1), event(1,3)) = 1; 
                    end
                end
            %Passenger moving backwards
            else
                %Shuttle moving forwards
                if (shuttlePosition(fasterShuttle(1),4) == 1)                    
                    %Determine last destination forwards of the shuttle
                    posaux2 = 0;
                    for i=1:numberOfStops
                        if (shuttleDest1(fasterShuttle(1),i) == 1)
                            posaux2 = i;
                        end
                    end
                    %Passenger_pos >= last_shuttle_forwards_dest
                    if (posaux2 >= event(1,2) && posaux2 == numberOfStops)
                        if (posaux2 ~= event(1,2))
                            shuttleDest2(fasterShuttle(1), event(1,2)) = 1;
                            shuttlePickPass2(fasterShuttle(1), event(1,2),1) = shuttlePickPass2(fasterShuttle(1), event(1,2),1) + 1;
                            shuttlePickPass2(fasterShuttle(1), event(1,2),1+shuttlePickPass2(fasterShuttle(1), event(1,2),1)) = time;
                            error2(fasterShuttle(1), event(1,2),1+shuttlePickPass2(fasterShuttle(1), event(1,2),1)) = 22;
                            shuttleDest2(fasterShuttle(1), event(1,3)) = 1;
                        else
                            shuttlePickPass1(fasterShuttle(1), event(1,2),1) = shuttlePickPass1(fasterShuttle(1), event(1,2),1) + 1;
                            shuttlePickPass1(fasterShuttle(1), event(1,2),1+shuttlePickPass1(fasterShuttle(1), event(1,2),1)) = time;
                            error1(fasterShuttle(1), event(1,2),1+shuttlePickPass1(fasterShuttle(1), event(1,2),1)) = 19;
                            shuttleDest2(fasterShuttle(1), event(1,3)) = 1;    
                        end
                    %Passenger_pos > last_shuttle_forwards_dest (don't happen)
                    else
                        shuttleDest1(fasterShuttle(1), event(1,2)) = 1;
                        shuttlePickPass1(fasterShuttle(1), event(1,2),1) = shuttlePickPass1(fasterShuttle(1), event(1,2),1) + 1;
                        shuttlePickPass1(fasterShuttle(1), event(1,2),1+shuttlePickPass1(fasterShuttle(1), event(1,2),1)) = time;                            
                        error1(fasterShuttle(1), event(1,2),1+shuttlePickPass1(fasterShuttle(1), event(1,2),1)) = 13;
                        shuttleDest2(fasterShuttle(1), event(1,3)) = 1;
                    end
                end
                %Shuttle moving backwards
                if (shuttlePosition(fasterShuttle(1),4) == -1) 
                    %Shuttle moving backwards and shuttle_pos >= passenger_pos
                    if (stopsPosition(event(1,2)) <= pos)
                        shuttleDest2(fasterShuttle(1), event(1,2)) = 1;
                        shuttlePickPass2(fasterShuttle(1), event(1,2),1) = shuttlePickPass2(fasterShuttle(1), event(1,2),1) + 1;
                        shuttlePickPass2(fasterShuttle(1), event(1,2),1+shuttlePickPass2(fasterShuttle(1), event(1,2),1)) = time;
                        error2(fasterShuttle(1), event(1,2),1+shuttlePickPass2(fasterShuttle(1), event(1,2),1)) = 23;
                        shuttleDest2(fasterShuttle(1), event(1,3)) = 1;
                        if (abs(stopsPosition(event(1,2)) - pos) < abs(stopsPosition(shuttlePosition(fasterShuttle(1),3)) - pos))
                            shuttlePosition(fasterShuttle(1),3) = event(1,2); 
                        end
                    %Shuttle moving backwards and shuttle_pos < passenger_pos (don't happen)
                    else
                        test = 2;
                    end
                end
            end
            shuttleArrivalTime(1,fasterShuttle(1)) = time + abs(stopsPosition(shuttlePosition(fasterShuttle(1),3)) - pos)/shuttleAverageSpeed;
        %Shuttle not moving, but serving passengers
        elseif (shuttlePosition(fasterShuttle(1),4) == 0 && shuttleNumberPass(1,fasterShuttle(1)) ~= 0) 
            shuttlePosition(fasterShuttle(1),2) = time;
            shuttleIdleTime(1,fasterShuttle(1)) = 0;
            shuttleArrivalTime(1,fasterShuttle(1)) = 0;
            shuttleNumberPass(1,fasterShuttle(1)) = shuttleNumberPass(1,fasterShuttle(1)) + 1; %[shuttle totalNumberPass] 
            %Passenger moving forwards
            if (event(1,2) < event(1,3))
                %Shuttle heading forwards
                if (pos < stopsPosition(shuttlePosition(fasterShuttle(1),3)))
                    %Shuttle_pos < Passenger_pos
                    if (pos < stopsPosition(event(1,2)))
                        shuttleDest1(fasterShuttle(1), event(1,2)) = 1; %binary variable with every destination of each shuttle
                        shuttlePickPass1(fasterShuttle(1), event(1,2),1) = shuttlePickPass1(fasterShuttle(1), event(1,2),1) + 1;
                        shuttlePickPass1(fasterShuttle(1), event(1,2),1+shuttlePickPass1(fasterShuttle(1), event(1,2),1)) = time;                            
                        error1(fasterShuttle(1), event(1,2),1+shuttlePickPass1(fasterShuttle(1), event(1,2),1)) = 14;
                        shuttleDest1(fasterShuttle(1), event(1,3)) = 1;
                        if (abs(stopsPosition(event(1,3)) - pos) < abs(stopsPosition(shuttlePosition(fasterShuttle(1),3)) - pos))
                            shuttlePosition(fasterShuttle(1),3) = event(1,3); 
                        end
                    %Shuttle_pos > Passenger_pos (don't happen)
                    elseif (pos > stopsPosition(event(1,2)))
                        test = 3;
                    %Shuttle_pos = Passenger_pos
                    elseif (pos == stopsPosition(event(1,2)))
                        aux = size(waitingTime,1) + 1;
                        waitingTime(aux,1) = time;
                        waitingTime(aux,2) = time;
                        waitingTime(aux,3) = waitingTime(aux,2) - waitingTime(aux,1);
                        waitingTime(aux,4) = event(1,2);
                        shuttleLeavingTime(1,fasterShuttle(1)) = time + stopTime;
                        shuttleDest1(fasterShuttle(1), event(1,3)) = 1; %binary variable with every destination of each shuttle
                        if (abs(stopsPosition(event(1,3)) - pos) < abs(stopsPosition(shuttlePosition(fasterShuttle(1),3)) - pos))
                            shuttlePosition(fasterShuttle(1),3) = event(1,3); 
                        end
                    end
                %Shuttle heading backwards
                else
                    %Determine last destination backwards of the shuttle
                    posaux2 = 0;
                    auxbin = 0;
                    while (auxbin == 0 && posaux2 < numberOfStops)
                        posaux2 = posaux2 + 1;
                        if (shuttleDest2(fasterShuttle(1),posaux2) == 1)
                            auxbin = 1;
                        end
                    end
                    %Passenger_pos >= last_shuttle_backwards_dest
                    if (posaux2 <= event(1,2) && posaux2 == 1)
                        if (posaux2 ~= event(1,2))
                            shuttleDest1(fasterShuttle(1), event(1,2)) = 1; %binary variable with every destination of each shuttle
                            shuttlePickPass1(fasterShuttle(1), event(1,2),1) = shuttlePickPass1(fasterShuttle(1), event(1,2),1) + 1;
                            shuttlePickPass1(fasterShuttle(1), event(1,2),1+shuttlePickPass1(fasterShuttle(1), event(1,2),1)) = time;
                            error1(fasterShuttle(1), event(1,2),1+shuttlePickPass1(fasterShuttle(1), event(1,2),1)) = 15;
                            shuttleDest1(fasterShuttle(1), event(1,3)) = 1;
                        else
                            shuttlePickPass2(fasterShuttle(1), event(1,2),1) = shuttlePickPass2(fasterShuttle(1), event(1,2),1) + 1;
                            shuttlePickPass2(fasterShuttle(1), event(1,2),1+shuttlePickPass2(fasterShuttle(1), event(1,2),1)) = time;
                            error2(fasterShuttle(1), event(1,2),1+shuttlePickPass2(fasterShuttle(1), event(1,2),1)) = 30;
                            shuttleDest1(fasterShuttle(1), event(1,3)) = 1;
                        end
                    %Passenger_pos < last_shuttle_backwards_dest (don't happen)
                    else
                        shuttleDest2(fasterShuttle(1), event(1,2)) = 1; %binary variable with every destination of each shuttle
                        shuttlePickPass2(fasterShuttle(1), event(1,2),1) = shuttlePickPass2(fasterShuttle(1), event(1,2),1) + 1;
                        shuttlePickPass2(fasterShuttle(1), event(1,2),1+shuttlePickPass2(fasterShuttle(1), event(1,2),1)) = time;
                        error2(fasterShuttle(1), event(1,2),1+shuttlePickPass2(fasterShuttle(1), event(1,2),1)) = 24;
                        shuttleDest1(fasterShuttle(1), event(1,3)) = 1;
                    end                        
                end
            %Passenger moving backwards
            else
                %Shuttle heading backwards
                if (pos > stopsPosition(shuttlePosition(fasterShuttle(1),3)))
                    %Shuttle_pos > Passenger_pos
                    if (pos > stopsPosition(event(1,2)))
                        shuttleDest2(fasterShuttle(1), event(1,2)) = 1; %binary variable with every destination of each shuttle 
                        shuttlePickPass2(fasterShuttle(1), event(1,2),1) = shuttlePickPass2(fasterShuttle(1), event(1,2),1) + 1;
                        shuttlePickPass2(fasterShuttle(1), event(1,2),1+shuttlePickPass2(fasterShuttle(1), event(1,2),1)) = time;
                        error2(fasterShuttle(1), event(1,2),1+shuttlePickPass2(fasterShuttle(1), event(1,2),1)) = 25;
                        shuttleDest2(fasterShuttle(1), event(1,3)) = 1;
                        if (abs(stopsPosition(event(1,3)) - pos) < abs(stopsPosition(shuttlePosition(fasterShuttle(1),3)) - pos))
                            shuttlePosition(fasterShuttle(1),3) = event(1,3); 
                        end
                    %Shuttle_pos < Passenger_pos (don't happen)
                    elseif (pos < stopsPosition(event(1,2)))
                        test = 4;
                    %Shuttle_pos = Passenger_pos
                    elseif (pos == stopsPosition(event(1,2)))
                        aux = size(waitingTime,1) + 1;
                        waitingTime(aux,1) = time;
                        waitingTime(aux,2) = time;
                        waitingTime(aux,3) = waitingTime(aux,2) - waitingTime(aux,1);
                        waitingTime(aux,4) = event(1,2);
                        shuttleLeavingTime(1,fasterShuttle(1)) = time + stopTime;
                        shuttleDest2(fasterShuttle(1), event(1,3)) = 1; %binary variable with every destination of each shuttle 
                        if (abs(stopsPosition(event(1,3)) - pos) < abs(stopsPosition(shuttlePosition(fasterShuttle(1),3)) - pos))
                            shuttlePosition(fasterShuttle(1),3) = event(1,3); 
                        end                            
                    end
                %Shuttle heading forwards
                else
                    %Determine last destination forwards of the shuttle
                    for j=1:numberOfStops
                        if (shuttleDest1(fasterShuttle(1),j) == 1)
                            posaux2 = j;
                        end
                    end
                    %Passenger_pos < last_shuttle_forwards_dest
                    if (posaux2 >= event(1,2)) 
                        shuttleDest2(fasterShuttle(1), event(1,2)) = 1; %binary variable with every destination of each shuttle
                        shuttlePickPass2(fasterShuttle(1), event(1,2),1) = shuttlePickPass2(fasterShuttle(1), event(1,2),1) + 1;
                        shuttlePickPass2(fasterShuttle(1), event(1,2),1+shuttlePickPass2(fasterShuttle(1), event(1,2),1)) = time;
                        error2(fasterShuttle(1), event(1,2),1+shuttlePickPass2(fasterShuttle(1), event(1,2),1)) = 26;
                        shuttleDest2(fasterShuttle(1), event(1,3)) = 1;
                    %Passenger_pos < last_shuttle_forwards_dest
                    else 
                        shuttleDest1(fasterShuttle(1), event(1,2)) = 1; %bina25ry variable with every destination of each shuttle
                        shuttlePickPass1(fasterShuttle(1), event(1,2),1) = shuttlePickPass1(fasterShuttle(1), event(1,2),1) + 1;
                        shuttlePickPass1(fasterShuttle(1), event(1,2),1+shuttlePickPass1(fasterShuttle(1), event(1,2),1)) = time;
                        error1(fasterShuttle(1), event(1,2),1+shuttlePickPass1(fasterShuttle(1), event(1,2),1)) = 16;
                        shuttleDest2(fasterShuttle(1), event(1,3)) = 1;
                    end
                end
            end
        %Shuttle serving no passenger
        elseif (shuttleNumberPass(1,fasterShuttle(1)) == 0)                
            shuttleNumberPass(1,fasterShuttle(1)) = shuttleNumberPass(1,fasterShuttle(1)) + 1; %[shuttle totalNumberPass] 
            %Shuttle IDLE where passenger arrived| No need for a new shuttle
            if (pos == stopsPosition(event(1,2)))
                aux = size(waitingTime,1) + 1;
                waitingTime(aux,1) = time;
                waitingTime(aux,2) = time;
                waitingTime(aux,3) = waitingTime(aux,2) - waitingTime(aux,1);
                waitingTime(aux,4) = event(1,2);
                shuttlePosition(fasterShuttle(1),2) = time;
                shuttlePosition(fasterShuttle(1),3) = event(1,3); 
                shuttleIdleTime(1,fasterShuttle(1)) = 0;
                shuttleArrivalTime(1,fasterShuttle(1)) = 0;
                shuttleLeavingTime(1,fasterShuttle(1)) = time + stopTime;
                if (event(1,2) < event(1,3))
                    shuttleDest1(fasterShuttle(1), event(1,3)) = 1; %binary variable with every destination of each shuttle
                else
                    shuttleDest2(fasterShuttle(1), event(1,3)) = 1; %binary variable with every destination of each shuttle
                end
            %Shuttle IDLE in a different position of where passenger arrived| No need for a new shuttle    
            else
                shuttlePosition(fasterShuttle(1),1) = pos; 
                shuttlePosition(fasterShuttle(1),2) = time;
                shuttlePosition(fasterShuttle(1),3) = event(1,2);
                if (event(1,2) < pos)
                    shuttlePosition(fasterShuttle(1),4) = -1; %binary variable with every destination of each shuttle
                else
                    shuttlePosition(fasterShuttle(1),4) = 1; %binary variable with every destination of each shuttle
                end
                shuttleIdleTime(1,fasterShuttle(1)) = 0;
                shuttleArrivalTime(1,fasterShuttle(1)) = time + abs(stopsPosition(event(1,2)) - pos)/shuttleAverageSpeed;
                shuttleLeavingTime(1,fasterShuttle(1)) = 0;
                if (event(1,2) < pos)                        
                    shuttleDest2(fasterShuttle(1), event(1,2)) = 1; %binary variable with every destination of each shuttle
                    shuttlePickPass2(fasterShuttle(1), event(1,2),1) = shuttlePickPass2(fasterShuttle(1), event(1,2),1) + 1;
                    shuttlePickPass2(fasterShuttle(1), event(1,2),1+shuttlePickPass2(fasterShuttle(1), event(1,2),1)) = time;
                    error2(fasterShuttle(1), event(1,2),1+shuttlePickPass2(fasterShuttle(1), event(1,2),1)) = 27;
                    if (event(1,2) < event(1,3))
                        shuttleDest1(fasterShuttle(1), event(1,3)) = 1; %binary variable with every destination of each shuttle;
                    else
                        shuttleDest2(fasterShuttle(1), event(1,3)) = 1; %binary variable with every destination of each shuttle;
                    end
                end
            end
        end
    end
            
    % Event 2 (Pick/Drop Passenger)
    if (event(1,1) == 2)
        shuttlePosition(event(1,2),1) = stopsPosition(shuttlePosition(event(1,2),3));  
        posNumber = shuttlePosition(event(1,2),3);
        shuttlePosition(event(1,2),2) = time;
        shuttleArrivalTime(1,event(1,2)) = 0;
        if (shuttlePosition(event(1,2),4) == 1)
            shuttleDest1(event(1,2),shuttlePosition(event(1,2),3)) = 0; %binary variable with every destination of each shuttle
        else
            shuttleDest2(event(1,2),shuttlePosition(event(1,2),3)) = 0; %binary variable with every destination of each shuttle
        end
        %Determine next destination
        auxbin = 0;
        nextDest=0;
        if (shuttlePosition(event(1,2),4) == 1)
            shuttleNumberPass(1,event(1,2)) = shuttleNumberPass(1,event(1,2)) - shuttleDestPass1(event(1,2),shuttlePosition(event(1,2),3),1); %[shuttle totalNumberPass]
            aux1 = shuttleDestPass1(event(1,2),shuttlePosition(event(1,2),3),1);
            aux2 = shuttlePickPass1(event(1,2),shuttlePosition(event(1,2),3),1);
            for i=1:aux1
                aux3 = size(totalTime,1) + 1;
                totalTime(aux3,1) = shuttleDestPass1(event(1,2),shuttlePosition(event(1,2),3),1+i);
                shuttleDestPass1(event(1,2),shuttlePosition(event(1,2),3),1+i) = 0;
                totalTime(aux3,2) = time;
                totalTime(aux3,3) = totalTime(aux3,2) - totalTime(aux3,1);
                totalTime(aux3,4) = auxTrack1(event(1,2),shuttlePosition(event(1,2),3),1+i);
                auxTrack1(event(1,2),shuttlePosition(event(1,2),3),1+i) = 0;
                totalTime(aux3,5) = posNumber;
            end
            for i=1:aux2
                aux4 = size(waitingTime,1) + 1;
                waitingTime(aux4,1) = shuttlePickPass1(event(1,2),shuttlePosition(event(1,2),3),1+i);
                shuttlePickPass1(event(1,2),shuttlePosition(event(1,2),3),1+i) = 0;
                waitingTime(aux4,2) = time;
                waitingTime(aux4,3) = waitingTime(aux4,2) - waitingTime(aux4,1);
                waitingTime(aux4,4) = posNumber;
            end
            shuttleDestPass1(event(1,2),shuttlePosition(event(1,2),3),1) = 0; %number of passenger going to each destination of each shuttle
            shuttlePickPass1(event(1,2),shuttlePosition(event(1,2),3),1) = 0;
            while (auxbin == 0 && nextDest < numberOfStops)
                nextDest = nextDest+1;
                if (shuttleDest1(event(1,2),nextDest) == 1)
                    auxbin = 1;
                end
            end    
            if (auxbin == 0)
                for i=1:numberOfStops
                    if (shuttleDest2(event(1,2),i) == 1)
                        nextDest = i;
                    end
                end
            end                
        else
            shuttleNumberPass(1,event(1,2)) = shuttleNumberPass(1,event(1,2)) - shuttleDestPass2(event(1,2),shuttlePosition(event(1,2),3),1); %[shuttle totalNumberPass]
            aux1 = shuttleDestPass2(event(1,2),shuttlePosition(event(1,2),3),1);
            aux2 = shuttlePickPass2(event(1,2),shuttlePosition(event(1,2),3),1);
            for i=1:aux1
                aux3 = size(totalTime,1) + 1;
                totalTime(aux3,1) = shuttleDestPass2(event(1,2),shuttlePosition(event(1,2),3),1+i);
                shuttleDestPass2(event(1,2),shuttlePosition(event(1,2),3),1+i) = 0;
                totalTime(aux3,2) = time;
                totalTime(aux3,3) = totalTime(aux3,2) - totalTime(aux3,1);
                totalTime(aux3,4) = auxTrack2(event(1,2),shuttlePosition(event(1,2),3),1+i);
                auxTrack2(event(1,2),shuttlePosition(event(1,2),3),1+i) = 0;
                totalTime(aux3,5) = posNumber;
            end
            for i=1:aux2
                aux4 = size(waitingTime,1) + 1;
                waitingTime(aux4,1) = shuttlePickPass2(event(1,2),shuttlePosition(event(1,2),3),1+i);
                shuttlePickPass2(event(1,2),shuttlePosition(event(1,2),3),1+i) = 0;
                waitingTime(aux4,2) = time;
                waitingTime(aux4,3) = waitingTime(aux4,2) - waitingTime(aux4,1);
                waitingTime(aux4,4) = posNumber;
            end
            shuttleDestPass2(event(1,2),shuttlePosition(event(1,2),3),1) = 0; %number of passenger going to each destination of each shuttle
            shuttlePickPass2(event(1,2),shuttlePosition(event(1,2),3),1) = 0;
            for i=1:numberOfStops
                if (shuttleDest2(event(1,2),i) == 1)
                    nextDest = i;
                end
            end
            if (nextDest == 0)
                while (auxbin == 0 && nextDest < numberOfStops)
                    nextDest = nextDest+1;
                    if (shuttleDest1(event(1,2),nextDest) == 1)
                        auxbin = 1;
                    end
                end 
            end
        end
        shuttlePosition(event(1,2),3) = nextDest; 
        shuttlePosition(event(1,2),4) = 0;  
        shuttleIdleTime(1,event(1,2)) = 0;
        shuttleLeavingTime(1,event(1,2)) = time + stopTime;
        if (shuttleNumberPass(1,event(1,2)) == 0)
            %Decide next destination for the shuttle to stay IDLE
            contaux = (numberOfStops + 1 - mod((numberOfStops+1),2))/2;
            cont = zeros(1, contaux); %counter of the number of shuttles per garage
            for i=1:numberShuttles
                for j=1:contaux
                    pos = shuttlePosition(i,1) + shuttlePosition(i,4)*shuttleAverageSpeed*(time - shuttlePosition(i,2));
                    if(pos == stopsPosition(j + numberOfStops) && shuttleNumberPass(1,i) == 0)
                        cont(1,j) = cont(1,j) + 1;
                    elseif (shuttlePosition(i,3) == j + numberOfStops && shuttleNumberPass(1,i) == 0)
                        cont(1,j) = cont(1,j) + 1;
                    end
                end
            end
            nextDest = numberOfStops + 1;
            contaux2 = stopsPosition(numberOfStops);
            for i=1:contaux
                if (cont(1,i) <= cont(1,nextDest - numberOfStops))
                    if (cont(1,i) < cont(1,nextDest - numberOfStops))
                        contaux2 = abs(stopsPosition(shuttlePosition(event(1,2),3)) - stopsPosition(i + numberOfStops));
                        nextDest = i + numberOfStops;
                    elseif (abs(stopsPosition(shuttlePosition(event(1,2),3)) - stopsPosition(i + numberOfStops)) < contaux2)
                        contaux2 = abs(stopsPosition(shuttlePosition(event(1,2),3)) - stopsPosition(i + numberOfStops));    
                        nextDest = i + numberOfStops;
                    end
                end
            end              
            %Shuttle already in the garage
            if (shuttlePosition(event(1,2),1) == stopsPosition(nextDest))
                shuttlePosition(event(1,2),3) = 0; 
                shuttlePosition(event(1,2),4) = 0;  
                shuttleIdleTime(1,event(1,2)) = time + stopTime;
                shuttleLeavingTime(1,event(1,2)) = 0;
            %Shuttle isn't in the garage
            else
                shuttlePosition(event(1,2),3) = nextDest; 
                shuttlePosition(event(1,2),4) = 0;  
                shuttleIdleTime(1,event(1,2)) = 0;
                shuttleLeavingTime(1,event(1,2)) = time + stopTime;
            end
        end        
    end
    
    % Event 3 (Shuttle IDLE)
    if (event(1,1) == 3)
        if (shuttlePosition(event(1,2),3) ~= 0)
            shuttlePosition(event(1,2),1) = stopsPosition(shuttlePosition(event(1,2),3)); 
        end
        shuttlePosition(event(1,2),2) = time;
        shuttlePosition(event(1,2),3) = 0; 
        shuttlePosition(event(1,2),4) = 0; 
        shuttleIdleTime(1,event(1,2)) = 0;
        shuttleArrivalTime(1,event(1,2)) = 0;
        shuttleLeavingTime(1,event(1,2)) = 0;
    end
    
    % Event 4 (Shuttle leaving station)
    if (event(1,1) == 4)  
        %Going to the garage 
        if (shuttleNumberPass(1,event(1,2)) == 0) 
            shuttlePosition(event(1,2),2) = time;
            shuttlePosition(event(1,2),4) = (shuttlePosition(event(1,2),3) - shuttlePosition(event(1,2),1))/(abs(shuttlePosition(event(1,2),3) - shuttlePosition(event(1,2),1))); 
            shuttleIdleTime(1,event(1,2)) = time + abs(shuttlePosition(event(1,2),1) - shuttlePosition(event(1,2),3))/shuttleAverageSpeed;
            shuttleArrivalTime(1,event(1,2)) = 0;
            shuttleLeavingTime(1,event(1,2)) = 0;
        %Shuttle going to pick/drop more passengers
        else
            shuttlePosition(event(1,2),2) = time;
            shuttlePosition(event(1,2),4) = (stopsPosition(shuttlePosition(event(1,2),3)) - shuttlePosition(event(1,2),1))/(abs(stopsPosition(shuttlePosition(event(1,2),3)) - shuttlePosition(event(1,2),1))); 
            shuttleIdleTime(1,event(1,2)) = 0;
            shuttleArrivalTime(1,event(1,2)) = time + abs(shuttlePosition(event(1,2),1) - stopsPosition(shuttlePosition(event(1,2),3)))/shuttleAverageSpeed;
            shuttleLeavingTime(1,event(1,2)) = 0;
        end
    end    
end

aux1 = size(totalTime,1);
auxsum1 = 0;
cont1 = 0;
for i=2:aux1
    if (totalTime(i,1) >= 0.5 && totalTime(i,1) <= 1.5)
        auxsum1 = auxsum1 + totalTime(i,3);
        cont1 = cont1 + 1;
    end
end
averageTotalTime = auxsum1/(cont1)*60;

aux2 = size(waitingTime,1);
auxsum2 = 0;
cont2 = 0;
for i=2:aux2
    if (waitingTime(i,1) >= 0.5 && waitingTime(i,1) <= 1.5)
        auxsum2 = auxsum2 + waitingTime(i,3);
        cont2 = cont2 + 1;
    end
end
averageWaitingTime = auxsum2/(cont2)*60;
    
travelAux = 1;
waitingAux = waitingTime;
for i=2:aux1
    cont = 0;
    j = 2;
    while (j < aux2 && cont == 0)
        if (waitingAux(j,1) == totalTime(i,1) && waitingAux(j,4) == totalTime(i,4))
            travelAux = travelAux + 1;
            travelTime(travelAux,1) = waitingAux(j,2);
            travelTime(travelAux,2) = totalTime(i,2);
            travelTime(travelAux,3) = totalTime(i,3) - waitingAux(j,3);
            travelTime(travelAux,4) = totalTime(i,4);
            travelTime(travelAux,5) = totalTime(i,5);
            waitingAux(j,1) = -1;
            cont = 1;            
        end
        j = j + 1;
    end
end
aux3 = size(travelTime,1);
auxsum3 = 0;
cont3 = 0;
for i=2:aux3
    if (travelTime(i,1) >= 0.5 && travelTime(i,1) <= 1.5)
        auxsum3 = auxsum3 + travelTime(i,3);
        cont3 = cont3 + 1;
    end
end
averageTravelTime = auxsum3/(cont3)*60;


output = [numberShuttles contSimulation contSimulation/numberShuttles averageWaitingTime averageTravelTime averageTotalTime];
%fprintf('Number of Shuttles:            %-8.0f\n',output(1));
%fprintf('Number of Passengers:          %-8.0f\n',output(2));
%fprintf('Passengers/Shuttle:            %-8.2f\n',output(3));
%fprintf('Average Waiting Time (min):    %-8.2f\n',output(4));
%fprintf('Average Travel Time (min):     %-8.2f\n',output(5));
%fprintf('Average Total Time (min):      %-8.2f\n',output(6));

