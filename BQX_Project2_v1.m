function [ output ] = BQX_Project2_v1 ( fleetSize )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EZ10 Simulation for the rush hour
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all

contSimulation = 0;

% Inputs
lambda = xlsread('StreetcarData', 'Data1'); %arrival rates (passenger per hour)
stopsPosition = xlsread('StreetcarData', 'Data2');  %miles
numberOfStops = size(stopsPosition,2);
shuttleAverageSpeed = 10.6; %miles/hour
stopTime = 15/3600;  %hour
shuttleCapacity = 12;
numberShuttles = fleetSize; %size of the fleet (even number)
simulationTime = 2; %hour

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

passengerServed = zeros(numberOfStops, numberOfStops); %total number of passengers served per station
timeStation = zeros(1,numberOfStops); %time last shuttle arrived in each station
shuttlePosition = zeros(numberShuttles,4);
shuttleArrivalTime = zeros(1,numberShuttles); %time event 2 will occur
shuttleLeavingTime = zeros(1,numberShuttles); %time event 3 will occur
shuttleDestPass1 = zeros(numberShuttles,numberOfStops,1+shuttleCapacity); %number of passenger going to each destination for each shuttle moving fowards
shuttleDestPass2 = zeros(numberShuttles,numberOfStops,1+shuttleCapacity); %number of passenger going to each destination for each shuttle moving backwards
auxTrack1 = zeros(numberShuttles,numberOfStops,1+shuttleCapacity); %auxiliary variable to help track the passengers
auxTrack2 = zeros(numberShuttles,numberOfStops,1+shuttleCapacity); %auxiliary variable to help track the passengers
shuttleNumberPass = int8(zeros(1,numberShuttles)); %[shuttle totalNumberPass]
travelTime = zeros(1,5);
waitingTime = zeros(1,4);
totalTime = zeros(1,5);

for i=1:numberShuttles/2
    shuttlePosition(i,1) = stopsPosition(numberOfStops)*(i-1)/(numberShuttles/2); %position o the shuttle at a specified moment 
    auxbin = 0;
    nextDest = 0;
    while (auxbin == 0 && nextDest < numberOfStops - 1)
        nextDest = nextDest + 1;
        if (stopsPosition(nextDest) <= shuttlePosition(i,1) && stopsPosition(nextDest+1) > shuttlePosition(i,1))
            auxbin = 1;
        end
    end
    nextDest = nextDest + 1;      
    shuttlePosition(i,2) = 0; %time in which shuttlePosition was accounted for
    shuttlePosition(i,3) = nextDest; %destination (next stop of the shuttle)
    shuttlePosition(i,4) = 1; %direction ( -1 means shuttle moving backwards; 0 means that the shuttle isn't moving; 1 means shuttle moving forwards)
    shuttleArrivalTime(1,i) = abs(stopsPosition(nextDest) - shuttlePosition(i,1))/shuttleAverageSpeed; %time event 2 will occur
    shuttleLeavingTime(1,i) = 0; %time event 4 will occur
    shuttleDestPass1(i,:,1) = 0; 
    shuttleDestPass2(i,:,1) = 0;
    auxTrack1(i,:,1) = 0; 
    auxTrack2(i,:,1) = 0;
    shuttleNumberPass(1,i) = 0; %[shuttle totalNumberPass]   
end

for i=(numberShuttles/2 + 1):numberShuttles
    shuttlePosition(i,1) = stopsPosition(numberOfStops)*(i-numberShuttles/2)/(numberShuttles/2); %position o the shuttle at a specified moment 
    auxbin = 0;
    nextDest = 0;
    while (auxbin == 0 && nextDest < numberOfStops-1)
        nextDest = nextDest + 1;
        if (stopsPosition(nextDest) <= shuttlePosition(i,1) && stopsPosition(nextDest+1) > shuttlePosition(i,1))
            auxbin = 1;
        end
    end
    if (auxbin == 0)
        nextDest = numberOfStops;
    end
    if (stopsPosition(nextDest) == shuttlePosition(i,1))
        nextDest = nextDest - 1;    
    end
    shuttlePosition(i,2) = 0; %time in which shuttlePosition was accounted for
    shuttlePosition(i,3) = nextDest; %destination (next stop of the shuttle)
    shuttlePosition(i,4) = -1; %direction ( -1 means shuttle moving backwards; 0 means that the shuttle isn't moving; 1 means shuttle moving forwards)
    shuttleArrivalTime(1,i) = abs(stopsPosition(nextDest) - shuttlePosition(i,1))/shuttleAverageSpeed; %time event 2 will occur
    shuttleLeavingTime(1,i) = 0; %time event 4 will occur
    shuttleDestPass1(i,:,1) = 0; 
    shuttleDestPass2(i,:,1) = 0;
    auxTrack1(i,:,1) = 0; 
    auxTrack2(i,:,1) = 0;
    shuttleNumberPass(1,i) = 0; %[shuttle totalNumberPass]   
end

time = 0;

% Events Simulation
while time < simulationTime
    %Choose event
    aux = simulationTime;
    event = zeros(1,4);    
    for i=1:numberShuttles
        %Passenger pick/drop
        if (shuttleArrivalTime(1,i) < aux && shuttleArrivalTime(1,i) > 0)
            aux = shuttleArrivalTime(1,i);
            event = [1 i 0 0];
        end
        %Shuttle leaving station
        if (shuttleLeavingTime(1,i) < aux && shuttleLeavingTime(1,i) > 0)
            aux = shuttleLeavingTime(1,i);
            event = [2 i 0 0];
        end
    end        
    time = aux;    
            
    % Event 1 (Pick/Drop Passenger)
    if (event(1,1) == 1)
        direction = shuttlePosition(event(1,2),4);
        posNumber = shuttlePosition(event(1,2),3);
        shuttlePosition(event(1,2),1) = stopsPosition(shuttlePosition(event(1,2),3));  
        shuttlePosition(event(1,2),2) = time;
        shuttlePosition(event(1,2),4) = 0;
        shuttleArrivalTime(1,event(1,2)) = 0;
        shuttleLeavingTime(1,event(1,2)) = time + stopTime;
        %Determine next destination
        if (direction == 1)
            if (posNumber < numberOfStops)
                shuttlePosition(event(1,2),3) = shuttlePosition(event(1,2),3) + 1;
            else
                shuttlePosition(event(1,2),3) = shuttlePosition(event(1,2),3) - 1;
            end
        else            
            if (posNumber > 1)
                shuttlePosition(event(1,2),3) = shuttlePosition(event(1,2),3) - 1;
            else
                shuttlePosition(event(1,2),3) = shuttlePosition(event(1,2),3) + 1;
            end
        end
        %Drop Passengers
        if (direction == 1)
            shuttleNumberPass(1,event(1,2)) = shuttleNumberPass(1,event(1,2)) - shuttleDestPass1(event(1,2),posNumber,1); %[shuttle totalNumberPass]
            aux1 = shuttleDestPass1(event(1,2),shuttlePosition(event(1,2),3),1);
            for i=1:aux1
                aux3 = size(travelTime,1) + 1;
                travelTime(aux3,1) = shuttleDestPass1(event(1,2),shuttlePosition(event(1,2),3),1+i);
                shuttleDestPass1(event(1,2),shuttlePosition(event(1,2),3),1+i) = 0;
                travelTime(aux3,2) = time;
                travelTime(aux3,3) = travelTime(aux3,2) - travelTime(aux3,1);
                travelTime(aux3,4) = auxTrack1(event(1,2),shuttlePosition(event(1,2),3),1+i);
                auxTrack1(event(1,2),shuttlePosition(event(1,2),3),1+i) = 0;
                travelTime(aux3,5) = posNumber;
            end
            shuttleDestPass1(event(1,2),shuttlePosition(event(1,2),3),1) = 0; %number of passenger going to each destination of each shuttle            
        else
            shuttleNumberPass(1,event(1,2)) = shuttleNumberPass(1,event(1,2)) - shuttleDestPass2(event(1,2),shuttlePosition(event(1,2),3),1); %[shuttle totalNumberPass]
            aux1 = shuttleDestPass2(event(1,2),shuttlePosition(event(1,2),3),1);
            for i=1:aux1
                aux3 = size(travelTime,1) + 1;
                travelTime(aux3,1) = shuttleDestPass2(event(1,2),shuttlePosition(event(1,2),3),1+i);
                shuttleDestPass2(event(1,2),shuttlePosition(event(1,2),3),1+i) = 0;
                travelTime(aux3,2) = time;
                travelTime(aux3,3) = travelTime(aux3,2) - travelTime(aux3,1);
                travelTime(aux3,4) = auxTrack2(event(1,2),shuttlePosition(event(1,2),3),1+i);
                auxTrack2(event(1,2),shuttlePosition(event(1,2),3),1+i) = 0;
                travelTime(aux3,5) = posNumber;
            end
            shuttleDestPass2(event(1,2),shuttlePosition(event(1,2),3),1) = 0; %number of passenger going to each destination of each shuttle
        end
        %Pick Passengers
        auxbin = 0;
        auxcont = 0;
        auxtime1 = simulationTime;
        auxvector = zeros(1,4);
        seatsAvailable = shuttleCapacity - shuttleNumberPass(1,event(1,2));
        while (auxcont < seatsAvailable)
            for i=1:numberOfStops
                if (passengerServed(posNumber,i)+1 <= numberPassenger(posNumber,i))
                    if (arrivalTime(posNumber,i,passengerServed(posNumber,i)+1) <= shuttleLeavingTime(1,event(1,2)) && arrivalTime(posNumber,i,passengerServed(posNumber,i)+1) <= auxtime1)
                        auxtime1 = arrivalTime(posNumber,i,passengerServed(posNumber,i)+1);
                        auxbin = 1;
                        auxvector = [1 auxtime1 i passengerServed(posNumber,i)+1];
                    end
                end
            end
            if (auxbin == 1)
                auxbin = 0;
                auxcont = auxcont + 1;
                auxtime1 = simulationTime;
                passengerServed(posNumber,auxvector(1,3)) = auxvector(1,4);
                timeStation(1,posNumber) = auxvector(1,2);
                shuttleNumberPass(1,event(1,2)) = shuttleNumberPass(1,event(1,2)) + 1;
                %Update the passengers destinations of the shuttle
                if (posNumber < auxvector(1,3))
                    shuttleDestPass1(event(1,2), auxvector(1,3),1) = shuttleDestPass1(event(1,2), auxvector(1,3),1) + 1; 
                    shuttleDestPass1(event(1,2), auxvector(1,3),1+shuttleDestPass1(event(1,2), auxvector(1,3),1)) = time;
                    auxTrack1(event(1,2), auxvector(1,3),1+shuttleDestPass1(event(1,2), auxvector(1,3),1)) = posNumber;
                else
                    shuttleDestPass2(event(1,2), auxvector(1,3),1) = shuttleDestPass2(event(1,2), auxvector(1,3),1) + 1; 
                    shuttleDestPass2(event(1,2), auxvector(1,3),1+shuttleDestPass2(event(1,2), auxvector(1,3),1)) = time;
                    auxTrack2(event(1,2), auxvector(1,3),1+shuttleDestPass2(event(1,2), auxvector(1,3),1)) = posNumber;
                end
                %Calculate waiting time of the passenger
                aux4 = size(waitingTime,1) + 1;
                waitingTime(aux4,1) = auxvector(1,2);
                waitingTime(aux4,2) = time;   
                waitingTime(aux4,3) = waitingTime(aux4,2) - waitingTime(aux4,1);
                if (waitingTime(aux4,3) <= 0)
                    waitingTime(aux4,3) = 0;
                end
                waitingTime(aux4,4) = posNumber;
            else
                auxcont = shuttleCapacity + 1;
            end
        end
    end
    
    % Event 2 (Shuttle leaving station)
    if (event(1,1) == 2)  
        shuttlePosition(event(1,2),2) = time;
        shuttlePosition(event(1,2),4) = (stopsPosition(shuttlePosition(event(1,2),3)) - shuttlePosition(event(1,2),1))/(abs(stopsPosition(shuttlePosition(event(1,2),3)) - shuttlePosition(event(1,2),1))); 
        shuttleArrivalTime(1,event(1,2)) = time + abs(shuttlePosition(event(1,2),1) - stopsPosition(shuttlePosition(event(1,2),3)))/shuttleAverageSpeed;
        shuttleLeavingTime(1,event(1,2)) = 0;        
    end    
end

aux1 = size(travelTime,1);
auxsum1 = 0;
cont1 = 0;
for i=2:aux1
    if (travelTime(i,1) >= 0.5 && travelTime(i,1) <= 1.5)
        auxsum1 = auxsum1 + travelTime(i,3);
        cont1 = cont1 + 1;
    end
end
averageTravelTime = auxsum1/(cont1)*60;

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
    
totalAux = 1;
waitingAux = waitingTime;
for i=2:aux1
    cont = 0;
    j = 2;
    while (j < aux2 && cont == 0)
        if (waitingAux(j,2) == travelTime(i,1) && waitingAux(j,4) == travelTime(i,4))
            totalAux = totalAux + 1;
            totalTime(totalAux,1) = waitingAux(j,1);
            totalTime(totalAux,2) = travelTime(i,2);
            totalTime(totalAux,3) = travelTime(i,3) + waitingAux(j,3);
            totalTime(totalAux,4) = travelTime(i,4);
            totalTime(totalAux,5) = travelTime(i,5);
            waitingAux(j,2) = 0;
            cont = 1;            
        end
        j = j + 1;
    end
end
aux3 = size(totalTime,1);
auxsum3 = 0;
cont3 = 0;
for i=2:aux3
    if (totalTime(i,1) >= 0.5 && totalTime(i,1) <= 1.5)
        auxsum3 = auxsum3 + totalTime(i,3);
        cont3 = cont3 + 1;
    end
end
averageTotalTime = auxsum3/(cont3)*60;


output = [numberShuttles contSimulation contSimulation/numberShuttles averageWaitingTime averageTravelTime averageTotalTime];
%fprintf('Number of Shuttles:            %-8.0f\n',output(1));
%fprintf('Number of Passengers:          %-8.0f\n',output(2));
%fprintf('Passengers/Shuttle:            %-8.2f\n',output(3));
%fprintf('Average Waiting Time (min):    %-8.2f\n',output(4));
%fprintf('Average Travel Time (min):     %-8.2f\n',output(5));
%fprintf('Average Total Time (min):      %-8.2f\n',output(6));


