function simulationBQX1_v2

for i=1:10
    for j=1:10
        aux = BQX_Project1FINAL_v2(50*i);
        waitingTime(j,i) = aux(4);
        travelTime(j,i) = aux(5);
        totalTime(j,i) = aux(6);
    end
    aux2(i) = 50*i;
end

figure (1);
boxplot(waitingTime,aux2);
xlabel('Fleet Size');
ylabel('Average Waiting Time (min)');
title('Average waiting time for different fleet sizes');

figure (2);
boxplot(travelTime,aux2);
xlabel('Fleet size');
ylabel('Average travel time (min)');
title('Average travel time for different fleet sizes');

figure (3);
boxplot(totalTime,aux2);
xlabel('Fleet size');
ylabel('Average total time (min)');
title('Average total time for different fleet sizes');


end


