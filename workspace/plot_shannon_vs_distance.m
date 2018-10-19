addpath(genpath('helper_functions'));
clearvars; hold on; close all;

c_data = [500 10; 503 10; 207 20; 194 20];
capacity = c_data(:,1);
distance = c_data(:,2);

p = polyfit(distance,capacity,2);
x1 = linspace(min(distance),max(distance));
f1 = polyval(p,x1);

figure('Name', 'Capacity vs Distance');
title( 'Capacity vs Distance' ); xlabel( 'Distance(cm)' ); ylabel( 'Channel capacity (/)' );
hold on
plot(distance,capacity,'o')
plot(x1,f1,'r--')
legend('data','function')