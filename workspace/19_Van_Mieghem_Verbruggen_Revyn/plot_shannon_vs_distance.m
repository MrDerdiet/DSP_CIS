addpath(genpath('helper_functions'));
clearvars; hold on; close all;

% volume 100%
% fs = 16000;
% dftsize = 128;
% noverlap = 4;  

c_data = [((1.310688545773893e+05)+(1.292250002811995e+05)+(1.315119748840791e+05))/3 0;
    ((1.136812395781240e+05)+(1.136051908115632e+05)+(1.133296514747234e+05))/3 5;
    ((1.009937377319820e+05)+(1.005348748028261e+05)+(1.002723927183364e+05))/3 10;
    ((8.893546174718878e+04)+(8.951044424842444e+04)+(8.553841234537572e+04))/3 20;
    ((8.684170894328842e+04)+(8.689789193129580e+04)+(8.680191767807242e+04))/3 25;
    ((6.181187536832097e+04)+(6.041927382383488e+04)+(6.104544624616882e+04))/3 55;
    ((5.147106259752672e+04)+(5.089394013463632e+04)+(5.067047814401254e+04))/3 95;
    ((4.800256053985221e+04)+(4.633677408691219e+04)+(4.502204018283998e+04))/3 120];
capacity = c_data(:,1);
capacity_log = log(c_data(:,1));
distance = c_data(:,2);

x1 = linspace(min(distance)-2,max(distance)+1);
p1 = polyfit(distance,capacity,2);
f1 = polyval(p1,x1);
p2 = polyfit(distance,capacity,3);
f2 = polyval(p2,x1);
p3 = polyfit(distance,capacity_log,3);
f3 = polyval(p3,x1);

figure('Name', 'Capacity vs Distance');
title( 'Capacity vs Distance' ); xlabel( 'Distance(cm)' ); ylabel( 'Channel capacity (/)' );
hold on
plot(distance,capacity,'o')
plot(x1,f1,'r-','LineWidth',1)
plot(x1,f2,'b-','LineWidth',1)
plot(x1,exp(f3),'m-','LineWidth',2)
legend('data','2nd order','3rd order','logarhitmic')

