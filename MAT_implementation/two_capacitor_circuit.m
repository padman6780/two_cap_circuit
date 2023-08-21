% Code to simulate the two-capacitor circuit
% Gnd --)|--o--\/\/\/--o--|(---Gnd
%       C1  V1    R    V2 C2
% One of the capacitance (C1) is varied with a user defined function
% Steady state is V1 = V2
% Intuition is that (Q1=C1*V1) + (Q2=C2*V2) is conserved

load('v1v2.mat');

tbegin = 1; % 1 ms
tdur = 2.75; % 2.75 ms
convert = 0;
if convert
    abs_cap = 1e8/pi;
    abs_R = 4e-3/pi;
else
    abs_cap = 1;
    abs_R = 1;
end
c1 = @(t) abs_cap*(1 + 0.5 * (t > tbegin & t < (tbegin+tdur) ) .* sin(2*pi*t));
dc1 = @(t)  abs_cap*(0.5 * 2 * pi * (t > tbegin & t < (tbegin+tdur) ) .* cos(2*pi*t));
c2 = abs_cap*(1);
R = abs_R * (1e0);
v10 = -65;
v20 = -65;

dt = 1/(40);
tSim = 0:dt:5;

plot(tSim, c1(tSim));
% figure(2);
% plot(tSim, dc1(tSim));


X0 = [v10 v20];
OdeOpts=odeset('MaxStep',dt,'AbsTol',1e-3,'RelTol',1e-3);
[t,X] = ode113(@(t,X) cap_v1v2(t,X(1),X(2),c1, c2, R, dc1),[0 5],X0,OdeOpts); % --> External function cap_v1v2


figure(3);
plot(t, X(:,1));
hold on
plot(v1nrn.t, v1nrn.v1, 'r');
ylim([-150, 0]);
legend('matlab', 'neuron');
ylabel('V1');

figure(4);
plot(t, X(:,2));
hold on;
plot(v2nrn.t, v2nrn.v2, 'r');
ylim([-150, 0]);
legend('matlab', 'neuron');
ylabel('V2');


% Total Charge
Q = c1(t).*X(:,1) + X(:,2);
Q_nrn = c1(v1nrn.t).*v1nrn.v1 + v2nrn.v2;
figure(5);
plot(t, Q);
hold on;
plot(v1nrn.t, Q_nrn, 'r');
legend('matlab', 'neuron');
ylabel('Charge (Q1 + Q2)');

figure(6);
plot(t, Q);
ylim([-130.1 -129.9]);
hold on;
plot(v1nrn.t, Q_nrn, 'r');
legend('matlab', 'neuron');
ylabel('Charge (Q1 + Q2)');


%t, v1, v2, c1, c2, R, dc1

