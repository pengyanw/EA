clear; clc; close all;

% initial condition (large angle)
x0 = [0.5; 0.0];  % theta (rad), dtheta

% define polynomial controller as a function handle
ufun = @(x) poly_ctrl(x);

% simulate
tspan = [0 10];
[t, X] = ode45(@(t,x) pendulum_dyn(t, x, ufun), tspan, x0);

% compute control for plotting
U = arrayfun(@(i) ufun(X(i,:)'), 1:length(t));

% plot
figure;
subplot(3,1,1); plot(t, X(:,1)); ylabel('\theta (rad)'); grid on;
subplot(3,1,2); plot(t, X(:,2)); ylabel('d\theta/dt (rad/s)'); grid on;
subplot(3,1,3); plot(t, U); ylabel('u'); xlabel('time (s)'); grid on;
title('Polynomial Feedback Control on Nonlinear Pendulum');
