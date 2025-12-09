function dx = pendulum_dyn(x, u)
    g = 9.81;
    L = 0.5;
    m = 0.15;
    b = 0.1;

    dx = zeros(2,1);
    dx(1) = x(2);
    dx(2) = (m*g*L*sin(x(1)) - b*x(2))/(m*L^2) + u/(m*L^2);
end
