function u = poly_ctrl(x)
    % polynomial coefficients
    k1 = -20;
    k2 = -5;
    k3 = -50;
    k4 = -2;

    theta = x(1);
    dtheta = x(2);

    u = k1*theta + k2*dtheta + k3*(theta^3) + k4*(dtheta^3);
end
