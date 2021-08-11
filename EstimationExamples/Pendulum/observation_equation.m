function z = observation_equation(x, m)
    offset = m(1);
    theta = x(1);
    
    z = offset + theta;
end