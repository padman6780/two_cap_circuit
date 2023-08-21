function Out = cap_v1v2(t, v1, v2, c1, c2, R, dc1)

% The differential equations to solve the two-capacitor circuit is given 
dv1 = 1/c1(t) * ( -(dc1(t) + 1/R)* v1 + v2/R );
dv2 = (v1 - v2) / (c2*R);

% The derivatives of state variables (voltages at nodes: v1 and v2)
Out = [dv1; dv2];


end