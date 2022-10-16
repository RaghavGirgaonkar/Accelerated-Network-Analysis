function y = genqc(x, A, coeffs, initial_phase)
%Create a simulated quadratic chirp y = A*sin(ax + bx^2 + cx^3) for a
%specified snr
phase = coeffs(1)*x + coeffs(2)*x.^2 + coeffs(3)*x.^3 + initial_phase;
sigy = sin(2*pi*phase);
y = A*sigy/norm(sigy);
