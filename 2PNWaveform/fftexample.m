a = 1:10;
A = zeros(1,6);

for i = 1:6
    sum = 0;
    for j = 0:9
        sum = sum + a(j+1)*exp(-2*pi*1i*(i-1)*j/10);
    end
    A(i) = sum;
end

% Aneg = zeros(1,4);

Aneg = conj(A(5:-1:2));

finalA = [A, Aneg]

ifft(finalA)