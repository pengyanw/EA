function out = normalmat(K)
[Nu, Nx] = size(K);
ratio = floor(Nx/Nu);
X = zeros(Nx, Nu);
for i = 1:Nu
    X((ratio*i-1):(ratio*(i+1)-2),i) = ones(ratio,1);
end
out = K*X;
end