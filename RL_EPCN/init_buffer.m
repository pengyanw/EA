function buf = init_buffer(N, sDim, dAct)
    buf.s = zeros(sDim, N);
    buf.z = zeros(dAct, N);
    buf.r = zeros(1, N);
    buf.count = 0;
    buf.maxN = N;
end


