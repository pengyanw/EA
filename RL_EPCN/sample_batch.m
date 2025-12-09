function [sBatch, zBatch, rBatch] = sample_batch(buf, B)
    N = min(buf.count, buf.maxN);
    ind = randi(N, [B,1]);
    sBatch = buf.s(:,ind);
    zBatch = buf.z(:,ind);
    rBatch = buf.r(ind);
end