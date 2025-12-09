function add_transition(buf, s, z, r)
    idx = mod(buf.count, buf.maxN) + 1;
    buf.s(:,idx) = s;
    buf.z(:,idx) = z;
    buf.r(idx)   = r;
    buf.count = buf.count + 1;
end

