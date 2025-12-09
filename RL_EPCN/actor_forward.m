function Z = actor_forward(actor, S)
    % S: sDim × batch
    % Z: dAct × batch

    h1 = actor.W1 * S + actor.b1;     % (128 × batch)
    a1 = max(h1, 0);

    h2 = actor.W2 * a1 + actor.b2;    % (128 × batch)
    a2 = max(h2, 0);

    Z  = actor.W3 * a2 + actor.b3;    % (dAct × batch)
end
