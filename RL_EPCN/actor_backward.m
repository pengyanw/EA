function grad = actor_backward(actor, critic, S, Z)
    batch = size(S,2);

    % 初始化所有字段的梯度为 0（和 actor 的参数结构相同）
    grad = zero_like_actor(actor);

    for t = 1:batch
        s = S(:,t);
        z = Z(:,t);

        % compute ∂Q/∂z for sample t
        dQdz = finite_diff_q(critic, s, z);

        % backprop through actor network
        grad_t = backprop_actor_single(actor, s, z, dQdz);

        % accumulate
        grad = add_grad(grad, grad_t);
    end

    % average
    grad = scale_grad(grad, 1/batch);
end
