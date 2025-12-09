function model = adam_update(model, grad, opt, lr)
    beta1=0.9; beta2=0.999; eps=1e-8;

    fields = fieldnames(grad);
    for k = 1:length(fields)
        name = fields{k};
        if startsWith(name,'W') || startsWith(name,'b')
            g = grad.(name);
            m = getfield(model,['m',name]);
            v = getfield(model,['v',name]);

            m = beta1*m + (1-beta1)*g;
            v = beta2*v + (1-beta2)*(g.^2);

            m_hat = m/(1-beta1);
            v_hat = v/(1-beta2);

            theta = getfield(model,name);
            theta = theta - lr * m_hat ./ (sqrt(v_hat)+eps);

            model = setfield(model,['m',name],m);
            model = setfield(model,['v',name],v);
            model = setfield(model,name,theta);
        end
    end
end

