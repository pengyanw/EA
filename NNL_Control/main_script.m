% Lyapunov NN: V_\theta(x)
layersV = [
    featureInputLayer(2)
    fullyConnectedLayer(32)
    tanhLayer
    fullyConnectedLayer(32)
    tanhLayer
    fullyConnectedLayer(1)
];

dlnetV = dlnetwork(layersV);

% Controller NN: u_\phi(x)
layersU = [
    featureInputLayer(2)
    fullyConnectedLayer(32)
    tanhLayer
    fullyConnectedLayer(32)
    tanhLayer
    fullyConnectedLayer(1)
];
%% 
numEpoch = 2000;
lr = 1e-3;
mb = 64;

trA_V = []; trB_V = [];
trA_U = []; trB_U = [];

for epoch = 1:numEpoch
    
    % Batch
    X = randn(2, mb);
    dlX = dlarray(X,'CB');

    % Compute loss inside dlfeval
    loss = dlfeval(@lyapunov_loss, dlnetV, dlnetU, dlX);

    % === Combine learnables ===
    allLearn = [dlnetV.Learnables; dlnetU.Learnables];

    % === Compute gradients wrt all learnables ===
    grads = dlgradient(loss, allLearn);

    % === Split gradients back ===
    nV = numel(dlnetV.Learnables);
    gradV = grads(1:nV);
    gradU = grads(nV+1:end);

    % === Update both networks ===
    [dlnetV, trA_V, trB_V] = adamupdate(dlnetV, gradV, trA_V, trB_V, epoch, lr);
    [dlnetU, trA_U, trB_U] = adamupdate(dlnetU, gradU, trA_U, trB_U, epoch, lr);

    if mod(epoch,100)==0
        fprintf("Epoch %d: Loss = %.5f\n", epoch, extractdata(loss));
    end
end

