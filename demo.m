function demo()

%% Parameters
nodes = [0 1 1 0;0 0 1 1];
femm_opt = struct('deg', 1, 'qdeg',4, 'min_area', 1e-4, 'edge', nodes);
opt = struct('femm_opt', femm_opt, 'reg', 1e-4, 'beta', 0.02);

OPTIMIZE =0;

% Acousto-Eletric-Modulation object.
aem_obj = aem(opt);

%% Two data sets. [GRADIENT, MEASUREMENT]
[v1, v1g, m1]=aem_obj.measurement(@neumann);

if OPTIMIZE
    [v2, ~, ~] = selectNeumann(aem_obj, v1); % compute the suitable v2 instead of prescribed one.
    v2 = v2 - mean(v2);              % project to the subspace.
    v2g = aem_obj.gradient(v2);
    m2 = aem_obj.measurement_from_grad(v2g);
else
    [v2, v2g, m2] = aem_obj.measurement(@neumann2); % prescribed.
end

% Adding noise to each of measurements.
noise_level = 0.05;

m1 = m1 .* (1 + (2 * rand(size(m1)) - 1) * noise_level);
m2 = m2 .* (1 + (2 * rand(size(m2)) - 1) * noise_level);

% Determinant of linear system (Check to avoid blowing system error)

% Should optimize the min(abs(det X)), X = [grad v1 grad v2] 
% with grad v1, grad v2 normalized. Basically means the two vectors must be
% having a relatively large angle in between.
theta  = (v1g(:,2) .* v2g(:,1) - v2g(:,2) .* v1g(:, 1));
system_det = [ abs(   v2g(:, 1)  ./ theta )   ...
               abs(   v2g(:, 2)  ./ theta )   ...
               abs(   v1g(:, 1)  ./ theta )   ...
               abs(   v1g(:, 2)  ./ theta )   ...
               ];

flag = 0;

if all(abs(system_det) < 1e12)
    fprintf('Pass. The worst determinant is %6.2e.\n', max(max(system_det)));
    flag = 1;
else
    fprintf('The system fails to be stably solved.\n');
end

% Reconstruct the current J0
if flag
    J0 = aem_obj.reconstruction(v1g, v2g, m1, m2);
    
    % Check reconstruction.
    diff = J0 - aem_obj.current;
    L2error = sqrt(diff(:,1)' * aem_obj.cache.m * diff(:,1) + ...
        diff(:,2)' * aem_obj.cache.m * diff(:,2));
    
    fprintf('L2 error of reconstruction is %6.2e.\n', L2error);
    
    figure('Renderer', 'painters', 'Position', [10 10 900 650]);
    subplot(2, 2, 1);
    aem_obj.plot(J0(:,1));
    title('1st component of reconstructed J_0');
    
    subplot(2, 2, 2);
    aem_obj.plot(J0(:,2));
    title('2nd component of reconstructed J_0');
    
    subplot(2, 2, 3);
    aem_obj.plot(aem_obj.current(:,1));
    title('1st component of exact J_0');
    
    subplot(2, 2, 4);
    aem_obj.plot(aem_obj.current(:,2));
    title('2nd component of exact J_0');
else
    fprintf('abort.\n');
end

end