% Parameters
femm_opt = struct('deg', 2, 'qdeg',4, 'min_area', 1e-4, 'edge', nodes);
opt = struct('femm_opt', femm_opt, 'reg', 1e-4, 'beta', 0.02);

% Acousto-Eletric-Modulation object.
aem_obj = aem(opt);

%% Two data sets. [GRADIENT, MEASUREMENT]
[vg1, m1]=aem_obj.measurement(@neumann);
[vg2, m2]=aem_obj.measurement(@neumann2);

% Adding noise to each of measurements.
noise_level = 0.02;

m1 = m1 .* (1 + (2 * rand(size(m1)) - 1) * noise_level);
m2 = m2 .* (1 + (2 * rand(size(m2)) - 1) * noise_level);

% Determinant of linear system (Check to avoid blowing system error)

% Should optimize the min(abs(det X)), X = [grad v1 grad v2] 
% with grad v1, grad v2 normalized. Basically means the two vectors must be
% having a relatively large angle in between.
system_det = vg1(:,2) ./ vg1(:,1) - vg2(:,2) ./ vg2(:, 1);

flag = 0;

if all(abs(system_det) > 1e-12)
    fprintf('Pass. The worst determinant is %6.2e.\n', min(abs(system_det)));
    flag = 1;
else
    fprintf('The system fails to be stably solved.\n');
end

% Reconstruct the current J0
if flag
    J0 = aem_obj.reconstruction(vg1, vg2, m1, m2);
    
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