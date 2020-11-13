%% Partial Data

%% Parameters
nodes = [0.1 0.9 0.9 0.1;0 0 1 1];
femm_opt = struct('deg', 1, 'qdeg',4, 'min_area', 2e-5, 'edge', nodes);
opt = struct('femm_opt', femm_opt, 'reg', 1e-4, 'beta', 0.02);

% Acousto-Eletric-Modulation object.
aem_obj = aem(opt);

noise = (2 * rand(aem_obj.cache.n, 2) - 1);


mask = @(x)(1-(x(2,:) == 0));

mask_nodes = mask(aem_obj.model.space.nodes);
edge_nodes = unique(aem_obj.model.space.edges);
mask_nodes = mask_nodes(edge_nodes);

%% MAIN PROGRAM
OPTIMIZE = 1;
MAXITER = 10;
% Two data sets. [GRADIENT, MEASUREMENT]

n1 = @(x)(neumann_partial(x, mask));
n2 = @(x)(neumann2_partial(x, mask));


if OPTIMIZE
    [v1, ~, ~] = aem_obj.measurement(n1);
    
    [v2_new, v2_new_bd, ~] = selectNeumann_partial(aem_obj, v1, n2, mask_nodes);% compute the suitable v2 instead of prescribed one.
    v2_new = v2_new - mean(v2_new);                                % project to the subspace.
    [v1_new, v1_new_bd, ~] = selectNeumann_partial(aem_obj, v2_new, n1, mask_nodes); % compute the suitable v1 instead of prescribed one.
    v1_new = v1_new - mean(v1_new);        
    
    v1_new = v1_new / norm(v1_new_bd);
    v1_new_bd = v1_new_bd/norm(v1_new_bd);
    
    v2_new = v2_new / norm(v2_new_bd);
    v2_new_bd = v2_new_bd/norm(v2_new_bd);
    
    % project to the subspace.
 
    iter = 0;
    while norm(v1_new - v1)/norm(v1) > 1e-3 && iter < MAXITER
        
        iter = iter + 1;
        v1 = v1_new;
        
        [v2_new, v2_new_bd, ~] = selectNeumann_partial(aem_obj, v1, v2_new_bd, mask_nodes); % compute the suitable v2 instead of prescribed one.
        v2_new = v2_new - mean(v2_new); 

        v1_new_bd = v1_new_bd/norm(v1_new_bd);                          % project to the subspace.
        
        [v1_new, v1_new_bd, ~] = selectNeumann_partial(aem_obj, v2_new, v1_new_bd, mask_nodes); % compute the suitable v1 instead of prescribed one.
        v1_new = v1_new - mean(v1_new);                                 % project to the subspace.
        
        v1_new = v1_new / norm(v1_new_bd);
        v1_new_bd = v1_new_bd/norm(v1_new_bd);

        v2_new = v2_new / norm(v2_new_bd);
        v2_new_bd = v2_new_bd/norm(v2_new_bd);
           
    end
    
    v1g = aem_obj.gradient(v1_new);
    m1 = aem_obj.measurement_from_grad(v1g);  
    v2g = aem_obj.gradient(v2_new);
    m2 = aem_obj.measurement_from_grad(v2g);
    
else
    [v1, v1g, m1] = aem_obj.measurement(n1);
    [v2, v2g, m2] = aem_obj.measurement(n2); % prescribed.
end

% Adding noise to each of measurements.
noise_level = 0.05;

m1 = m1 .* (1 + noise(:,1) * noise_level);
m2 = m2 .* (1 + noise(:,2) * noise_level);

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

if all(abs(system_det) < 1e30)
    fprintf('Pass. The worst determinant is %6.2e.\n', max(max(system_det)));
    flag = 1;
else
    fprintf('The system fails to be stably solved.\n');
end


% figure(2);
% aem_obj.plot(aem_obj.parameter.sigma);

% Reconstruct the current J0
if flag
    J0 = aem_obj.reconstruction(v1g, v2g, m1, m2);
    
    % Check reconstruction.
    diff = J0 - aem_obj.current;
    L2error = sqrt(diff(:,1)' * aem_obj.cache.m * diff(:,1) + ...
        diff(:,2)' * aem_obj.cache.m * diff(:,2))/ ...
        sqrt(J0(:,1)' * aem_obj.cache.m *  J0(:,1) + ...
        J0(:,2)' * aem_obj.cache.m * J0(:,2));
    
    fprintf('L2 error of reconstruction is %6.2e.\n', L2error);
    
    figure('Renderer', 'painters', 'Position', [10 10 900 650]);
    subplot(2, 2, 1);
    aem_obj.plot(J0(:,1)-aem_obj.current(:,1));
    title('Error of 1st component from reconstructed J_0');
    
    subplot(2, 2, 2);
    aem_obj.plot(J0(:,2)-aem_obj.current(:,2));
    title('Error of 2nd component from reconstructed J_0');
    
    subplot(2, 2, 3);
    aem_obj.plot(aem_obj.current(:,1));
    title('1st component of exact J_0');
    
    subplot(2, 2, 4);
    aem_obj.plot(aem_obj.current(:,2));
    title('2nd component of exact J_0');
else
    fprintf('abort.\n');
end

