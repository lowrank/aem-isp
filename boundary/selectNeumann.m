function [v2, x, fval] = selectNeumann(obj, v1, neumann_bc)

n = obj.cache.n;
edge_nodes = unique(obj.model.space.edges);
if isempty(obj.cache.t)
    obj.cache.t = (obj.cache.A - 1E-8*speye(n))\ obj.cache.e(:, edge_nodes);
end

m = size(obj.cache.t, 2);
T = (obj.cache.t)' * (obj.cache.s * v1);
T = T * T';

Q = T + 1e-10* obj.cache.s(edge_nodes, edge_nodes); % add a small regularization.
f = zeros(m,1);
c = 0;

% normalization
H{1} = speye(m);
k{1} = zeros(m,1);
d{1} = -m/2;

%%
v= ver('MATLAB'); 
if strcmp(v.Release, '(R2014b)')
    
    options = optimoptions(@fmincon,'Display', 'iter', 'Algorithm','interior-point',...
        'GradObj','on','GradConstr','on',...
        'Hessian', 'bfgs',...
        'TolFun', 1e-8, 'TolX', 1e-8);
else
    
    options = optimoptions(@fmincon,'Display', 'iter','Algorithm','interior-point',...
        'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
        'HessianApproximation', 'bfgs' ,...
        'TolFun', 1e-8, 'TolX', 1e-8);

end

fun = @(x)quadobj(x,Q,f,c);
nonlconstr = @(x)quadconstr(x,H,k,d);

% initial guess from a prescribed boundary condition.

if isa(neumann_bc,'function_handle')
    x0 = neumann_bc(obj.model.space.nodes)'; % Column vector
    x0 = x0(edge_nodes);
else
    x0 = neumann_bc;
end
x0 = sqrt(m) * x0 / (x0'*x0);

[x,fval] = fmincon(fun,x0,...
    [],[],[],[],[],[],nonlconstr,options);

v2 = obj.cache.t * x;
end