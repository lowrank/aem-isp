function [v2, x, fval] = selectNeumann(obj, v1)

n = obj.cache.n;
E = obj.model.build('e', 1, 'all');

T = (obj.cache.s * v1);
T = T * T';

edge_nodes = unique(obj.model.space.edges);
m = length(edge_nodes);

M = (obj.cache.A - 1e-6*speye(n))\ E(:, edge_nodes);

Q = M' * T * M + 0 * obj.cache.s(edge_nodes, edge_nodes); % add a small regularization.
f = zeros(m,1);
c = 0;

% normalization
H{1} = speye(m);
k{1} = zeros(m,1);
d{1} = -m/2;

%%
options = optimoptions(@fmincon,'Display', 'iter', 'Algorithm','interior-point',...
    'GradObj','on','GradConstr','on',...
    'HessMult',@(x,lambda)quadhess(x,lambda,Q,H),...
    'TolFun', 1e-6, 'TolX', 1e-8);

fun = @(x)quadobj(x,Q,f,c);
nonlconstr = @(x)quadconstr(x,H,k,d);

% initial guess from a prescribed boundary condition.
x0 = neumann2(obj.model.space.nodes)'; % Column vector
x0 = x0(edge_nodes);
x0 = sqrt(m) * x0 / (x0'*x0);

[x,fval] = fmincon(fun,x0,...
    [],[],[],[],[],[],nonlconstr,options);


v2 = M * x;
end