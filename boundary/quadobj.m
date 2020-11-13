function [y,grady] = quadobj(x,Q,f,c, mask)
y = 1/2*x'*Q*x + f'*x + c;

if nargin < 5
    mask = ones(size(x))';
end

if nargout > 1
    grady = Q*x + f;
    grady = grady .* mask';
end


