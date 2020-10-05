classdef aem < handle
    %   Acousto-Electric-Modulation Class Object
    %   It is used to solve the inverse source problem for
    %   Acousto-Electric-Modulation.
    
    properties (Access = public)
        
        dimension
        cache
        model
        parameter
        reg
        current
        u_0
        grad
        load_vector
        beta
    end
    
    methods (Access = public)
        
        function obj = aem(opt)
            assert(isfield(opt, 'femm_opt'));
            assert(isfield(opt, 'beta'));
            assert(isfield(opt, 'reg'));
            
            obj.cache = struct ( ...
                'A', [], ...   % linear system (model)
                's', [], ...   % stiffness matrix (regularization use)
                'm', [], ...   % mass matrix      (regularization use)
                'n', [], ...   % number of nodes
                'dof', [], ... % degree of freedom
                'ndof', [], ...   % non degree of freedom
                'e', [], ... 
                't', []);
            
            obj.parameter = struct (...
                'sigma', []); ... % conductivity
                
            obj.beta = opt.beta;
            
            obj.model = femm(opt.femm_opt);
            
            % The cached matrices, for regularization use.
            obj.cache.s = obj.model.build('s', 1);
            obj.cache.m = obj.model.build('m', 1);
            obj.cache.e = obj.model.build('e', 1, 'all');
            
            % Nodes
            obj.cache.n = size(obj.model.space.nodes, 2);
            obj.cache.ndof = unique(obj.model.space.edges);
            obj.cache.dof = setdiff(1:obj.cache.n, obj.cache.ndof);
            
            
            % Parameter (user supplied)
            obj.parameter.sigma = conductivity(obj.model.space.nodes)';
            
            % Regularization in case of noise
            obj.reg = opt.reg;
            
            % create the exact current J_0 (user supplied, vector valued)
            obj.current = current_func(obj.model.space.nodes)'; % note here transpose
            
            % compute the div J_0
            div_current = obj.divergence(obj.current);
                
            
            % map the conductivity to quadrature points.
            q_sigma = obj.mapping(obj.parameter.sigma, ...
                                  obj.model.space.elems, ...
                                  obj.model.facet.ref');
                              
            q_div_current = obj.mapping(div_current, ...
                                  obj.model.space.elems, ...
                                  obj.model.facet.ref');
                          
                              
            obj.cache.A = -obj.model.build('s', q_sigma); % matrix A, but singular.
            
            obj.load_vector = obj.model.build('l', q_div_current); % vector b
            
            % Solve A u = b by augmenting one more row to A and one more
            % row to b.
            
            A = obj.cache.A; %A(end,:) = 1; % still sparse format.
            b = obj.load_vector; %b(end) = 0;
           
            % Or use iterative methods
 
            [obj.u_0, ~] = bicgstab(A, b - mean(b), 1e-10, 1000);
            obj.u_0 = obj.u_0 - mean(obj.u_0);

            
            obj.grad = obj.gradient(obj.u_0);
            
           
           
      
        end
        
        function [v, v_grad] = forward_solve(obj, neumann_func)

            q = obj.model.builder.qnodes1D(...
                obj.model.space.nodes, ...
                obj.model.edge.qnodes, ...
                obj.model.space.edges);
            
            load = neumann_func(q);

            b = obj.model.build('g', load, 'all');
            
            A = obj.cache.A; 
                 
            [v, ~]  = bicgstab(A, b - mean(b), 1e-10, 1000);
            v = v - mean(v);
            
            v_grad = obj.gradient(v);

            
        end
        
        
        function [v, v_grad, m]= measurement(obj, neumann_func)
            [v, v_grad] = obj.forward_solve(neumann_func);
            
            m = (obj.beta * obj.parameter.sigma .* obj.grad(:,1) - obj.current(:,1)) .* ...
                v_grad(:, 1) + ...
                (obj.beta * obj.parameter.sigma .* obj.grad(:,2) - obj.current(:,2)) .* ...
                v_grad(:, 2);
        end
        
        function [m]= measurement_from_grad(obj, v_grad)
            
            m = (obj.beta * obj.parameter.sigma .* obj.grad(:,1) - obj.current(:,1)) .* ...
                v_grad(:, 1) + ...
                (obj.beta * obj.parameter.sigma .* obj.grad(:,2) - obj.current(:,2)) .* ...
                v_grad(:, 2);
        end
        
        function J0 = reconstruction(obj, vg1, vg2, m1, m2)
            
            K = [spdiags(vg1(:,1), 0, obj.cache.n, obj.cache.n) ...
                spdiags(vg1(:,2), 0, obj.cache.n, obj.cache.n);...
                spdiags(vg2(:,1), 0, obj.cache.n, obj.cache.n)...
                spdiags(vg2(:,2), 0, obj.cache.n,obj.cache.n)];
            
            

            z =  [m1;m2];
            % make sure K is non-singular.
            H = K \ z;
            
            H = reshape(H, obj.cache.n, 2);
            
            % make sure beta is not 1.
            divJ0 = obj.divergence(H) / (obj.beta - 1);
            
            
            q_divJ0 = obj.mapping(divJ0, ...
                                  obj.model.space.elems, ...
                                  obj.model.facet.ref');
            b = obj.model.build('l', q_divJ0);
            
            
            [U0, ~] = bicgstab(obj.cache.A, b - mean(b), 1e-10, 1000);
            U0 = U0 - mean(U0);

            
            gradU0 = obj.gradient(U0);
            
            J0 = [obj.beta * obj.parameter.sigma .* gradU0(:,1) - H(:,1) ...
                obj.beta * obj.parameter.sigma .* gradU0(:,2) - H(:,2)];
            
            % Impose zero boundary constraint on J0
            J0(unique(obj.model.space.edges),:) = 0;
            
            
        end
        
        
        function val = divergence(obj, u)
            % Input must be N x 2
            assert(size(u, 1) == obj.cache.n);
            assert(size(u, 2) == 2);

            [DX, DY] = obj.model.builder.reference_grad(obj.model.rffs.nodes);
            [val_x, ~] = obj.model.gradient(u(:,1), DX, DY);
            [~, val_y] = obj.model.gradient(u(:,2), DX, DY);

            val = zeros(obj.cache.n, 1);
            for i = 1:size(obj.model.space.elems, 2)
                val(obj.model.space.elems(:, i), 1) = val_x(:, i) + val_y(:, i);
            end
        end      


        function [grad] = gradient(obj, u)
            [DX, DY] = obj.model.builder.reference_grad(obj.model.rffs.nodes);
            [val_x, val_y] = obj.model.gradient(u, DX, DY);
            
            grad = zeros(obj.cache.n, 2);
            for i = 1:size(obj.model.space.elems, 2)
                grad(obj.model.space.elems(:, i), 1) = val_x(:, i);
                grad(obj.model.space.elems(:, i), 2) = val_y(:, i);
            end
            
        end
        
        
        function plot(obj, f)
            trisurf(obj.model.space.elems(1:3,:)', ...
                obj.model.space.nodes(1,:), obj.model.space.nodes(2,:), ...
                f);colormap jet;view(2);colorbar;shading interp;
        end
    end
    
    methods(Static)

        function [interpolate] = mapping(func, elems, trans_ref)
            numberofqnodes = size(trans_ref, 1);
            interpolate = zeros(numberofqnodes, size(elems, 2));
            for i = 1: size(elems, 2)
                interpolate(:, i) = trans_ref * func(elems(:, i));
            end
        end
        
        
        function r = normsq(v)
            r = sum(v.^2);
        end


    end
    
end

