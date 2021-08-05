classdef MatrixSpline < handle
    properties
        a {mustBeNumeric}
        b
        num_segments
        output_dim
        segment_length
        num_knots
        num_basis
        knot_sequence
        coefficients
    end
    methods
        function obj = MatrixSpline()
            
        end
        
        function setup(obj, a, b, num_segments, output_dim)
            obj.a = a;
            obj.b = b;
            obj.num_segments = num_segments;
            obj.output_dim = output_dim;
            
            obj.segment_length = (b - a)/num_segments;
        
            obj.num_knots = num_segments + 7;
            obj.num_basis = num_segments + 3;
            
            obj.knot_sequence = (a - 3*obj.segment_length):obj.segment_length:(b + 3*obj.segment_length);
            obj.coefficients = zeros(output_dim*obj.num_basis,1);
        end

        function set_coefficients(obj, coefficients)
            if length(obj.coefficients) == length(coefficients)
                obj.coefficients = coefficients;
            else
                error('Incorrect length of vector coefficients');
            end
        end
        
        function x = evaluate(obj, t, deriv_order)
            x = zeros(obj.output_dim, length(t));

            for iii = 1:length(t)
                if (t(iii) <= obj.a) || (t(iii) >= obj.b)
                    x(:,iii) = nan(obj.output_dim, 1);
                else
                    phi = get_phi(obj, t(iii), deriv_order);
                    x(:,iii) = phi*obj.coefficients;
                end
            end
        end
        
        function phi = get_phi(obj, t, deriv_order)
            persistent B D
            
            B = 1/6*[1  4  1  0
                    -3  0  3  0
                     3 -6  3  0
                    -1  3 -3  1];
    
            D =     [0  1  0  0
                     0  0  2  0
                     0  0  0  3
                     0  0  0  0];
    
            phi = zeros(obj.output_dim, obj.num_basis*obj.output_dim);
            interval = find_interval(obj, t);
    
            t_ii = obj.knot_sequence(interval);
            t_ip1 = obj.knot_sequence(interval + 1);
            u = (t - t_ii)/(t_ip1 - t_ii);
            U = [1, u, u^2, u^3];
    
            for iii = 1:obj.output_dim
                b1 = obj.output_dim*(interval - 4) + iii;
                b2 = obj.output_dim*(interval - 3) + iii;
                b3 = obj.output_dim*(interval - 2) + iii;
                b4 = obj.output_dim*(interval - 1) + iii;

                P_i = zeros(4, obj.num_basis*obj.output_dim);
                P_i(1, b1) = 1;
                P_i(2, b2) = 1;
                P_i(3, b3) = 1;
                P_i(4, b4) = 1;

                phi(iii,:) = U*(D^deriv_order)*B*P_i;
            end
            
            phi = phi/(obj.segment_length^deriv_order);
        end
        
        function interval = find_interval(obj, t)
            if (t < obj.a) || (t > obj.b)
                warning('t not within proper interval');
            end
            
            interval = find(diff(t >= obj.knot_sequence));
        end
        
        function fit(obj, t, x, data_per_coeff)
            if nargin == 3
                data_per_coeff = 1;
            end
            
            obj.num_basis = floor(length(t)/data_per_coeff);
            obj.num_segments = obj.num_basis - 3;
            obj.setup(t(1), t(end) + 1e-6, obj.num_segments, size(x,1));
            
            x = x(:);
            
            % Build stacked phi matrix
            stacked_phi = nan(length(t)*obj.output_dim, obj.num_basis*obj.output_dim);
            
            for iii = 1:length(t)
                phi = get_phi(obj, t(iii), 0);
                lo = (iii - 1)*obj.output_dim + 1;
                hi = lo + obj.output_dim - 1;
                stacked_phi(lo:hi,:) = phi;
            end
            
            % Now we have x = stacked_phi*c. Let's solve for the coeff.
            c = stacked_phi \ x;
            
            obj.set_coefficients(c);
        end
    end
end