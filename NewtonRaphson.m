classdef NewtonRaphson
    %NewtonRaphson Some parameters used for NR
    %   Detailed explanation goes here
    
    properties (SetAccess=public)
        tol;
        max_iter;
        mult;
        iter
        residual
    end
    properties (Hidden)
        step
    end
    properties (SetAccess=private, Hidden)
        numsteps
    end
    
    methods
        function obj = NewtonRaphson(NR_tol, max_iter, numsteps)
            obj.tol   = NR_tol;
            obj.max_iter = max_iter;
            obj.numsteps = numsteps;
        end
        function value = get.mult(obj)
            value = obj.step/obj.numsteps;
        end
    end
end

