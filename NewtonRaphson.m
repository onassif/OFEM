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
        time=0;
        fctr=0;
    end
    
    methods
        function obj       = NewtonRaphson(NR_tol, max_iter, num, varargin)
            obj.tol        = NR_tol;
            obj.max_iter   = max_iter;
            obj.numsteps   = num.steps;
            mult_length    = (num.BC>=num.FORCE)*num.BC + (num.FORCE>num.BC)*num.FORCE;
            obj.mult       = zeros(mult_length,2);
            switch nargin
                case 4
                    obj.time     =[zeros(1,size(varargin{1},2)); varargin{1}];
                    obj.time(:,2)=cumsum(obj.time(:,2));
                case 5
                    obj.time     =[zeros(1,size(varargin{1},2)); varargin{1}];
                    obj.time(:,2)=cumsum(obj.time(:,2));
                    obj.fctr=[zeros(mult_length,1) varargin{2}(:,:,1), zeros(mult_length,1), varargin{2}(:,:,2)];
                    obj.fctr=reshape(obj.fctr, mult_length , size(varargin{2},2)+1,2);
            end
        end
        function obj = set.mult(obj,value)
            if (~obj.time | length(obj.fctr)==1)%if user didn't pass time or factors or Initialization
                obj.mult(:,:) = value;
            else
                row = find(obj.step>obj.time(:,2),1,'last')+1;
                ncrmnt = squeeze(...
                    (obj.fctr(:,row,:)-obj.fctr(:,row-1,:))/...
                    (obj.time(  row,2)-obj.time(  row-1,2)) );
                
                obj.mult(:,:) = obj.mult(:,:) + ncrmnt;
            end
        end
    end
end
