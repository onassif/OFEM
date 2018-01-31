classdef ConstantMixedHardening
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        iso_modulus
        kin_modulus
        mixedHardParam
        initYieldStress
    end
    
    methods
        function obj = ConstantMixedHardening(hardProps)

            for i=1:length(hardProps)
                switch hardProps{i,1}
                    case 'K'
                        obj.iso_modulus    = hardProps{i,2};
                    case 'H'
                        obj.kin_modulus    = hardProps{i,2};
                    case 'M'
                        obj.mixedHardParam = hardProps{i,2};
                    case 'sig0'
                        obj.initYieldStress= hardProps{i,2};
                    otherwise
                        error("You've chosen Plane Strain material but specified incompatible material properties, I'm disapponted");
                end
            end
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

