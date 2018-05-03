classdef FCC
   %UNTITLED4 Summary of this class goes here
   %   Detailed explanation goes here
   
   properties
      nslip = 12;
      
      b;
      n;
   end
   
   methods
      function obj = FCC()
         f = 1/sqrt(2);
         obj.b = [...
            .0 -f  f
            +f  0 -f
            -f  f  0
            .0  f  f
            +f  0  f
            +f -f  0
            .0 -f  f
            -f  0 -f
            +f  f  0
            .0  f  f
            +f  0 -f
            -f -f  0];
         
         g = 1/sqrt(3);
         obj.n = [...
            +g  g  g
            +g  g  g
            +g  g  g
            -g -g  g
            -g -g  g
            -g -g  g
            -g  g  g
            -g  g  g
            -g  g  g
            +g -g  g
            +g -g  g
            +g -g  g];
      end
   end
end

