classdef MTSob
   %UNTITLED Summary of this class goes here
   %   Detailed explanation goes here
   
   properties
      tau_a
      tau_ht_y
      g0_y
      q_y
      p_y
      eps0_dt_y
      tau_ht_v
      g0_v
      q_v
      p_v
      eps0_dt_v
      k0
      boltz
      b
      theta0
      mu0
      D0
      t0
      hardenExp
      E
      nu
      voche_m
      elasticType
      numCrystals
      angleConvention
      angleType
      angles
      
      slip
   end
   
   methods
      function obj = MTS(hardProps, slip)
         
         for i=1:length(hardProps)
            switch hardProps{i,1}
               case 'tau_a'
                  obj.tau_a            = hardProps{i,2};
               case 'tau_ht_y'
                  obj.tau_ht_y         = hardProps{i,2};
               case 'g0_y'
                  obj.g0_y             = hardProps{i,2};
               case 'q_y'
                  obj.q_y              = hardProps{i,2};
               case 'p_y'
                  obj.p_y              = hardProps{i,2};
               case 'eps0_dt_y'
                  obj.eps0_dt_y        = hardProps{i,2};
               case 'tau_ht_v'
                  obj.tau_ht_v         = hardProps{i,2};
               case 'g0_v'
                  obj.g0_v             = hardProps{i,2};
               case 'q_v'
                  obj.q_v              = hardProps{i,2};
               case 'p_v'
                  obj.p_v              = hardProps{i,2};
               case 'eps0_dt_v'
                  obj.eps0_dt_v        = hardProps{i,2};
               case 'k0'
                  obj.k0               = hardProps{i,2};
               case 'boltz'
                  obj.boltz            = hardProps{i,2};
               case 'b'
                  obj.b                = hardProps{i,2};
               case 'theta0'
                  obj.theta0           = hardProps{i,2};
               case 'mu0'
                  obj.mu0              = hardProps{i,2};
               case 'D0'
                  obj.D0               = hardProps{i,2};
               case 't0'
                  obj.t0               = hardProps{i,2};
               case 'hardenExp'
                  obj.hardenExp        = hardProps{i,2};
               case 'E'
                  obj.E                = hardProps{i,2};
               case 'nu'
                  obj.nu               = hardProps{i,2};
               case 'voche_m'
                  obj.voche_m          = hardProps{i,2};
               case 'elasticType'
                  obj.elasticType      = hardProps{i,2};
               case 'numCrystals'
                  obj.numCrystals      = hardProps{i,2};
               case 'angleConvention'
                  obj.angleConvention  = hardProps{i,2};
               case 'angleType'
                  obj.angleType        = hardProps{i,2};
               case 'angles'
                  obj.angles           = hardProps{i,2};
               otherwise
                  error("Unknown MTS hardening parameter");
            end
            
            obj.slip = slip;
         end
      end
   end
end