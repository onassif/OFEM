options ={...
'Tensile  5 steps  x-dir', 'Tensile 10 steps x-dir','Tensile 20 steps x-dir',...
'Tensile   5 steps y-dir', 'Tensile 10 steps y-dir','Tensile 20 steps y-dir',...
'Shear     5 steps'      , 'Shear   10 steps'      , 'Shear   20 steps', ...
'Tensile, 2 nodes fixed,  5 steps  x-dir (for debugging)'};

[Selection,ok] = listdlg('PromptString','Select A problem to solve',...
    'Name', 'Problem Selection                    ','SelectionMode',...
    'single','ListString',options,'ListSize',[265 132]);

if ok
    inputParameters
     if     (Selection ==1 || Selection ==2 || Selection ==3)        
        BC=BC_T; FORCE=FORCE_Tx;
            
        if     (Selection == 1)     n_steps = 5;
        elseif (Selection == 2)     n_steps = 10;
        elseif (Selection == 3)     n_steps = 20;
        end 
        
    elseif  (Selection ==4 || Selection ==5 || Selection ==6)
        BC=BC_T; FORCE=FORCE_Ty;
        
        if     (Selection == 4)     n_steps = 5;
        elseif (Selection == 5)     n_steps = 10;
        elseif (Selection == 6)     n_steps = 20;
        end
        
    elseif  (Selection ==7 || Selection ==8 || Selection ==9)
        BC=BC_sh; FORCE=FORCE_sh;
        
        if     (Selection == 7)     n_steps = 5;
        elseif (Selection == 8)     n_steps = 10;
        elseif (Selection == 9)     n_steps = 20;
        end
        
     elseif (Selection ==10)
         BC=BC_Tp; FORCE=FORCE_Tp;   n_steps = 5;
    end
else
    error('>>>>> Error!! No option selected, exiting...');
end

clear('BC_sh', 'BC_T', 'BC_Tp', 'FORCE_sh', 'FORCE_x', 'FORCE_y',...
 'FORCE_Tp', 'Selection', 'ok', 'options')