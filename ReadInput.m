function scriptStr = ReadInput()
  inputPath = [pwd filesep 'Input' filesep];

  list = dir(inputPath);

  indx = 0;
  names = cell(size(list,1),1);
  for i = 1:size(list,1)
    if size(list(i).name,2) >1
      if strcmp(list(i).name(end-1:end),'.m')
        indx = indx+1;
        names{indx} = list(i).name;
      end
    end
  end

  height = 15*indx*(15*indx<1000) + 1000*(15*indx>=1000);
  [Selection,ok] = listdlg(...
    'ListString'   , names(1:indx),...
    'SelectionMode', 'single',...
    'ListSize'     , [300 height],...
    'Name'         , 'Input Selection',...
    'PromptString' , 'Select one input file:',...
    'OKString'     , 'Yay',...
    'CancelString' , 'Nay');

  if ok
    scriptStr = [inputPath names{Selection}];
    fid = fopen([inputPath 'recent.txt'],'w');
    fprintf(fid, scriptStr);
    fclose(fid);
  else
    try
      fid = fopen(recent,'r');
      scriptStr = fscanf(fid,'%s');
      fclose(fid);
    catch
      error(" You haven't selected a file and there's no recent file saved")
    end
  end
             
end
