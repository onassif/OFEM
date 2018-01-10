function scriptStr = ReadInput()
cdir   = sprintf('%s%sInput',pwd,filesep);
recent = sprintf('%s%sInput%srecent.txt', pwd, filesep, filesep);

list = dir(cdir);

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

height = 20*indx*(20*indx<1000) + 1000*(20*indx>=1000);
[Selection,ok] = listdlg('ListString',names(1:indx), 'SelectionMode','single',...
    'ListSize',[300 height], 'Name', 'Input Selection',...
    'PromptString', 'Select one input file:', 'OKString', 'Yay', 'CancelString', 'Nay');

if ok
    scriptStr = names{Selection};
    fid = fopen(recent,'w');
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
