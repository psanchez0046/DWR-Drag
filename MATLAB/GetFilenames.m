function [ names_cell ] = GetFilenames(pattern, inputDataFilepath)
% This function returns a cell array with the names of the files 
% having the specific pattern contained in the specified filepath

dinfo = dir(inputDataFilepath);
names_cell = {dinfo.name};
aux=0;
names_cell_aux = names_cell;
for i=1:length(names_cell)
    a=strfind(names_cell(i),pattern);
    if a{1,1}~=0
        aux=aux+1;
        names_cell_aux(aux)=names_cell(i);
    end
end    
names_cell=names_cell_aux(1:aux);
end