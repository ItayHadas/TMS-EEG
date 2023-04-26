function [fileNames, pathName]=Z_getSetsFileNames(ext,pathi)

if nargin < 1
    [fileNames,pathName] = uigetfile({['*.*']},'Choose Files','multiselect','on');
elseif nargin == 1
    [fileNames,pathName] = uigetfile({['*.' ext]},'Choose Files','multiselect','on');
elseif nargin == 2
    [fileNames,pathName] = uigetfile({[pathi '*.' ext]},['Choose ' ext ' Files'],'multiselect','on');
end

if ~iscell(fileNames)
    fileNames={fileNames};
end

if size(fileNames, 2)>1
    fileNames=fileNames';
end

end