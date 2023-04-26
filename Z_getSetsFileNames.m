function [fileNames, pathName]=Z_getSetsFileNames(ext)

if nargin < 1
    [fileNames,pathName] = uigetfile({'*.*'},'Choose Files','multiselect','on');
else
    [fileNames,pathName] = uigetfile({['*.' ext]},['Choose ' ext ' Files'],'multiselect','on');
end

if ~iscell(fileNames)
    fileNames={fileNames};
end

if size(fileNames, 2)>1
    fileNames=fileNames';
end

end