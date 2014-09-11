%Read .Tiff files and save them into a single file.
function [proj,mypath] = readImg(projFilePath,filetype)
%% projFilePath:Path of data  
%% filetpe: Image type
%% proj: 3D volume  Width*Height*projnum   uint16
%% [proj,datapath] = readImg('F:\data\source Data','tif')

path(path,projFilePath);
AbsoPath = strcat(projFilePath,'\*.',filetype);
%cd(projFilePath);
%Imgetype= strcat('*.',filetype);
%ImgeFiles = dir(Imgetype);
ImgeFiles = dir(AbsoPath);
numfiles = length(ImgeFiles);
%numfiles = 10; %test

try
ImgeInfo= imfinfo(ImgeFiles(1).name);  %remember to add the source data to your path, 
                                        %or you cannot access those files.
catch exception1
    error('Error: Remember to add source data to your path');
end

%default setting is double percision, use uint16 to save space.
proj = zeros(ImgeInfo.Width,ImgeInfo.Height,numfiles,'uint16'); 


for k = 1:numfiles
    %myfilename = sprintf('%d.tif',k);
    myfilename = ImgeFiles(k).name;
    proj(:,:,k) = imread(myfilename);
end

mypath=fullfile(pwd,date,'data');
if ~isdir(mypath)
    mkdir(mypath);
    path(path,mypath);
end

save(fullfile(mypath,'proj.mat'),'proj');


end


