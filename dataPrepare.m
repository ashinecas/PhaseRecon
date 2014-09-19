function  dataPrepare(projFilePath)
%dataPrepare convert .dat format data to .tif format data. Output img in a
%new folder.
%projFilePath:  Path of proj data  
%Syntax:
%dataPrepare('F:\data\source Data')

path(path,projFilePath); % add projFilePath to the end of search path.
AbsoPath = strcat(projFilePath,'\*.dat');
ImgFiles = dir(AbsoPath);
numfiles = length(ImgFiles);


%% Test path.
try
projSingle=readDat(ImgFiles(1).name);
   %ImgeInfo= imfinfo(ImgFiles(1).name);  %remember to add the source data to your path, 
                                        %or you cannot access those files.
catch exception1
    error('Error: cannot Read original image file, Please Check your path again. ');
end
%% Create a new file folder to save converted img.
newpath=fullfile(pwd,'data',date);
if ~isdir(newpath)
    mkdir(newpath);
    path(path,newpath);
end
%%
for i = 1:numfiles
    projSingle = readDat(ImgFiles(i).name);
    projSingle = projSingle(1000:1010,:);
    filename = fullfile(newpath,strcat(num2str(i),'.tif'));
    % filename = strcat(num2str(i),'.tif');
    imwrite(projSingle,filename);
end


% %default setting is double percision, use uint16 to save space.
% proj = zeros(ImgeInfo.Width,ImgeInfo.Height,numfiles,'uint16'); 
% 
% 
% for k = 1:numfiles
%     %myfilename = sprintf('%d.tif',k);
%     myfilename = ImgeFiles(k).name;
%     proj(:,:,k) = imread(myfilename);
% end
% 
% mypath=fullfile(pwd,date,'data');

% 
% save(fullfile(mypath,'proj.mat'),'proj');
% 
% 
% end
% 
% 
% 
% proj_down = projDownSample(proj,0.25);
% 
% outImgpath = fullfile(pwd,date,'postproc_proj')
% mkdir(outImgpath);
% save(fullfile(outImgpath,'proj_down.mat'),'proj_down');
% cd(outImgpath);
% outImg( proj_down );
