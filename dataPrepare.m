
[proj,projdatapath] = readImg('F:\data\mouse0705','tif');
proj_down = projDownSample(proj,0.25);

outImgpath = fullfile(pwd,date,'postproc_proj')
mkdir(outImgpath);
save(fullfile(outImgpath,'proj_down.mat'),'proj_down');
cd(outImgpath);
outImg( proj_down );
