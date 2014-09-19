function [ proj_new ] = projDownSample( proj , downsample, method)
%projDownsample 
%   Down sampling projections, in order to save memory and time for test.
%   Method = inrerpolation method;
%   'nearest','bilinear','bicubic',... refer to imresize() manual.
%   proj_down = projDownSample(proj,0.25)
if nargin == 1;
    method = 'bicubic';
    downsample = 1.0;
end 

if nargin == 2;
    method = 'bicubic';
end 

[width,height,num] = size(proj);
proj_new = zeros(width*downsample,height*downsample,num,'uint16');

for i = 1: num
   proj_new(:,:,i) = imresize(proj(:,:,i),downsample,method);
end


end



