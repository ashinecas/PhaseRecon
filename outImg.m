function [  ] = outImg( proj )
%UNTITLED 此处显示有关此函数的摘要
%   outImg(proj)

[~,~,num] = size(proj);
for i =1:num
    filename = strcat(num2str(i),'.tif');
    imwrite(proj(:,:,i),filename);
end

