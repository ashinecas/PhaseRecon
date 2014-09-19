% Agray=mat2gray(A);
% grayI=uint8(zeros(231,231,231));
% for i=1:188
%     grayI(:,:,i)=uint8(round(Agray(:,:,i)*255));
% end
grayIre= grayI(:,50:180,:);
fid=fopen('mouse0718-3.raw','w+');
cnt=fwrite(fid,grayIre,'uint8');
fclose(fid);