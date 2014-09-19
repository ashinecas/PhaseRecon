function projSingle = readDat( dataDectory )
%readDat Read single projection from .dat format data. 
%   Syntax: proj = readDat( dataDectory )
fileID = fopen(dataDectory);
projSingle = fread(fileID,[2154,1326],'double');
%proj = proj(850:2154,:);
end

