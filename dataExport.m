fid = fopen('testfile','w+');
fwrite(fid,uint32(size(Z,1)),'uint32');
fwrite(fid,uint32(size(Z,2)),'uint32');
fwrite(fid,Z,'double');
fclose(fid);

system('mv testfile ~/Documents/Development/Python/testfile');

vert = vertexEnumeration(expLambda,explambda);
rf = convhull(vert(:,1),vert(:,2));

fid = fopen('testfile1','w+');
fwrite(fid,uint32(size(vert,2)),'uint32');
fwrite(fid,uint32(length(rf)),'uint32');
fwrite(fid,vert(rf,:)','double');
fclose(fid);

system('mv testfile1 ~/Documents/Development/Python/testfile1');


vert = vertexEnumeration(sLambda,slambda);
rf = convhull(vert(:,1),vert(:,2));

fid = fopen('testfile2','w+');
fwrite(fid,uint32(size(vert,2)),'uint32');
fwrite(fid,uint32(length(rf)),'uint32');
fwrite(fid,vert(rf,:)','double');
fclose(fid);

system('mv testfile2 ~/Documents/Development/Python/testfile2');

vert = vertexEnumeration(eLambda,elambda);
rf = convhull(vert(:,1),vert(:,2));

fid = fopen('testfile3','w+');
fwrite(fid,uint32(size(vert,2)),'uint32');
fwrite(fid,uint32(length(rf)),'uint32');
fwrite(fid,vert(rf,:)','double');
fclose(fid);

system('mv testfile3 ~/Documents/Development/Python/testfile3');


vert = vertexEnumeration(iLambda,ilambda);
rf = convhull(vert(:,1),vert(:,2));

fid = fopen('testfile4','w+');
fwrite(fid,uint32(size(vert,2)),'uint32');
fwrite(fid,uint32(length(rf)),'uint32');
fwrite(fid,vert(rf,:)','double');
fclose(fid);

system('mv testfile4 ~/Documents/Development/Python/testfile4');