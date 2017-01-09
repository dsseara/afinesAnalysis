clear
afile = fopen('actins.txt','r');
mfile = fopen('amotors.txt','r');
pfile = fopen('pmotors.txt','r');
i=1;

% read actin file
while feof(afile) == 0
    id = fgetl(afile);
    if (isempty(strfind(id,'='))==0)
        strtemp=id;
        strtemp(strfind(strtemp, '=')) = [];
        Key1='t'; 
        Key2='N';
        Index1 = strfind(strtemp, Key1);
        Index2 = strfind(strtemp, Key2);
        timestep(i) = sscanf(strtemp(Index1(1) + length(Key1):end), '%g', 1);
        anum(i) = sscanf(strtemp(Index2(1) + length(Key2):end), '%g', 1);
    else   
           for j = 1 : 1: anum(i)-1
               adata(j,:,i) = str2num(fgetl(afile));
           end
           i=i+1;
   end

end

fclose(afile);

%read motor file
i=1;
while feof(mfile) == 0
    id = fgetl(mfile);
    if (isempty(strfind(id,'='))==0)
        strtemp=id;
        strtemp(strfind(strtemp, '=')) = [];
        Key1='t'; 
        Key2='N';
        Index1 = strfind(strtemp, Key1);
        Index2 = strfind(strtemp, Key2);
        timestep(i) = sscanf(strtemp(Index1(1) + length(Key1):end), '%g', 1);
        mnum(i) = sscanf(strtemp(Index2(1) + length(Key2):end), '%g', 1);
    else   
           for j = 1 : 1: mnum(i)-1
               mdata(j,:,i) = str2num(fgetl(mfile));
           end
           i=i+1;
   end

end

fclose(mfile);

%read crosslink data file
i=1;
while feof(pfile) == 0
    id = fgetl(pfile);
    if (isempty(strfind(id,'='))==0)
        strtemp=id;
        strtemp(strfind(strtemp, '=')) = [];
        Key1='t'; 
        Key2='N';
        Index1 = strfind(strtemp, Key1);
        Index2 = strfind(strtemp, Key2);
        timestep(i) = sscanf(strtemp(Index1(1) + length(Key1):end), '%g', 1);
        pnum(i) = sscanf(strtemp(Index2(1) + length(Key2):end), '%g', 1);
    else   
           for j = 1 : 1: pnum(i)-1
               pdata(j,:,i) = str2num(fgetl(pfile));
           end
           i=i+1;
   end

end

fclose(pfile);

save([pwd,'\simdata.mat']); % Added by IL - 10/25/16