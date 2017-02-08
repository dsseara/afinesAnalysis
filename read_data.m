clear
if ispc %change delimiter based on machine being used
    afile = fopen([pwd,'\txt_stack\actins.txt'], 'r');
    mfile = fopen([pwd,'\txt_stack\amotors.txt'],'r');
    pfile = fopen([pwd,'\txt_stack\pmotors.txt'],'r');
elseif isunix
    afile = fopen([pwd,'/txt_stack/actins.txt'], 'r');
    mfile = fopen([pwd,'/txt_stack/amotors.txt'],'r');
    pfile = fopen([pwd,'/txt_stack/pmotors.txt'],'r');
end

% read actin file
i=1;
ii = 1; % Separate counter for timestep
% Guarantees proper incrementing for simulations with zero motor or
% crosslinker densities (a_motor_density = 0 OR p_motor_density = 0)-
% (IAL:1/31/17)
while feof(afile) == 0
    id = fgetl(afile);
    if (isempty(strfind(id,'='))==0)
        strtemp=id;
        strtemp(strfind(strtemp, '=')) = [];
        Key1='t'; 
        Key2='N';
        Index1 = strfind(strtemp, Key1);
        Index2 = strfind(strtemp, Key2);
        timestep(ii) = sscanf(strtemp(Index1(1) + length(Key1):end), '%g', 1);
        anum(i) = sscanf(strtemp(Index2(1) + length(Key2):end), '%g', 1);
        
        ii = ii+1; % Increment timestep counter every loop  
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
ii = 1; % Separate counter for timestep
% Guarantees proper incrementing for simulations with zero motor or
% crosslinker densities (a_motor_density = 0 OR p_motor_density = 0)-
% (IAL:1/31/17)
while feof(mfile) == 0
    id = fgetl(mfile);
    if (isempty(strfind(id,'='))==0)
        strtemp=id;
        strtemp(strfind(strtemp, '=')) = [];
        Key1='t'; 
        Key2='N';
        Index1 = strfind(strtemp, Key1);
        Index2 = strfind(strtemp, Key2);
        timestep(ii) = sscanf(strtemp(Index1(1) + length(Key1):end), '%g', 1);
        mnum(i) = sscanf(strtemp(Index2(1) + length(Key2):end), '%g', 1);
        
        ii = ii+1; % Increment timestep counter every loop
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
ii = 1; % Separate counter for timestep
% Guarantees proper incrementing for simulations with zero motor or
% crosslinker densities (a_motor_density = 0 OR p_motor_density = 0)-
% (IAL:1/31/17)
while feof(pfile) == 0
    id = fgetl(pfile);
    if (isempty(strfind(id,'='))==0)
        strtemp=id;
        strtemp(strfind(strtemp, '=')) = [];
        Key1='t'; 
        Key2='N';
        Index1 = strfind(strtemp, Key1);
        Index2 = strfind(strtemp, Key2);
        timestep(ii) = sscanf(strtemp(Index1(1) + length(Key1):end), '%g', 1);
        pnum(i) = sscanf(strtemp(Index2(1) + length(Key2):end), '%g', 1);
        
        ii = ii+1; % Increment timestep counter every loop  
    else   
           for j = 1 : 1: pnum(i)-1
               pdata(j,:,i) = str2num(fgetl(pfile));
           end
           i=i+1;
    end
   
end

fclose(pfile);
if ispc
    save([pwd,'\simdata.mat']); % Added by IL - 10/25/16
elseif isunix
    save([pwd,'/simdata.mat']);
end