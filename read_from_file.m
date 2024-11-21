function [p1, p2, v1, v2] = read_from_file(path, type, n)
    if(strcmp(type,"triangle"))
        num_cp = (n+1)*(n+2)/2;
    elseif(strcmp(type,"patch"))
        num_cp = (n+1)*(n+1);
    else
        assert(false)
    end
    fileID = fopen(path,'r');
    p1 = fread(fileID, [3 num_cp], 'double');
    p2 = fread(fileID, [3 num_cp], 'double');
    v1 = fread(fileID, [3 num_cp], 'double');
    v2 = fread(fileID, [3 num_cp], 'double');
    fclose(fileID);
end