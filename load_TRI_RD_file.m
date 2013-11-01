function [ c1a,ac2,c1c2,id,n_c1a,n_ac2,n_c1c2,cover ] = load_TRI_RD_file(file_name)
% load triplet junction read distribution file
%
% Leo J Lee, University of Toronto, 2013

    fid=fopen(file_name);
    temp=textscan(fid,'%s%n%n%n%s%s%s%s','Delimiter','\t','CommentStyle','#');
    fclose(fid);
    id = temp{1}; n = numel(id);
    n_c1a = temp{2}; n_ac2 = temp{3}; n_c1c2 = temp{4};
    cover = (n_c1a + n_ac2)/2 + n_c1c2;
    temp2 = textscan(temp{5}{1},'%d','Delimiter',','); p = numel(temp2{1});
    c1a = zeros(n,p); ac2 = c1a; c1c2 = c1a;
    for i=1:n
        temp2=textscan(temp{6}{i},'%n','Delimiter',',');
        c1a(i,:)=temp2{1};
        temp2=textscan(temp{7}{i},'%n','Delimiter',',');
        ac2(i,:)=temp2{1};
        temp2=textscan(temp{8}{i},'%n','Delimiter',',');
        c1c2(i,:)=temp2{1};
    end
    
end