%% Format the controller such that it can be loaded into Labview
% clear, close all

% load K

num = K.num; den = K.den;

Knum_arr = [];
for i=1:size(num,1)
    for j=1:size(num,2)
        Knum_arr(end+1,:) = num{i,j};
    end
end

Kden_arr = [];
for i=1:size(den,1)
    for j=1:size(den,2)
        Kden_arr(end+1,:) = den{i,j};
    end
end

formatstring = '%0.6f';
for i=2:size(den{1,1},2)
    formatstring = [formatstring ' \t %0.6f']; 
end    
formatstring = [formatstring '\n'];
fid = fopen('Knum.txt','wt');
fprintf(fid,formatstring,Knum_arr');
fclose(fid);

fid = fopen('Kden.txt','wt');
fprintf(fid,formatstring,Kden_arr');
fclose(fid);