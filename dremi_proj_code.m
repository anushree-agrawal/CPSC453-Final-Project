load adrenaldata.csv
%% Transpose
adrenaldata = adrenaldata';
%% Write FCS
writeFCS('adrenalfcs', adrenaldata)

%% Dremi
simpledremi

%%
load heart_data.csv

%% Transpose
heart_data = heart_data';
%% Write FCS
writeFCS('heartdata.fcs', heart_data);

%dremi
%% Get data from image
D=get(gca,'Children');
CData = get(D, 'CData');

%% sort and find top values
[index, value] = sort(CData(:));
index = flipud(index);
value = flipud(value);

%% Find best pos
[m,n] = size(CData);
pos = [];
for i=1:20
    pos(1, i) = ceil((value(i)/m));
    pos(2, i) = mod(value(i), m);
    if (pos(2, i) == 0)
        pos(2, i) = m;
    end
end
pos