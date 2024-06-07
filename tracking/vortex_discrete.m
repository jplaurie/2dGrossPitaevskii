
th = 10.0;
xlocs =[];
ylocs =[];
pol=[];
Gamma=load('Circ.00000');

h=fspecial('gaussian',size(Gamma),1);
Gamma = imfilter(Gamma,h);

negareas = bwlabel(Gamma < -th);
posareas = bwlabel(Gamma > th);

for i =1:max(max(posareas))
    [r,c] = find(posareas == i);
    if(length(r) >1)
        xlocs = [xlocs;c];
        ylocs = [ylocs;r];
        pol = [pol;1];
    end
end

% for i =1:max(max(negareas))
%    [r,c] = find(negareas ==i);
%    if(length(r) >1)
%        xlocs = [xlocs;c];
%        ylocs = [ylocs;r];
%        pol = [pol;-1];
%    end
% end
scatter(xlocs,ylocs);