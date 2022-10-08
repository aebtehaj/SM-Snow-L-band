function cmap = cmapgen(CMat,steps)
% creats colormaps
% interpolates linearly colors between the rows of CMat
% steps are the number of inrepolation points
% length(steps) vector need to be equal to size(CMat,1)-1
cfun =@(x,y,stps) [(linspace(x(1),y(1),stps)') (linspace(x(2),y(2),stps)') (linspace(x(3),y(3),stps)')];
cmap = [];
for i=1:size(CMat,1)-1
    cmap = [cmap(1:end-1,:);cfun(CMat(i,:),CMat(i+1,:),steps(i))];
end
end