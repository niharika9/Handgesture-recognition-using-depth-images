function [femd_value] = fingerEMD(area1, new_points1, area2, new_points2)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

featurevector1 = [];
featurevector2 = [];
for i = 1:2:length(new_points1(:,1))
featurevector1 = [[featurevector1] ; [new_points1(i,1) , new_points1(i+1,1) ] ];
end

for i = 1:2:length(new_points2(:,1))
featurevector2 = [[featurevector2] ; [new_points2(i,1) , new_points2(i+1,1)] ];
end

[x, fval] = emd(featurevector1,featurevector2 ,area1',area2', @gdf);


%% ground distance matrix

d = gdm(featurevector1,featurevector2,@gdf);

%% Finger earth mover's distance 
l1 = length(new_points1(:,1))/2;
l2 = length(new_points2(:,1))/2;
df = 0;
fv = 0;
for i = 1 : l1 * l2
       df = d(i,1) * x(i,1) + df;
       fv = fv + x(i,1);
end

wsum1 = sum(area1');
wsum2 = sum(area2');

beta = 0.5;
femd_value = ( beta * df + ((1 - beta) * abs(wsum1 - wsum2)))/ fv;



end

