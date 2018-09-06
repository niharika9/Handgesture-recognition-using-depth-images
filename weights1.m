function [area ,new_points3] = weights1( values )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% to find the  threshold lines of each fingers 
plot(values(:,1),values(:,2),'.k');
pointsx = [];
pointsy = [];

%sort_val = sortrows(values, [1,2]);
valsize = size(values);

for i = 1 : valsize(1,1)
    if ( values(i,2) > 1.6 )  && ( values(i,2) <= 1.65) && ( values(i,1) > 0.15) && (values(i,1) < 0.7 )
       pointsx = [pointsx values(i,1)];
       pointsy = [pointsy values(i,2)];
    end
end

points = [pointsx' pointsy'];
sort_p = sortrows(points, [1,2]);

for i = 1 : (length(sort_p(:,1))-1)
    for j = i+1 : (length(sort_p(:,1)))
        if abs(sort_p(i,1) - sort_p(i+1,1)) < 0.01       
            sort_p(i+1,1) = sort_p(i,1);
            sort_p(i+1,2) = sort_p(i,2);
        else
            break;
        end
    end
end

new_points = unique(sort_p,'rows');

hold on;


%% to find correct new points  without any errors

 
% rangeup = [1.5 1.55];
% rangedown = [1.45 1.5];
for i = 1:length(new_points(:,1))
    countup=0;
    countdown=0; 
    rangeup = [new_points(i,2) new_points(i,2) + 0.08];
    rangedown = [new_points(i,2)-0.08 new_points(i,2)];
    
    rangex = [new_points(i,1)-0.005  new_points(i,1)+0.005];
    pts= [];
    for j = 1 : valsize(1,1)
        if (values(j,1) > rangex(1,1) && values(j,1) < rangex(1,2)) && (values(j,2) > rangeup(1,1)) && (values(j,2) < rangeup(1,2))        
          countup = countup + 1;
          %pts = [ pts; [values(j,1) values(j,2)]] ; 
        elseif  (values(j,1) > rangex(1,1) && values(j,1) < rangex(1,2)) && (values(j,2) > rangedown(1,1)) && (values(j,2) < rangedown(1,2))
          countdown = countdown + 1;
        end
    end
    
    if countup==0 && (countdown > 0)
       new_points(i,:) = [0 0]; 
    end
    if (countup > 0) && (countdown== 0)
        new_points(i,:) = [0 0]; 
    end

end

% to eliminate non zero values 
new_points3 = new_points(any(new_points,2) , :);

for i = 1 : 2 : length(new_points3(:,1))
   plot([new_points3(i,1) new_points3(i+1,1)] , [new_points3(i,2) new_points3(i+1,2)], '-');
end

%% to find the area of the fingers 
area_points = cell(length(new_points3(:,1))/2,1);
 
px =[];
py =[];
for i = 1 : 2 : length(new_points3(:,1))
    for j = 1 : valsize(1,1)
      if values(j,2) > 1.6  && ( values(j,1) < new_points3(i+1,1)+ 0.02 && values(j,1) > new_points3(i,1)-0.02)
        px = [values(j,1) px];
        py = [values(j,2) py];
      end  
    end
    
    z = (i+1)/2;
    area_points{z,1} = [px' py'];
    area_points{z,1} = sortrows(area_points{z,1},[1,2]);
  
   px = [];
   py = [];
end

area = [];
for i = 1:2:length(new_points3(:,1))
    z = (i+1)/2;
    height = new_points3(i+1,1)-new_points3(i,1);
    sum_of_base = new_points3(i+1,2)+new_points3(i,2);
    area_remain = 0.5*height*sum_of_base;
   
    area(z) = trapz(area_points{z,1}(:,1),area_points{z,1}(:,2))-area_remain;
end

end

