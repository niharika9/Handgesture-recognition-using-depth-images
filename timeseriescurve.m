function [ values ] = timeseriescurve(ges, inum)

filename = strcat('E:\8th sem\Computer vision\project\Image\P1\G',num2str(ges),'\',num2str(inum),'.jpg');
filename1 = strcat('E:\8th sem\Computer vision\project\Depth\P1\G',num2str(ges),'\',num2str(inum),'.txt');


img = imread(filename);
hand_img = img(:,31:end,:);

im_double = im2double(hand_img);
hand_copy = zeros(size(im_double));

depthvalues = csvread(filename1);
xmin=610;
xmax=0;
ymin=480;
ymax=0;
for r=1:1:480
    for c=1:1:610
        if depthvalues(r,c)>=600 && depthvalues(r,c)<= 900
           if c < xmin 
               xmin=c;
           end
           
           if c > xmax 
              xmax=c;
           end
           
            if r < ymin 
               ymin=r;
            end
            
            if r > ymax 
               ymax=r;
            end
                
            hand_copy(r,c,:) = im_double(r,c,:) ;
        else
            hand_copy(r,c,:) = [1 1 1];
        end
    end
end

cropped_hand = imcrop(hand_copy,[xmin-20,ymin-20,xmax+50-xmin,ymax+20]);
%cropped_hand = hand_copy(ymin-20:ymax+20,xmin-20 : xmax+20-xmin,:);

%  figure;
%  imshow(hand_copy);

crophand_copy = cropped_hand;
%  figure; imshow(crophand_copy);

mask = ones(3,3);
%convert cropped to gray scale
gray_hand  = rgb2gray(cropped_hand);
size_cropped= size(gray_hand);



gray_hand1 = zeros(size_cropped(1)+2,size_cropped(2)+2);

gray_hand1(2:size_cropped(1)+1,2:size_cropped(2)+1) = gray_hand;


%% masking done to detect the belt
for r = 2:1:size_cropped(1)+1
    for c = 2:1:size_cropped(2)+1
        overlap = mask.*gray_hand1(r-1:r+1,c-1:c+1);
        if sum(sum(overlap))/9 < 0.10
            cropped_hand(r-1,c-1,:) = [1 1 1];
        else
            cropped_hand(r-1,c-1,:) = [0 0 0];
        end
    end
end
% figure;
% imshow(cropped_hand);


gray = im2bw(cropped_hand);

cm = bwlabel(gray);
%figure; imshow(gray);


labels = unique(cm);

nums = histc(cm(:),labels);

maxlabel = max(nums(2:length(nums)));

labelindex = find(maxlabel == nums);

crop_belt = ones(size_cropped);
for i = 1 : size_cropped(1)
for j = 1: size_cropped(2)
   if cm(i,j) == (labelindex-1)
       crop_belt(i,j) = 0;
   end
end
end

 %figure; imshow(crop_belt);

% canny edge detection
% edgeDetect = edge(crop_belt,'sobel');
% figure; imshow(edgeDetect);


% to find xmin , xmax ymin ,ymax for diagonal method
beltsize = size(crop_belt);
x = [];
y = [];
check =0;
for i = 1 : beltsize(1)
    for j = 1 : beltsize(2)
   if crop_belt(i,j) == 0
      y = [y i];
      startx = j;
      check = 1;
      break;
   end
   if(check == 1)
       break;
   end
    end
end

check=0;
for i =  beltsize(1):(-1):1
    for j = 1 : beltsize(2)
   if crop_belt(i,j) == 0
      y = [y i];
      check = 1;
      break;
   end
   if(check == 1)
       break;
   end
    end
end
check =0;

for j = 1 : beltsize(2)
   for i = 1 : beltsize(1)
   if crop_belt(i,j) == 0
      x = [x j];
      check = 1;
      break;
   end
   if(check == 1)
       break;
   end
    end
end
check =0;
 for j = beltsize(2) : (-1) : 1
   for i = 1:beltsize(1)
   if crop_belt(i,j) == 0
      x = [x j];
      check = 1;
      break;
   end
   if(check == 1)
       break;
   end
  end
 end

 
 %drawing line from xmin,ymin to xmax ymax to crop the hand
 
%  this method is using drawing a diagonal line
 dummy_crophand =crophand_copy;
 slope = (y(2) - y(1))/(x(2) - x(1));
 for r = 1 : size_cropped(1)
  for c = 1 : size_cropped(2)
 
      temp =  (r - y(1))- slope * (c - x(1));   
      if temp > 0
        dummy_crophand(r,c,:) = [0 0 0];
      end
  end
 end
  %figure; imshow(dummy_crophand);  
 
 % to get the binary image of only the boundary and the diagnol  
 dummy_crophand2 = zeros(size_cropped);
 for r = 1 : size_cropped(1)
  for c = 1 : size_cropped(2)
      if dummy_crophand(r,c,1) == 1 && dummy_crophand(r,c,2) == 1 && dummy_crophand(r,c,3) == 1        
       dummy_crophand2(r,c) = 0;
      elseif dummy_crophand(r,c,1) ~= 0 && dummy_crophand(r,c,2) ~= 0 && dummy_crophand(r,c,3) ~= 0 
          dummy_crophand2(r,c) = 1;
      end     
  end
 end
 
 %figure; imshow(dummy_crophand2);
final_crop = medfilt2(dummy_crophand2, [10 10]);
 
 bndist1 = bwdist(~final_crop);
 %figure, imshow(bndist1,[]), title('Distance transform of ~bw');
 boundary1 = ~(edge(final_crop,'Sobel'));
 %figure; imshow(boundary1);

 
 % to find red point
red = [y(1) x(1)];
% check12 =0;
% for r = 1 : size_cropped(1)
%     for c = 1: size_cropped(2)
% %        temp = (r - mid1(1)) - slope * (c - mid1(2));  This is for mid
% %        line
%           temp = (r - y(1)) - slope * (c - x(1));  
%      if boundary1(r,c) == 0 && temp ==0
%         red = [r c]; 
%         check12 =1;
%          break;
%      end
%     end
%     if check12 == 1;
%         break;
%     end
% end

%figure;imshow(final_crop);
 
 
 bndist = bwdist(~final_crop);
 %figure, imshow(bndist,[]), title('Distance transform of ~bw');
 
 
 %% to find the centre point
 cd =  max(max(bndist));

[m1, n1] = size(bndist);
c_x=[];
c_y=[];
for i = 1:m1
 for j = 1:n1
    if cd == bndist(i,j)
        %count = count +1;
        c_x = [c_x j];
        c_y = [c_y i];       
    end
 end
end

% This is the center point.
cyan = [c_y(1) c_x(1)];


% sobeledge for time representation curve
boundary = ~(edge(final_crop,'Sobel'));
%figure; imshow(boundary);

% time series curve 


% This red is the starting point
 %red = [y(1) x(1)]; 
% other_red is the ending point of line 
% other_red = [y(1) x(2)];



% check12 =0;
% for r =  size_cropped(1) : (-1) : 1
%     for c = 1: size_cropped(2)
%        temp =  (r - mid1(1))- slope * (c - mid1(2)); 
%      if boundary(r,c) == 0 && temp ==0
%         other_red = [r c]; 
%         check12 =1;
%          break;
%      end
%     end
%     if check12 == 1;
%         break;
%     end
% end


%% to draw the time series curve 
v1 = [red(1,1) - cyan(1,1)   red(1,2)-(cyan(1,2))];
 d1 = norm(v1); %sqrt(v1(1,1)^2 + v1(1,2)^2);    
  
valuex =  []; 
valuey =  [];
rpixel = [];
cpixel = [];


 slope_commonarm = (red(1)-cyan(1))/(red(2)-cyan(2));
for r = 1: size_cropped(1) 
for c = 1: size_cropped(2)
    temp =  (r - y(1))- slope * (c - x(1));
     if boundary(r,c) == 0 && temp < 0         
                  
         v2 = [r-cyan(1,1)  c-cyan(1,2)];
          d2 = norm(v2); %sqrt(v2(1,1)^2 + v2(1,2)^2);
          v = dot(v1,v2);            
            % dotproduct  = [dotproduct  , v/(d1*d2)]
            dotproduct = v/(d1*d2);
            if dotproduct > 1.00 
               dotproduct = 1.0;
            elseif dotproduct < (-1.00)
                dotproduct = -1.0;
            end
          ang = acosd(dotproduct);          
          commonline = (r - cyan(1)) - slope_commonarm*(c - cyan(2));
           if commonline > 0
               ang = 360 - ang;
           end     
              
          ang_n = ang/360;
         d = sqrt( (cyan(1,2) - c )^2 + (cyan(1,1) - r)^2 );
         d_n = d/cd; 
        valuex =  [valuex ang_n]; 
        valuey =  [valuey d_n];
        
     end
end
end

% slope = (mid1(1) - mid2(1)) / (mid1(2) - mid2(2));
dumb = zeros(size_cropped);
for r = 1:size_cropped(1)
    for c = 1 : size_cropped(2)
        temp = (r - y(1)) - slope * (c - x(1));
        if temp==0
            dumb(r,c)=1;
        end
    end
end
%figure;imshow(dumb);




%%  ts = time series
 values = [valuex',valuey'];
 sort_val = sortrows(values, [1,2]);

% figure;
% ts = plot(valuex,valuey,'.','Color','k');

end

