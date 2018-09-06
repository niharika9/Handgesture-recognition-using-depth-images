% Taking test images and finding the gesture recognition.


tic;
recall = [];
totalcorrect = 0;

for iter = 1 : 5

    
    %generating random values for getting test images 
    random = [];

    for i = 2:6
    if i==2 || i==3 || i == 6
        r = randperm(10,2);
        random = [random;i,r(1),r(2)];
    elseif i == 4
        r = randperm(9,2);
        random = [random;i,r(1),r(2)];
    elseif i == 5
        r = randperm(8,2);
        random = [random;i,r(1),r(2)];
    end
    end

% Retrieving the respective images as test images and then finding their
% time series curve and area values;
    count = 0;

    for i = 1:5
        ges = random(i,1);
        for j = 2:3
            count = count+1;
            img = random(i,j);
            [values] = timeseriescurve(ges,img);
            [area , newpoints] = weights1(double(values));
            save(['test_' num2str(count) '.mat'] ,'newpoints','ges','values','area');
        end
    end

    
   
% Calculating femd values for each test image with the rest of the images and then finding 
% the least femd value and thus, its corresponding gesture.

    
    output = cell(10,2);
    for i = 1:10
    output{i,1} = zeros(1,37);
    output{i,2} = zeros(1,37);
    end
    
    correctans = 0;
    for i = 1:10
        file  = [strcat('test_',num2str(i),'.mat')];
        load(file);
        area1 = area;
        newpoints1= newpoints;
        testges = ges;
        count = 0;
        for gesnum = 2:6
           if gesnum==2 || gesnum==3 || gesnum==6
               for imnum = 1 :10
                  if imnum ~= random(gesnum-1,2) && imnum ~= random(gesnum-1,3)
                      datapathfile = strcat('points_',num2str(gesnum),num2str(imnum),'.mat') 
                    load(datapathfile);
                    area2 = area;
                    newpoints2= newpoints;
                    fmval = fingerEMD(area1,newpoints1,area2,newpoints2) ; 
                 count = count +1;
                  output{i,1}(1,count) = fmval;
                  output{i,2}(1,count) = gesnum;
                  end  
               end


           elseif gesnum == 4
               for imnum = 1:9

                   if imnum ~= random(gesnum-1,2) && imnum ~= random(gesnum-1,3)
                    load(strcat('points_',num2str(gesnum),num2str(imnum),'.mat'));
                    area2 = area;
                    newpoints2= newpoints;
                     fmval = fingerEMD(area1,newpoints1,area2,newpoints2)  ;
                 count = count +1;
                  output{i,1}(1,count) = fmval;
                  output{i,2}(1,count) = gesnum;
                  end


               end
           elseif gesnum == 5
               for imnum = 1:8

                   if imnum ~= random(gesnum-1,2) && imnum ~= random(gesnum-1,3)
                    load(strcat('points_',num2str(gesnum),num2str(imnum),'.mat'));
                    area2 = area;
                    newpoints2= newpoints;
                    fmval = fingerEMD(area1,newpoints1,area2,newpoints2) ; 
                 count = count +1;
                  output{i,1}(1,count) = fmval;
                  output{i,2}(1,count) = gesnum;
                  end

               end
           end  
                    
        end

        final_ges = min(output{i,1}(1,:));
        index = find(final_ges == output{i,1}(1,:));
     % value of observed gesture = output{i,2}(1,index)-1;
      %value of actual gesture = testges-1;

      opts1=  optimset('display','off');
       
      if (output{i,2}(1,index)-1) == (testges-1)
          correctans = correctans + 1;
          totalcorrect = totalcorrect + 1; 
       end

    end

    recall = [recall ; correctans/10*100];

end
prog_recall = totalcorrect/(iter*10) * 100;
close all;
toc;

