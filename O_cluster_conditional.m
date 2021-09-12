clear all
% analyze_results
%cell_name2 is the csv file name of molecular list from thunder storm
%analysis
cell_name2='Gad_cell_day2__cld488_zo1_647_Channel_488_1017.csv'
%initial definition of ring size and the range size 
global max_R linind_rings
pixel_size = 20; % nm
max_R = 50;
max_y_gr_plot = 5; % used for plots
ms = 1; % MarkerSize
rings = rings2array(max_R);
linind_rings = linindrings(rings);
        nm_data2 = readnmdata(cell_name2);
         crop_nm_coordinates = cropnmcoordinates(nm_data2,nm_data2);
            crop_nm_x_min = crop_nm_coordinates(1);
    crop_nm_x_max = crop_nm_coordinates(2);
    crop_nm_y_min = crop_nm_coordinates(3);
    crop_nm_y_max = crop_nm_coordinates(4);
    croped_nm_data2 = crop_data(nm_data2,crop_nm_x_min,crop_nm_x_max,crop_nm_y_min,crop_nm_y_max);
 
D2=40;        %Distance parameter of matlab function "cluster"
n_inclust=10; %Define a big cluster: "n_inclust" is the minimum points in the big cluster

%C1=croped_nm_data1;   %the green ponts
C2=croped_nm_data2;    %the red points

Y = pdist(C2(:,1:2));   
clear T;
clear mr;
Z = linkage(Y);
T = cluster(Z,'cutoff',D2,'criterion','distance');
num_clusters2=max(T);


xy_clust=[];
xy_no_clust=[];
for m2=1:num_clusters2
    ncm2=find(T==m2);
    x2=C2(ncm2,1); y2=C2(ncm2,2);
    xc2=mean(x2); yc2=mean(y2);
        if length(ncm2)==2
        hold on; plot(x2,y2,'-m');
    end
      if length(ncm2)==1
        hold on; plot(x2,y2,'.b');
    end
        
        mr(m2,2)=length(ncm2);%sum(inpolygon(x2(k2),y2(k2),cr(:,1),cr(:,2)));
       if length(ncm2)>2
        k2 = convhull(x2,y2);
        mr(m2,1)=polyarea(x2(k2),y2(k2));
     
         hold on; plot(x2,y2,'.');

       end
       if length(ncm2)>n_inclust-1
                xy_clust=[xy_clust;[x2,y2]];
       else
           xy_no_clust=[xy_no_clust;[x2,y2]];
       end
      
end    
    

figure;hist(mr(:,2),max(mr(:,2)))

up_to_handred=[]

 for j = 1:length(mr)
     
  j = find(([mr(:,2)]  >=  10))
  up_to_handred(:,:) = mr(j,:) ;   

 end
 
save('Gad_cell_day2__cld488_zo1_647_Channel_647_1017_for_image_simolation.mat','D2','num_clusters2', 'mr','up_to_handred')