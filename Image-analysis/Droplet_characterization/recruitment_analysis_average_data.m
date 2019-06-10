%%% recruitment_analysis_average_data takes in two key variables: path to
%%% the folder, and the params variable which contains information on
%%% channels involved, channel names, droplet_size conditions, and
%%% intensity threshold. The output folder is where the excel data is
%%% stored.
        

function [IT,T,Q,T_total] = recruitment_analysis_average_data(folder_with_images, params,output_folder)


clc;
close all;


%   Store relevant channel names and input channel wavelengths
channel_names = params.channel_names;
wavelength_channels = params.wavelength_channels;
N_channels = length(channel_names);

%   Get list of images
list_of_files = dir(folder_with_images);
N_files = size(list_of_files,1);
list_of_image_files = {};
count =1;
write_file_names = {};


%   Identify a list of unique images with first channel name
for i=3:1:N_files
    if (~isempty(strfind(list_of_files(i).name,wavelength_channels{1})) && ~isempty(strfind(lower(list_of_files(i).name),'tif')))
        list_of_image_files{1,count} = list_of_files(i).name;
        write_file_names{count} = list_of_files(i).name;
        count = count+1;
    end
end


%   Gather other images with same unique prefix from other channels.
for c= 2:1:length(wavelength_channels)
    count = 1;
        
    %   This is the format for files that Ann et al have sent, now alsoo
    %   used by Ben Sabari for the three  component analysis. Unique ID is
    %   the prefix before _w in image file name
    for i=1:1:size(list_of_image_files,2)
        
        prefix = strsplit(list_of_image_files{1,i},'_w');
 
        for j=3:1:N_files
    
            if (~isempty(strfind(list_of_files(j).name,prefix{1})) && ~isempty(strfind(list_of_files(j).name,wavelength_channels{c})))
    
                list_of_image_files{c,i} = [folder_with_images list_of_files(j).name];
            end
        end
        
    end
end


%   Data structures to calculate average and total concentration
track_average_conc = zeros(size(list_of_image_files,2),2*length(wavelength_channels));
track_total_conc = zeros(size(list_of_image_files,2),2*length(wavelength_channels));

minimum_drop_size  = params.minimum_drop_size ; % In pixel^2


store_area = [];
store_ar = [];
store_COM_intensity = [];
image_file_temp = {};
store_bulk_intensity = [];

for i=1:1:size(list_of_image_files,2)
    
    %   First, create an average image across relevant channels
    list_of_image_files{1,i} = [folder_with_images list_of_image_files{1,i}];

    [avg_image,map] =  imread( list_of_image_files{1,i});
    C_store{1,i} =avg_image;
    for c = 2:1:length(wavelength_channels)

        [X,map] = imread(list_of_image_files{c,i});
        avg_image = avg_image + X;
        C_store{c,i} = X;
        
    end
    
    if params.average_scaffold
        X = avg_image/3;
    else
        X = C_store{params.scaffold_channel,i};
    end
    
    %   Binarize into black and white using global thresholds
    threshold_multiplier =params.threshold_multiplier;
    
%     bw = imbinarize(X,mean(mean(X))*threshold_multiplier/65535);
        
    if params.bg_subtract
        X  = X - params.bg_value;
        X(X<0) = 0;
        for c = 1:1:length(wavelength_channels)
           C_store{c,i}= (C_store{c,i}-params.bg_value) ;
           C_store{c,i}(C_store{c,i}<0) = 0 ; 

        end
    end
    flat_data = double(reshape(X,length(X)^2,1));
    [vals,bins] = histcounts(flat_data);
    [max_p,max_value] = max(vals);
    pos_of_bg_peak= bins(max_value);
    bw = imbinarize(X,(pos_of_bg_peak+ std(flat_data)*threshold_multiplier)/65535);
    bw1 = bw;
%     bw = imbinarize(X,mean(mean(X))*threshold_multiplier/65535);


    bw = bwareaopen(bw,minimum_drop_size);
    
    [B,L] = bwboundaries(bw,'noholes');
    
    %     Display the label matrix and draw each boundary
    if params.draw_figure
        figure;
        %     imshow(label2rgb(L, @jet, [.5 .5 .5]))
        imshow(X,[]);
        hold on
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2); hold on;
        end
        drawnow;
    end
    
    %   Get size & A.R. distribution based on "average" channel
    [ar_relevant,area_relevant,COM] = calculate_size_distribution(B,L,params);
    COM_temp = [];
    if ~isempty(COM)
        for channel =1:1:length(wavelength_channels)
           COM_temp = [COM_temp diag(C_store{channel,i}(floor(COM(:,2)),floor(COM(:,1)))) ];
        end
        store_COM_intensity = [store_COM_intensity; COM_temp];

    end
    store_area = [store_area area_relevant];
    store_ar = [ store_ar ar_relevant];
    o = cell(1,length(ar_relevant));
    o(:) = write_file_names(i);
    image_file_temp = [image_file_temp o];
    
    pixels_inside_droplets = find(bw>0);
    if ~params.bg_size_flag
        pixels_outside_droplets = find(bw==0);
    else
        pixels_outside_droplets = find(bw1==0);

    end
    
    %   Read each TIFF fike
    for c= 1:1:length(wavelength_channels)
        [X,map] = imread(list_of_image_files{c,i});
        track_average_conc(i,2*c-1) = mean(C_store{c,i}(pixels_inside_droplets));
        track_total_conc(i,2*c-1) = sum(double((C_store{c,i}(pixels_inside_droplets))));
        
        if isempty(pixels_inside_droplets)
            track_average_conc(i,2*c-1) = 0;
            track_total_conc(i,2*c-1) = 0;

        end
        track_average_conc(i,2*c) = mean(C_store{c,i}(pixels_outside_droplets));
        track_total_conc(i,2*c) = sum(double(C_store{c,i}(pixels_outside_droplets)));

    end
    N = ones(1,length(ar_relevant))*track_average_conc(i,2*params.scaffold_channel);
    store_bulk_intensity = [store_bulk_intensity N];
    
end
% close all;

%   Giving a unique_name to the output_excel file which has the names of
%   the unique channels.

output_file = [output_folder 'recruitment_size_'];

for c = 1:1:length(channel_names)
    output_file = [output_file channel_names{c} '_'];
end

output_file = [output_file 'analysis.xlsx'];



params_relevant = [];
variable_names = {};

params_relevant_total = [];

for c=1:1:length(wavelength_channels)

    % Assigning the information on average intensity inside and outside the
    % droplet
    params_relevant = [params_relevant track_average_conc(:,2*c-1) track_average_conc(:,2*c) track_average_conc(:,2*c-1)./track_average_conc(:,2*c) (track_average_conc(:,2*c-1)./track_average_conc(:,2*c))./(1+track_average_conc(:,2*c-1)./track_average_conc(:,2*c))];
    variable_names = [variable_names ,[channel_names{c} '_droplet' ],[channel_names{c} '_bulk' ],[channel_names{c} '_ratio' ],[channel_names{c} '_division_ratio']]; 

    
    %   Assigning the information on the total intensity inside and
    %   "outside" the droplets
    params_relevant_total = [params_relevant_total track_total_conc(:,2*c-1) track_total_conc(:,2*c) track_total_conc(:,2*c-1)./track_total_conc(:,2*c) (track_total_conc(:,2*c-1)./track_total_conc(:,2*c))./(1+track_total_conc(:,2*c-1)./track_total_conc(:,2*c))];

end

l = size(params_relevant,1);
params_relevant(l+1,:) = mean(params_relevant(1:l,:));
params_relevant(l+2,:) = std(params_relevant(1:l,:));

params_relevant_total(l+1,:) = mean(params_relevant_total(1:l,:));
params_relevant_total(l+2,:) = std(params_relevant_total(1:l,:));

T = array2table (params_relevant);
T_total =array2table(params_relevant_total);
if ~isempty(params_relevant)

    
    T(:,4*length(wavelength_channels)+1) = [write_file_names'; 'Mean_data' ; 'Std_data'];
    T_total(:,4*length(wavelength_channels)+1) = [write_file_names'; 'Mean_data' ;  'Std_data'];

  
    T.Properties.VariableNames = [variable_names, 'Image_file'];
    T_total.Properties.VariableNames = [variable_names, 'Image_file'];

    writetable(T, output_file,'Sheet',['Avg_I_size_' num2str(params.minimum_drop_size) '_TI_' num2str(params.threshold_multiplier)]);
    writetable(T_total, output_file,'Sheet',['Total_I_size_' num2str(params.minimum_drop_size) '_TI_' num2str(params.threshold_multiplier)]);

    IT = T;
    T = [];
    params_relevant = [store_area' store_ar' store_COM_intensity];
end

if ~isempty(params_relevant)
    T = array2table ([params_relevant store_bulk_intensity']);
    T(:,end+1) = [image_file_temp'];

    if length(wavelength_channels) ==3
        T.Properties.VariableNames = {'Area_in_pixel_square','Aspect_ratio','Intensity_C_488','Intensity_C_561','Intensity_C_640','Intensity_sc_dilute','Image_file'};
    elseif length(wavelength_channels) ==2
        T.Properties.VariableNames = {'Area_in_pixel_square','Aspect_ratio','Intensity_C_488','Intensity_C_561','Intensity_sc_dilute','Image_file'};
    else
        T.Properties.VariableNames = {'Area_in_pixel_square','Aspect_ratio','Intensity_C_488','Intensity_sc_dilute','Image_file'};

    end
    
    writetable(T, output_file,'Sheet',['Aspect_area_MDS_' num2str(params.minimum_drop_size) '_TI_' num2str(params.threshold_multiplier)]);
    Q = array2table([mean(store_area) mean(store_ar) length(store_area)]);
    Q.Properties.VariableNames = {'Mean_area_in_pixel_square','Mean_aspect_ratio','Number_of_droplets'};
    writetable(Q, output_file,'Sheet',['Aspect_area_MDS_' num2str(params.minimum_drop_size) '_TI_' num2str(params.threshold_multiplier)],'Range','H1');

    %
    %
    disp(['Number of ID droplets are ' num2str(size(store_area,2))]);
    disp(['Average area is ' num2str(mean(store_area))]);
    disp(['Average aspect ratio is ' num2str(mean(store_ar))]);
end

if isempty(params_relevant)
    T = [];
    Q = [];
end
end


%%  Function takes in identified boundaries and perimeter to generate statistics 
%   including aspect ratios, centroids, and areas.
function [ar_relevant,area_relevant,COM] = calculate_size_distribution(B,L,params)

    stats = regionprops(L,'Area','Centroid', 'MajorAxisLength','MinorAxisLength');
    track_aspect_ratio = [];
    track_area = [];
    track_COM = [];

    threshold = 0.94;
    
    % loop over the boundaries
    for k = 1:length(B)
        
        %Determine the aspect ration
        aspect_ratio(k) = stats(k).MajorAxisLength/stats(k).MinorAxisLength;
        track_aspect_ratio = [track_aspect_ratio aspect_ratio(k)];
        
        % obtain (X,Y) boundary coordinates corresponding to label 'k'
        boundary = B{k};
        
        % compute a simple estimate of the object's perimeter
        delta_sq = diff(boundary).^2;
        perimeter = sum(sqrt(sum(delta_sq,2)));
        
        % obtain the area calculation corresponding to label 'k'
        area = stats(k).Area;
        store_area(k) = area;
        track_area = [track_area store_area(k)];
        % compute the roundness metric
        metric = 4*pi*area/perimeter^2;
        store_metric(k) = metric;
        
        track_COM = [track_COM; stats(k).Centroid];
        
    end

max_area_threshold = 1000;
min_area_threshold = params.minimum_drop_size;
ar_threshold = 1.8;

idx_area_relevant = find(track_area<=max_area_threshold & track_area>=min_area_threshold);
idx_ar_relevant = find(track_aspect_ratio<=ar_threshold);
idx_relevant = intersect(idx_ar_relevant,idx_area_relevant);
ar_relevant = track_aspect_ratio(idx_relevant);
area_relevant = track_area(idx_relevant);
COM = track_COM(idx_relevant,:);




end