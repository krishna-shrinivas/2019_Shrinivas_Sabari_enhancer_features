%%  Takes in the path  to the excel file (file_name) which is an excel spreadsheet
%   containing the data regarding channel information and unique folders.
%   The function returns all_params --> a cell array of params for each run
%   of recruitment_analysis_average_data.m

function [all_params] = parse_data(file_name)

input_data = readtable(file_name);
all_params = {};

for i=1:1:height(input_data)
    all_params{i} = [];
    all_params{i}.channel_names = {};
    all_params{i}.wavelength_channels = {};
    all_params{i}.folder_with_images = strjoin([ fileparts(file_name) '/' input_data.folder_path(i) '/'],'');
    if iscellstr(input_data.Channel_488(i))
        all_params{i}.wavelength_channels = [all_params{i}.wavelength_channels '488'];
        all_params{i}.channel_names = [all_params{i}.channel_names input_data.Channel_488(i)];
        
    end
    
    
    if iscellstr(input_data.Channel_561(i))
        all_params{i}.wavelength_channels = [all_params{i}.wavelength_channels '561'];
        all_params{i}.channel_names = [all_params{i}.channel_names input_data.Channel_561(i)];
        
    end
    
    
    if iscellstr(input_data.Channel_640(i))
        all_params{i}.wavelength_channels = [all_params{i}.wavelength_channels '640'];
        all_params{i}.channel_names = [all_params{i}.channel_names input_data.Channel_640(i)];
        
    end
    
    %%  Default settings for image analysis
    %   Threshold multiplier for intensity thresholding
    %   Draw_figure for plotting outlines on identified droplets
    %   minimum_drop_size for size thresholds
    %   average_scaffold if the average of images is to be used
    %   bg_subtract for background subtraction
    %   bg_value is the value of the background signal that is subtracted
    all_params{i}.threshold_multiplier = 3;
    all_params{i}.draw_figure = 1;
    all_params{i}.minimum_drop_size = 9;
    all_params{i}.average_scaffold = 0;
    all_params{i}.scaffold_channel = 2;
    all_params{i}.bg_subtract = 0;
    all_params{i}.bg_value = 93;
    all_params{i}.bg_size_flag = 0;
end

end