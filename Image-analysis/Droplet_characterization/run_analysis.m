
%% Default initialization input parameters
input_parameters = [];
input_parameters.type_of_protein = [];
input_parameters.average_scaffold = [];
input_parameters.minimum_drop_size = [];
input_parameters.scaffold_channel = [];
input_parameters.bg_subtract = [];
input_parameters.bg_value = [];
input_parameters.draw_figure = [];
%% User intialization of input parameters
names = {'MED1_OCT4','ODNA_MED1_OCT4'};

for name=1:1:length(names)
%   Path to excel sheet with information
input_parameters.file_name = ['/media/krishna/Windows/MIT/Research/Transcription/How do SEs nucleate/In_vitro_analysis_Ben/three_component_droplets/Background_subtraction/20181126_plusminusDNA/data_information_' names{name} '.xlsx'];

%   Output will be stored under output_'type_of_protein', defaults to
%   Output_default. Used to identify different classes of outputs. Ideally each
%   excel file will have a unique output file.
input_parameters.type_of_protein = names{name};

%   These are parameters that the user can manually set to override default
%   manuals. Defaults can be found under parse_data. Briefly, the settings:
%   draw_figure for plotting outlines on identified droplets
%   minimum_drop_size for size thresholds
%   average_scaffold if the average of images is to be used
%   scaffold_channel decides which channel to use as scaffold : 1 -488, 2-
%   561, 3-640
%   bg_subtract for background subtraction
%   bg_value is the value of the background signal that is subtracted

input_parameters.average_scaffold = 0;
input_parameters.minimum_drop_size = 9;
input_parameters.scaffold_channel = 2;
input_parameters.bg_subtract =1;
input_parameters.bg_value = 80;
input_parameters.draw_figure = 0;
input_parameters.bg_size_flag = 1;

call_input_parse_data(input_parameters);
end