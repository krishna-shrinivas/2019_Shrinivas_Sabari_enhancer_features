%% call_analysis after parsing input_data_structures

function call_input_parse_data(input_parameters)

file_name = input_parameters.file_name;
all_params = parse_data(file_name);


if ~isempty(input_parameters.average_scaffold)
    average_scaffold = input_parameters.average_scaffold;
end

if ~isempty(input_parameters.bg_size_flag)
    bg_size_flag = input_parameters.bg_size_flag;
end

if ~isempty(input_parameters.minimum_drop_size)
    
    minimum_drop_size = input_parameters.minimum_drop_size;
end

if ~isempty(input_parameters.scaffold_channel)
    scaffold_channel = input_parameters.scaffold_channel;
end

if ~isempty(input_parameters.bg_subtract)
    bg_subtract = input_parameters.bg_subtract;
end

if ~isempty(input_parameters.bg_value)
    if bg_subtract
        bg_value = input_parameters.bg_value;
    else
        bg_value = 0;
    end
end

if isempty(input_parameters.type_of_protein)
    input_parameters.type_of_protein = 'default'
end

if ~isempty(input_parameters.draw_figure)
    draw_figure = input_parameters.draw_figure;
end


output_folder = [fileparts(file_name) '/Output_' input_parameters.type_of_protein '/'];

%   Create unique output folder for a set of conditions to store excel
%   files and scatter plots
if ~isdir(output_folder)
    mkdir(output_folder);
end

unique_date = date;

summary_statistics_file =  [output_folder 'summary_statistics_' unique_date '_' input_parameters.type_of_protein '.xlsx'];
Ncond = size(all_params,2);
summary_statistics = [];
write_file_names = {};




for i=1:1:Ncond
    close all;
    
    params  = all_params{i};
    params.bg_size_flag  =bg_size_flag;
    params.average_scaffold = average_scaffold;
    params.minimum_drop_size = minimum_drop_size;
    params.draw_figure =draw_figure;
    params.scaffold_channel = scaffold_channel;
    params.bg_subtract = bg_subtract;
    params.bg_value = bg_value;
    for c=1:1:length(params.channel_names)
        params.channel_names{c} = strrep(params.channel_names{c},'.','_');
        params.channel_names{c} = strrep(params.channel_names{c},' ','_');
        params.channel_names{c} = strrep(params.channel_names{c},'-','_');
        params.channel_names{c} = strrep(params.channel_names{c},'+','_');

        params.channel_names{c} = [params.channel_names{c} '_scaf_' (params.wavelength_channels{params.scaffold_channel})];
    end
    b = strsplit(params.folder_with_images,'/');
    write_file_names = [write_file_names; b(end-1)];
    disp(['The file being parsed is ' (b)]);
    
    if ~params.average_scaffold
        output_folder_sheets = [output_folder 'Scaf_' (params.wavelength_channels{params.scaffold_channel}) '_bgs_' num2str(params.bg_value) '_TI_' num2str(params.threshold_multiplier) '_MDS_' num2str(params.minimum_drop_size) '/'];
    else
        output_folder_sheets = [output_folder 'Scaf_avg_bgs_' num2str(params.bg_value) '_TI_' num2str(params.threshold_multiplier) '_MDS_' num2str(params.minimum_drop_size) '/'];
        
    end
    if ~isdir(output_folder_sheets)
        mkdir(output_folder_sheets)
    end
    
    [IT,T,Q,T_total]= recruitment_analysis_average_data(params.folder_with_images, params,output_folder_sheets);
    for c=1:1:length(params.channel_names)
        % Transformed parition coefficient
%         summary_statistics(i,6*c-5:6*c) = [table2array(IT(end-1,4*c)) table2array(IT(end,4*c)) table2array(T_total(end-1,4*c)) table2array(T_total(end,4*c)) table2array(T_total(end-1,4*c-3))+ table2array(T_total(end-1,4*c-2)) table2array(T_total(end,4*c-3))+ table2array(T_total(end,4*c-2))];
    
        % Raw partition coefficient
        summary_statistics(i,6*c-5:6*c) = [table2array(IT(end-1,4*c-1)) table2array(IT(end,4*c-1)) table2array(T_total(end-1,4*c)) table2array(T_total(end,4*c)) table2array(T_total(end-1,4*c-3))+ table2array(T_total(end-1,4*c-2)) table2array(T_total(end,4*c-3))+ table2array(T_total(end,4*c-2))];

    end
    
    if ~isempty(T)


        
        green_channel = table2array(T(:,3));
        if length(params.channel_names) >=2
            red_channel = table2array(T(:,4));
        end
        if length(params.channel_names) == 3
            DNA_channel = table2array(T(:,5));
        end
        size_droplet = table2array(T(:,1));
    end
    output_file = [output_folder_sheets 'scatter_'];
    
    for c = 1:1:length(params.channel_names)
        output_file = [output_file params.channel_names{c} '_'];
    end
    output_file = [output_file '_' num2str(params.minimum_drop_size) '_' num2str(params.scaffold_channel) '_'];
    
    if ~isempty(T)
        output_file = [output_file 'analysis.svg'];
        figure1=figure;
        set(gcf, 'Visible', 'off');

        axes1 = axes('Parent',figure1);
        
        if length(params.channel_names) == 3
            scatter3(green_channel,red_channel,DNA_channel,double(size_droplet*0.5)); hold on;
            xlabel(['GFP-' strrep(params.channel_names{1},'_',' ')]);
            ylabel(['RFP-' strrep(params.channel_names{2},'_',' ')]);
            zlabel(['YFP-' strrep(params.channel_names{3},'_',' ')]);
            view(axes1,[-37.5 30]);
            
        elseif length(params.channel_names) ==2
            scatter(green_channel,red_channel,double(size_droplet*0.5)); hold on;
            xlabel(['GFP-' strrep(params.channel_names{1},'_',' ')]);
            ylabel(['RFP-' strrep(params.channel_names{2},'_',' ')]);
        end
        % zlabel('640-DNA');
        
        grid(axes1,'on');
        % Set the remaining axes properties
        set(axes1,'FontSize',24);
        
        set(gcf, 'Position', get(0, 'Screensize'));
        
        saveas(gcf,output_file);
    end
end


variable_names = {};

for c=1:1:length(params.channel_names)
    variable_names = [variable_names, ['Partition_ratio_mean_' params.wavelength_channels{c} ],['Partition_ratio_std' params.wavelength_channels{c} ], [ 'Condensed_f_mean_' params.wavelength_channels{c}],[ 'Condensed_f_std_' params.wavelength_channels{c}],[ 'Total_raw_i_mean_' params.wavelength_channels{c}],[ 'Total_raw_i_std_' params.wavelength_channels{c}]];
end

save([output_folder_sheets 'example_params.mat'],'params');
Data = array2table(summary_statistics);
Data(:,6*length(params.wavelength_channels)+1) = [write_file_names];
Data.Properties.VariableNames = [variable_names, 'Image_file'];
if ~params.average_scaffold
    writetable(Data, summary_statistics_file,'Sheet',['SS_' num2str(params.minimum_drop_size) '_TI_' num2str(params.threshold_multiplier) '_scaf_' (params.wavelength_channels{params.scaffold_channel}) '_bg_' num2str(bg_value) ]);
else
    writetable(Data, summary_statistics_file,'Sheet',['SS_' num2str(params.minimum_drop_size) '_TI_' num2str(params.threshold_multiplier) '_scaf_avg_bg_' num2str(bg_value) ]);
end
end