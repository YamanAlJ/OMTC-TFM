%AVERAGED_OMTC_TFM: computes and plots the averaged TFM-OMTC across one CONDITION

%USER INPUT 1/9:OMTC data array name files (lines 5-27)and TFM data name
%files (lines 31-53)
clc; close all; 
clear all;

ARRAY = [ "LnCaP_10nMR1881_3kPa/", ...
           "PC3AR_10nMR1881_3kPa/"]
                   
%TFM Arrays             
%TFM_ARRAY = ["LNCaP_1nMR1881_12kPa/", ...
%              "PC3_1nMR1881_12kPa/", ...
%              "PC3AR_1nMR1881_12kPa/"]
                             

               

       
% USER INPUT 2/9: selects which OMTC and TFM arrays to plot
MEGA_OMTC= {ARRAY}
MEGA_TFM = {ARRAY}
               

%selection of OMTC attributes to plot
a = 3
h = 5
GStar = 6                 
G1 = 9
G2 = 10
tan_delta = 11
ChooseMean = 12
ChooseMedian = 13

%USER INPUT 3/9: limit on OMTC attribute
G_limit = 100000


% USER INPUT 4/9: choose OMTC attribute to plot:
OMTC_ATTRIBUTE = GStar


% USER INPUT 5/9: TFM attribute: RMST 1, strain energy 2, avergade displacment 3
TFM_ATTRIBUTE = 2

%USER INPUT 6/9: input 12 for the mean OMTC attribute, 13 for the median OMTC aaattribute
%IMPORTANT: change line 179 function to median or mean 
MedianOrMean = 12

switch MedianOrMean
    case ChooseMean
        mean_title = "Mean"
    case ChooseMedian
        mean_title = "Median"
end

%switches plot setting parameters based on chosen attribute
switch OMTC_ATTRIBUTE
     case a
        edges= [0:0.1:2];
        attribute_title ='AMPLITUDE (um)'
        y_scale = [0.0,2];
        
     case h
        ylim= [0:0.01:0.25];
        attribute_title ='H VALUE (s)';
        y_scale = [0.0,0.25];
    case GStar
        edges = [0:25:1000];
        attribute_title ='Apparent Modulus (Pa)';
        y_scale = [0.0,1000];
    case G1
        edges = [0:25:1200];
        attribute_title ='Storage Modulus ';
        y_scale = [0.0,1200];
    case G2
        edges = [0:25:1200];
        attribute_title ='Loss Modulus ';
        y_scale = [0.0,1200];
    case tan_delta
        edges = [0:0.1:2.5];
        attribute_title ='Tan Delta'
        y_scale = [0.0,2.5];
end

%USER INPUT 7/9: directory name where OMTC results are located
CONDITION_OMTC_FINAL = 'D:\RAW DATA\Compiled_OMTC'
%USER INPUT 8/9: directory name where TFM results are located
CONDITION_TFM_FINAL = 'D:\RAW DATA\TFM_results'
FILE_PATTERN = "*_SUMMARY_STATS_FINAL_GM.xlsx"

f1 = figure();
Ax(1) = axes(f1); 

colors = [1.0 0.0 0.0

    0.0 1.0 0.0
    
    0.0 1.0 1.0

    0.47 0.25 0.80
     ];  



%% AVERAGES STRAIN ENERGY AND CORRESPONDING OMTC ATTRIBUTE
for j= 1:length(MEGA_OMTC)
    
    label_array = []
    OMTC_array = []
    mean_strain_array =[]
    condition_array= []
    err_array= []
    err_strain_array=[]
    filtered_attribute_column = []
    strain_column=[]

    OMTC_column = MEGA_OMTC{j}

    % averages OMTC data for given cell type on given modulus
    for i = 1:length(OMTC_column)
            

            stiffness = OMTC_column{i};

            
            file_search_pattern_OMTC = strcat(CONDITION_OMTC_FINAL ,"/", stiffness, "*_FINAL/", FILE_PATTERN);
            file_directory = dir(file_search_pattern_OMTC);

            if (length(file_directory) ~= 1) 
                error ('Expected single excel file in ' + file_search_pattern_OMTC +" number:"+length(file_directory));
            end

            file_string = strcat(file_directory(1).folder,"/", file_directory(1).name);
            condition_table_NaN = readtable(file_string);

            condition_table = rmmissing(condition_table_NaN);

            attribute_column = table2array(condition_table(:, OMTC_ATTRIBUTE))
            for k = 1:length(attribute_column)
                if attribute_column(k) < G_limit;
                    filtered_attribute_column =   [filtered_attribute_column; attribute_column(k)] ;
                end
            end
            %USER INPUT 9/9: change function to median or mean 
            mean_G = mean(filtered_attribute_column);
            stderror = std(filtered_attribute_column) / sqrt( length(filtered_attribute_column));
            
            
             err_array = [err_array, (stderror/2)];
             OMTC_array = [OMTC_array; mean_G'];
               

    end

    TFM_ARRAY = MEGA_TFM{j}
    
    % averages corresponding TFM data for given cell type on given modulus

    for i = 1:length(TFM_ARRAY)

            
        TFM_folder = TFM_ARRAY{i};
        file_search_pattern_TFM = strcat(CONDITION_TFM_FINAL, "/", TFM_folder,"results.mat");
        file_directory = dir(file_search_pattern_TFM);

        if (length(file_directory) ~= 1) 
            error ('Expected single excel file in ' + file_search_pattern_TFM +" number:"+length(file_directory));
        end

        file_string = strcat(file_directory(1).folder,"/", file_directory(1).name);
        results_file = load(file_string);

        switch TFM_ATTRIBUTE
         case 1
            TFM_data = (results_file.rmst)'
            TFM_title ='RMST (Pa)'
            factor_variable = 1
         case 2
            TFM_data = (results_file.strain_energy)'
            TFM_title ='Strain Energy (pJ)'
            factor_variable = 1000000000000
         case 3
            TFM_data = (results_file.average_displacement)'
            TFM_title ='Average Displacement (um)'
            factor_variable = 1

        end


        strain_column = [strain_column ; TFM_data]

            
        factor_correction = repmat(factor_variable, length(strain_column), 1)
        strain_to_scale = (strain_column.*factor_correction)    
        
 
        mean_strain = mean(strain_to_scale);
        stderror = std(strain_to_scale) / sqrt( length(strain_to_scale));
        err_strain_array = [err_strain_array, (stderror/2)];
        mean_strain_array = [mean_strain_array; mean_strain'];


             
 
    end
%% PLOTS STRAIN ENERGY VS OMTC ATTRIBUTE 
    color = [ 0.0 0.0 0.0 ]
    labels = {'LnCaP'; 'PC3AR'};         
    markers='ospd^'
   
    
    scatterplot = gscatter(mean_strain_array,  OMTC_array, labels, color,markers, 4);
    
    %to fill markers
%     for g = 1:length(scatterplot)
%         scatterplot(g).MarkerFaceColor = scatterplot(g).Color;
%     end
   
    for n = 1:length(scatterplot)
        set(scatterplot(n), 'MarkerFaceColor', colors(n,:));
        set(scatterplot(n), 'MarkerSize', 15);
    end
    
    
    hold on
    lineplot = plot(mean_strain_array,OMTC_array, 'lineWidth', 2, 'Parent', Ax(1))
    lineplot.Color= colors(j,:)
 
    
    
    set(Ax(1), 'Box','off')
    errorplot = errorbar(mean_strain_array,OMTC_array,err_array,err_array,err_strain_array,err_strain_array,'LineStyle','none', 'Color', 'k','Parent',Ax(1))
    
    leg1= legend(Ax(1),[scatterplot],labels, 'Location', "southeast");set(leg1,'FontSize',12);
    title (leg1, strcat('Modulus ', mean_title) );
    leg1.TextColor= 'k'
    %Ax(2) = copyobj(Ax(1),gcf);
    %delete(Ax(2),'Children');
    
    %mlabel1= plot (1000,1,'_', 'MarkerSize',100,'Parent',Ax(2))
    %mlabel1.Color = [ 0.90 0.14 0.14]
    %mlabel2= plot (1000,1,'_','MarkerSize',100,'Parent',Ax(2))
    %mlabel2.Color = [ 1.00 0.54 0.00]
    %mlabel3= plot (1000,1,'_','MarkerSize',100,'Parent',Ax(2))
    %mlabel3.Color = [0.25 0.80 0.54]
    %mlabel4= plot (1000,1,'_','MarkerSize',100,'Parent',Ax(2))
    %mlabel4.Color = [ 0.47 0.25 0.80]
    
    %set(Ax(2), 'Color', 'none', 'XTick', [], 'YAxisLocation', 'right', 'Box', 'Off', 'Visible', 'off')

    %leg2 = legend([mlabel1,mlabel2,mlabel3,mlabel4],'22RV1', 'LNCaP', 'DU145', 'PC3', 'Location', "southwest");
    %set (leg2, "FontSize", 12, "LineWidth", 1.5);
    %title (leg2, 'Cell Type' );

    xlabel(TFM_title)
    set (gca, "FontSize", 12, "LineWidth", 1.5);
    ylabel(attribute_title)
    
    set (gca, "FontSize", 12, "LineWidth", 1.5);
    
    title( strcat(TFM_title," v.s."," ", mean_title, " ", attribute_title));
    set (gca, "FontSize", 15, "LineWidth", 1.5);
    
    xlim([0, 2])
    ylim([0,1000])


end 


%% TO OVERLAY IMAGE
%    mini=imread('mean_mini.tif'); axes('position',[.4 .45 .45 .45]); imagesc(mini);
%     
     hold on
