% Goal:  calculate the chemotactic index for each cell in a microchamber video
% Written by Lauren Hein and Kristen Loesel for organizing/ analyzing
% underagarose chemotaxis tracks - October 2021 // 2022
% adjusted and rewritten by Sam Collie for analyzing chemotaxis tracks from 
% microchmber or underagarose assays - February 2023
% MEASURES: Speed, XFMI, Accumulated Distance, Euclidean Distance, 
% Directionality and Angle of Deviation alongside X and Y
% FILTERS: Filters cells which travel less than 20 microns as these
% are mostly broken tracks or background. Additionally removes data for
% cells past 60 min of movie time since the gradient is often dispersed at
% this stage. Both filters can be adjusted.
% NOTES: can be adjusted for microchamber analysis by adjusting XFMI to
% YFMI. also spits out a plot of the xy of each tracked cell
%

% clear the workspace and the command window
clear
clc

close all;
clc;
clearvars;
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Saving Directory and File Output Names %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
% 
% The following code finds all the data files in the folder specified with 
% the ending name pattern of '*spots.csv'
% Set the starting folder where your files are stored 
prompt1 = 'Please indicate the folder where your data is stored.' %#ok<NOPTS> 
startingFolder = uigetdir();

% Put out an error message if the user cancels
if startingFolder == 0
    errorMessage = sprintf('Please specify a new folder.', startingFolder);
    uiwait(warndlg(errorMessage));
    startingFolder = uigetdir(); % Ask for a new one.
    if startingFolder == 0
         % User clicked Cancel
         return;
    end
end


% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(startingFolder, '*spots.csv'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
%Loop through each text file in the folder
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
% Remove characters after 'um' in file name. This will be used to create
% an additional column containing the name of the condition. This needs to
% be adjusted based on the names of your files
conditionName = eraseBetween(baseFileName, "um", ".csv");
condition = extractBefore(conditionName, '.csv');
% Get the file name without the extension
[~, baseNameNoExt, ~] = fileparts(baseFileName);
%reading your text file and saving the matrix
data = readmatrix(fullFileName);
%this is the pattern that you want to save the text files with 
baseFileName = sprintf('%s_results.csv', baseNameNoExt);
sumFileName = sprintf('%s_results_mean_med_std.csv', baseNameNoExt);
%this is where the text files will be exported to (your starting folder)
organized_dataset_save = fullfile(theFiles(k).folder, baseFileName);
%Getting rid of all the data columns we do not need 
data(:,[1:2,4,7,9:20])=[];
data(1,:)=[];
%reorganizing the columns so that column 2 = unique cell id; column 3 =
% frames; column 4 = x values; and column 5 = y values
data = data(:,[1 4 2 3]);
% remove rows that contain NaN/missing values
data = rmmissing(data);
% remove rows of cells for which there is only one data point. We do this  
% because this causes a bug with time interval calculations
for k = 1 : (length(data)-1)
    c_id = data(k,1);
    one_check = data((data(:,1)==c_id),1);
    if length(one_check) < 2
        data(k,:) = [];
    end
end
% sort rows based on column 2 (cellID)
data=sortrows(data,1);

% after sorting, the difference between the first two data points in the
% time (sec) column is the interval in sec which we convert to min
time_interval = abs(data(1,2) - data(2, 2))/60;

% filter out cells who travel further than specified time (60min or 3600sec)
data_filter = data(:,2) > 3600;
data(data_filter,:) = [];

% number of tracks/cell IDs in file (must add 1 because initial track index is 0)
n = max(data(:, 1)) + 1;
% number of rows in file
l = length(data);

% create outputs matrix of to hold 
% the 6 different variables. this will be a placeholder with an extra row 
% but we will remove this later
outputs = zeros(1,7);
output_con = string(zeros(l,1));
tc = 1;

% Randomly assign the 100 cells that will be plotted for each
% movie. this formula generates k random integers in the interval 
% 1 -> n. +1 for the bounds since 0 is a cell ID but not an integer
n = max(data(:,1))+1;
k = 100;
r = randperm(n,k);

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculating Chemotactic Measures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\

%%%%%%%%% this is a series of nested loops that take all the info for one
%%%%%%%%% cell id, determine if the cell has traveled for longer than the 
%%%%%%%%% desired threshold, normalize it to (0,0), plots the data, and
%%%%%%%%% calculates the XFMI, accumulated distance, speed, and Euclidian
%%%%%%%%% distance
% loop through all track number 
for i = 0:(n-1)
    % set/reset group variable to 0
    group = 0;
    % set/reset adjusted position to 0
    adjusted_pos = 0;
    % set/reset counter to 1
    cc = 1;
    
    % run a for loop to cycle through every row of the data matrix, which
    % contains the data sorted by cell id
    for j = 1:l
        
        % determine if each row belongs to the current cell ID (i)
        if data(j, 1) == i
            % assign cell id and track slice to new temporary matrix
            group(cc, 1) = data(j, 1);
            group(cc, 2) = data(j, 2);
            % assign x value of track
            group(cc, 3) = data(j, 3);
            % assign y value of track
            group(cc, 4) = data(j, 4);
            % update counter variable in order to save next data point to a
            % new row of the matrix (group)
            cc = cc + 1;
        end
        
    end
    
    % determine the dimensions of the group matrix
    % the size command returns (rows, columns) of a matrix
    s = size(group);
    % the number of rows (length) is equal to the first value of (size)
    lg = s(1);
    
            % calculate the accumulated distance the cell traveled by summing 
            % the distance between each data point
            % set/reset accumulated distance variable to 0
            acc_dist = 0;
            % loop through the adjusted position matrix for this cell ID
            for z = 1:(lg-1)
                % assign x1 and y1
                x1 = group(z, 3);
                y1 = group(z, 4);
                % assign the next point as x2 and y2
                x2 = group((z+1), 3);
                y2 = group((z+1), 4);
                % use the distance formula to calculate the distance between
                % two points
                dist = sqrt((x2 - x1)^2 + (y2 - y1)^2);
                % add the distance calculated to the accumulated distance
                % variable
                acc_dist = acc_dist + dist;
            end

    % filter out tracks in which cells begin past the point of constriction
    % and tracks who travel less than 20 microns
    % if the 'group' matrix contains data for a cell for which the first y
    % value is less than or equal to 100um, then the code will continue
    % through this if statement and calculate the metrics of interest.  If it
    % does not meet the distance criteria, it will exit this if statement
    % and go back to line 71
    if acc_dist > 20
            % create a new matrix to hold the values of the adjusted position
            % of the cells
            adjusted_pos = zeros(lg, 4);
            % assign the cell id and slice number to the adjusted position
            % matrix
            adjusted_pos(:, 1:2) = group(:, 1:2);
            % define the amount to adjust the x position by using the x
            % position from the the first slice of the cell
            adjustx = group(1, 3);
            % define the amount to adjust the y position by using the y
            % position from the the first slice of the cell
            adjusty = group(1, 4);
            
            % loop through all slices of the cell's data and subtract the
            % 'adjustx' value from the x position data and the 'adjusty' value
            % from the y position data (results in the origin being (0,0) and
            % everything else being adjusted accordingly)
            for ap = 1:lg
                % x position calculation, assign new value to the adjusted
                % position matrix
                adjusted_pos(ap, 3) = group(ap, 3) - adjustx;
                % y position calculation, assign new value to the adjusted
                % position matrix
                adjusted_pos(ap, 4) = group(ap, 4) - adjusty;
            end
            
            % add adjusted position matrix for this particular track to a new
            % combined data vector
            data = vertcat(data, adjusted_pos);
            
            % plot data that has been rescaled to start at the origin for 
            % 100 random cells in the dataset
            if ismember(group(1,1)+1, r)
                X = adjusted_pos(:,3); % x positions of cell
                Y = adjusted_pos(:,4); % y positions of cell
                Z = adjusted_pos(:,2)/59.9; % time of cell position in minutes
                surf([X(:) X(:)], [Y(:) Y(:)], [Z(:) Z(:)], ...  % Reshape and replicate data
                'FaceColor', 'none', ...    % Don't bother filling faces with color
                'EdgeColor', 'interp', ...  % Use interpolated color for edges
                'LineWidth', 0.85);            % adjust line thickness
                view(2);   % Default 2-D view
                colorbar;  % Add a colorbar
                axis([-200 1000 -500 500]);
                hold on
            end
        
        %%%%%%%%%% calculate XFMI for each cell
        % first calculate the accumulated distance by summing the distance
        % between each data point
        % set/reset accumulated distance variable to 0
        acc_dist = 0;
        % loop through the adjusted position matrix for this cell ID
        for z = 1:(lg-1)
            % assign x1 and y1
            x1 = adjusted_pos(z, 3);
            y1 = adjusted_pos(z, 4);
            % assign the next point as x2 and y2
            x2 = adjusted_pos((z+1), 3);
            y2 = adjusted_pos((z+1), 4);
            % use the distance formula to calculate the distance between
            % two points
            dist = sqrt((x2 - x1)^2 + (y2 - y1)^2);
            % add the distance calculated to the accumulated distance
            % variable
            acc_dist = acc_dist + dist;
        end
        
        % then calcualte XFMI using accumulated distance and final x
        % position
        x_final = adjusted_pos(lg, 3);
        XFMI = x_final / acc_dist;
            % if the cell does not travel, the accumulated distance will be
            % equal to 0, resulting in an XFMI of NaN.  If this is the
            % case, change the XFMI to zero because the cell did not
            % chemotax in any direction
            if isnan(XFMI) == 1
                % change value to 0
                XFMI = 0;
            end
            
        % calculate total time the cell traveled
        % determine the first slice in which the cell appeared (first row,
        % second column of "group" matrix)
        slice_initial = group(1, 2);
        % determine the last slice in which the cell appeared (first row,
        % last column of "group" matrix)
        slice_final = group(lg, 2);
        % calculate the difference between the first and last slices
        slice_lapse = slice_final - slice_initial;

        % time the cell traveled in minutes
        time_lapse = slice_lapse/60;

        % calculate average speed in um/min (total time the cell traveled
        % is calculated earlier in the code to determine which tracks to
        % include based on time of track; accumulated distance is also
        % calculated earlier)
        speed = (acc_dist / time_lapse);
        
        % calculate the Euclidian distance (initial position is 0
        % because all data points have been normalized to origin)
        y_final = adjusted_pos(lg, 4);
        % use the distance formula to calculate the distance between the
        % final and initial points of the cell
        Euc_dist = sqrt((x_final - 0)^2 + (y_final - 0)^2);
        
        % calculate Directedness based on accumulated and euclidean
        % distance for each cell
        direct = (Euc_dist/acc_dist);

        % calculate Angle of Deviation from the x-axis for each cell
        ang_dev = atan2d(y_final,x_final);

        % add cell number, XFMI, speed, and accumulated distance
        % to temporary cell output matrix
        cell_outputs = [adjusted_pos(1, 1) XFMI speed Euc_dist acc_dist direct ang_dev];
        
        % append cell track number, XFMI, speed, and accumulated distance
        % to outputs summary matrix
        outputs = vertcat(outputs, cell_outputs);
        
    end
    % repeat back to line 71 to continue loop for next cell ID/track number
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculating Mean, Median, STD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\

outputs_val = outputs;
outputs_val(:,1) = [];
output_mean = mean(outputs_val);
output_median = median(outputs_val);
output_std = std(outputs_val);
outputs_sum = [output_mean 
    output_median 
    output_std];

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Creating Plot and Saving Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\

%%%%% format plot

% axis labels
xlabel('distance (microns)')
ylabel('distance (microns)')


% % specify the starting folder as the output folder
output_file = append(startingFolder,'\', baseFileName); 
output_summary_file = append(startingFolder,'\',sumFileName);
output_fig1 = append(startingFolder,'\', condition);
output_fig2 = append(output_fig1,'.eps');


% save figure and export as eps
saveas(gcf, output_fig1, 'png');
exportgraphics(gcf,output_fig2, 'ContentType','vector');
close all;

% remove placeholder row
outputs(1,:) = [];

% original number of tracks
original_tracks = max(data(:, 1))
% included number of tracks
final_size = size(outputs);
final_tracks = final_size(1)


% export the outputs table to an Excel spreadsheet
% add labels to the the columns in the outputs file
% do the same for the summary stats table

outputs_labeled = array2table(outputs, 'VariableNames', {'cell ID',  ...
    'XFMI', 'Speed (um/min)', 'Euclidean Distance', 'Accumulated Distance', ...
    'Directionality', 'Angle of Deviation'});
outputs_labeled.Condition = repmat(string(condition), final_tracks, 1);

outputs_summary_labeled = array2table(outputs_sum, 'VariableNames', {
   'XFMI', 'Speed (um/min)', 'Euclidean Distance', 'Accumulated Distance', ...
    'Directionality', 'Angle of Deviation'});

% write the data to a csv file
writetable(outputs_labeled, output_file);
writetable(outputs_summary_labeled, output_summary_file);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

