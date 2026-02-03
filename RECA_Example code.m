clear all
clc

%This is an example code to show how to use RECA anthropogenic carbon
%algorithms. Any question/bug report please email xinyuli@ude.edu

%% load data and organize for calculations

ECOA_CODAP = readtable("CODAP_ECOA.xlsx"); %replace your own data

% Initialize structure to store each cruise's data
CruiseData = struct();

% List of required columns (now including Latitude and Longitude)
cols = { ...
    'EXPOCODE', 'Cruise_ID', 'Year_UTC', 'Month_UTC', 'Depth', ...
    'Latitude', 'Longitude', ...
    'CTDTEMP_ITS90', 'CTDTEMP_ITS90_flag', ...
    'CTDSAL_PSS78', 'CTDSAL_PSS78_flag', ...
    'Oxygen', 'Oxygen_flag_1', ...
    'DIC', 'DIC_Flag', ...
    'TA', 'TA_flag', ...
    'Phosphate', 'Phosphate_flag', ...
    'Silicate', 'Silicate_flag', ...
    };

% Read only the selected columns
opts = detectImportOptions(file, 'Sheet', 'CTD DATA');
opts.SelectedVariableNames = cols; clear cols
rdata = readtable(file, opts);

% Get unique Cruise IDs
cruise_list = unique(rdata.Cruise_ID);

for i = 1:length(cruise_list)
    cruise_name = cruise_list{i};
    idx = strcmp(rdata.Cruise_ID, cruise_name);
    cruise_data = rdata(idx, :);
    safe_name = matlab.lang.makeValidName(cruise_name); % Ensures valid field name
    CruiseData.(safe_name) = cruise_data;
end
clear rdata

 %% Filter data for use:
% List of variable/flag pairs
% Get all cruise names
cruiseNames = fieldnames(CruiseData);

for c = 1:length(cruiseNames)
    tbl = CruiseData.(cruiseNames{c});
    % Get all variable names in the current table
    colNames = tbl.Properties.VariableNames;
    
    % Find all columns whose names end with '_flag' or '_Flag' (case insensitive)
    isFlagCol = cellfun(@(x) ~isempty(regexp(x, '_flag$', 'ignorecase')), colNames);
    flag_cols = colNames(isFlagCol);
    
    % Only proceed if there are flag columns!
    if ~isempty(flag_cols)
        % Build a logical matrix (rows x flags)
        flag_matrix = false(height(tbl), length(flag_cols));
        for f = 1:length(flag_cols)
            col = tbl.(flag_cols{f});
            % Only flags 2 or 6 are good
            flag_matrix(:,f) = (col == 2 | col == 6);
        end
        
        % Keep only rows where all flags pass
        mask = all(flag_matrix, 2);
        CruiseData.(cruiseNames{c}) = tbl(mask, :);
    end
end

%% Match the cruise data to the function name
% Get all cruise names and Vertically concatenate all tables (all columns should match now)

cruiseNames = fieldnames(CruiseData);
tbls = cellfun(@(x) CruiseData.(x), cruiseNames, 'UniformOutput', false);
USEC_flat = vertcat(tbls{:});
oldVars = {'Longitude', 'Latitude', 'Depth', 'Year_UTC', 'DIC', 'CTDOXY', 'Phosphate', ...
           'Silicate', 'CTDSAL_PSS78', 'TALK', 'CTDTEMP_ITS90', 'DIC'}; % Example
USEC_All = sortrows(USEC_All, 'Year_UTC');
newVariableNames = {'lon', 'lat', 'depth', 'time', 'disinc_A','o2_A','po4_A','si_A','so_A','talk_A','thetao_A','dissic_A'}; %must be remaed by these name to use RECA
USEC_All.Properties.VariableNames = newVariableNames;

%% Subregion for each coastal region (This is not always necessary, depending on the coasts)
%Subset region
Subregion_Criteria{1} = [35,41.5]; %Not include the gulf of maine region because no data max(USEC_All.lon)
Subregion_Direction{1} = 'lat';
SData = LME_subregion(USEC_All, Subregion_Criteria{1}, Subregion_Direction{1});

%% Start RECA calculation
RECA_struct = struct();
AllerrorResults = [];
ComRegionResults = [];

% Loop over each subregion in SData
for n = 1:length(SData)
    
    tdCanthro_RECA = [];
    % Extract the subsample table for the current subregion
     tdata = SData{n};
            
    % Extract longitude and latitude vectors from the table
     Lons = tdata.lon;
     Lats = tdata.lat;
            
    %Data below mixed layer depth (or cetertain depth)      
     
    rdata = tdata((tdata.depth>20&tdata.depth<2000),:);
        
     FData = table();
     % Define VIF threshold
     vifThreshold = 10; % You can change this to 10 or another value based on your criteria
     %Get the inputs parameters
     X = [rdata(:,["thetao_A","so_A"]),...
          rdata(:,["o2_A","po4_A","si_A"])]; % Replace with your variable names
       
     % Start the iterative process
        continueRemoving = true;
        
            while continueRemoving
                numVars = width(X);
                VIF = zeros(1, numVars);
        
                   % Calculate VIF for each variable
                   for k = 1:numVars
                        tempX = X(:, setdiff(1:numVars, k));
                        [~,~,~,~,stats] = regress(X{:,k}, [ones(size(tempX, 1), 1) tempX{:,:}]);
                        R2 = stats(1);
                        VIF(k) = 1 / (1 - R2);
                    end
        
                % Find variables with high VIF
                  highVIF = VIF > vifThreshold;
        
                 if any(highVIF)
                    % Exclude the variable with the highest VIF
                    [~, idx] = max(VIF);
                    X(:, idx) = [];
                else
                    % If no variables exceed the threshold, exit the loop
                    continueRemoving = false;
                end
            end
        % Save the input tables
        FData    = X;
        VIFtable = join(string(X.Properties.VariableNames), ', ');
    
    Year = table2array(unique(round(rdata(:,4),1)));  % Drop the decimal part and duplicates and make it a array
    
     run RECAInter.m

     %Organizing output
     for z = 1:length(dCanthro_RECA)
         tdCanthro_RECA = [tdCanthro_RECA; dCanthro_RECA{z}];
     end
             
     dCRECA = array2table(tdCanthro_RECA);
     dCantdata = [rdata,dCRECA];
     ComRegionResults = [ComRegionResults;dCantdata];
            
     % Clean up temporary variables for this iteration
     clear dCanthro_RECA modeltops
     toc
     clear MLRModelS2 dCanthro_RECA dCanthro_RECAtemp dCanthro_stdRECA MLRModels2 MLRModelST2
     clear MLRModelNamesST2 MLRModelS2R2 MLRModelS2ref MLRModelS2R MLRModelS2AIC MLRModelS2AICref
end  % end for n
         
%Save as RECA Table
RECATable    = ComRegionResults;
% Save the RECA struct to the base workspace with a dynamic name,
assignin('base', ['ECOA_RECA'], RECATable);
USEC_CODAP.(['dCanthro_RECA']) = RECATable;

clear fulldata dCanthro_RECA dCRECA dCantdata

filename = 'USEC_CODAP.mat';
save(filename, 'USEC_CODAP');

