%% 


   %note VData contains key %Strategy 1
   for z= 1:numel(Year)

        tempyear = find(rdata.time==Year(z));
        X = FData(tempyear,:); % Replace with your variable names
        Y = rdata(tempyear,"dissic_A"); % Assuming 'DIC' is your dependent variable
        
        % Save the input tables
        InputXtable{z} = X.Properties.VariableNames;

        % Extract the names of the predictor variables dynamically
        predictorNames = X.Properties.VariableNames; % Exclude the response variable name if it's the first column
        responsorNames = Y.Properties.VariableNames{1};
        
        % Construct the formula string for stepwiselm, including all predictors without interactions
        mdl = stepwiselm([X,Y],'interactions', 'Upper', 'linear', 'Criterion', 'bic', 'Verbose', 0);

        MLRModelName = sprintf('%s%s_%s','USEC',string(floor(Year(z))),'MLR');
        % eval([MLRModelName ' = model;']);
        
       % Assuming you have models named model1, model2, model3, etc.
       %MLRmodels.USWC_Subsample = model; % Add your models to this array
        MLRModelST2{z} =   mdl;
        MLRmodelsST2.USWC_Subsample = MLRModelST2';
        MLRModelNamesST2{z} =   MLRModelName;
 
    % Extract statistics from each model
        R_squared(z) = MLRmodelsST2.USWC_Subsample{z}.Rsquared.Ordinary;
        Adj_R_squared(z) = MLRmodelsST2.USWC_Subsample{z}.Rsquared.Adjusted;
        RMSE(z) = MLRmodelsST2.USWC_Subsample{z}.RMSE;
        % Concatenate all TermNames into a single string, separated by commas
        Formula{z} = join(string(MLRmodelsST2.USWC_Subsample{z}.Formula.TermNames), ', ');
        MLRModelName = sprintf('%s%s%s_%s','USEC',string(floor(Year(z))),'MLR');
        MLRModelNamesST2{z} =   MLRModelName;
   
   end
   
       % % Create a summary table
       %  modelComparison.USWC_Subsample = table(R_squared', Adj_R_squared', RMSE', Formula',...
       %                      'RowNames', MLRModelNamesST2); % Add more row names if you have more models
       % 
       % 
        clear  allCombinationsFormulas models2 models2ref enR_squaredS2 enAdj_R_squaredS2 MLRModelNamesS2_Formula enR_squaredS2ref enAdj_R_squaredS2ref  MLRModelNamesS2ref_Formula
        clear  R_squaredS2 Adj_R_squaredS2 RMSE_S2 R_squaredS2ref Adj_R_squaredS2ref RMSE_S2ref MLRModelNamesS2

        
        %Merge all the terms into one cell and get ready for the uniform equation
       allTerms = {};
        for i = 1:length(Formula)
            % Split the current string into individual terms and trim whitespace
            terms = strtrim(strsplit(Formula{i}, ','));
            % Append these terms to the list of all terms
            allTerms = [allTerms, terms];
        end
        
        % Count the occurrence of each unique term
        [uniqueTerms, ~, idx] = unique(allTerms); % 'idx' are indices for reconstructing original rdata from 'uniqueTerms'
        frequencyCounts = accumarray(idx, 1); % Count occurrences
        
        % Ensure 'uniqueTerms' is a column vector to match the dimensions of 'frequencyCounts'
        uniqueTerms = uniqueTerms(:);
        
        % Now, 'uniqueTerms' and 'frequencyCounts' should have the same length
        % Combine the terms and their counts into a table for easy viewing/sorting
        termTable = table(uniqueTerms, frequencyCounts, 'VariableNames', {'Term', 'Frequency'});
        termsTable = termTable(2:end,:);
        % Assuming 'termsTable' is already defined and contains a 'Term' column and a 'Frequency' column

        % Filter 'termsTable' for rows where 'Frequency' is greater than 30 and extract the 'Term' column
        % Convert the table column to a cell array of character vectors if it's not already
        termsTable.Frequency(termsTable.Term == "so_A") = length(Year);
        termsTable.Frequency(termsTable.Term == "thetao_A") = length(Year);
        termsTable.Frequency(termsTable.Term == "o2_A") = length(Year);

        MostcommonTerms = termsTable.Frequency > length(Year)-1; 
        commonTerms = termsTable.Frequency > length(Year).*0.7 & termsTable.Frequency < length(Year) ;%If this is a reasonable criteria
                
        %Add an ensemble to the case
        termsToJoin = termsTable.Term(MostcommonTerms);
        termsToEnsemble = termsTable.Term(commonTerms);

        
        % Ensure termsToJoin is a cell array of character vectors for strjoin
        % If termsToJoin is a string array, you might need to convert it to a cell array
        if isstring(termsToJoin)
            termsToJoin = cellstr(termsToJoin); % Convert string array to cell array of character vectors
        end
        
                % Ensure termsToJoin is a cell array of character vectors for strjoin
        % If termsToJoin is a string array, you might need to convert it to a cell array
        if isstring(termsToEnsemble)
            termsToEnsemble = cellstr(termsToEnsemble); % Convert string array to cell array of character vectors
        end


        % Construct the formula
        allCombinationsFormulas = {};
        allCombinationsFormulas{1} = sprintf('%s ~ %s', responsorNames, strjoin(termsToJoin, ' + '));
        for comboSize = 1:length(termsToEnsemble)
             combinations = nchoosek(termsToEnsemble, comboSize);
             % Loop through each combination and construct the formula
             for i = 1:size(combinations, 1)
                 % Get the current combination and convert it into a string for the formula
                  currentCombination = strjoin(combinations(i, :), ' + ');
        
                 % Construct the formula by adding the combination to the existing most common terms
                  formula = sprintf('%s ~ %s + %s', responsorNames, strjoin(termsToJoin, ' + '), currentCombination);
        
                 % Store the formula in the cell array
                  allCombinationsFormulas{end+1} = formula;
             end
        end


        %concentrat X and Y

        % Perform stepwise linear regression using stepwiselm, specifying 'linear' to avoid interactions
        %mdl = stepwiselm(horzcat(Y,X), formula, 'Upper', 'linear', 'Criterion', 'bic', 'Verbose', 2);'

        MLRModelS2AIC = cell(length(Year), length(allCombinationsFormulas));
        MLRModelS2AICref = cell(length(Year), length(allCombinationsFormulas));
        MLRModelS2R2 = cell(length(Year), length(allCombinationsFormulas));
        MLRModelS2RMSE = cell(length(Year), length(allCombinationsFormulas));
        
   for z= 1:numel(Year)

        Year = table2array(unique(round(rdata(:,4),1)));  % Drop the decimal part and duplicates and make it a array
        tempyear = find(rdata.time==Year(z));
        refyear  = find(rdata.time==Year(1));

        %[VIFtable, FData] = RECA_DataVIF(data);
        data = [rdata(:,["thetao_A","so_A"]),...
             rdata(:,["o2_A","po4_A","si_A","talk_A","dissic_A"])]; % Replace with your variable names

        X = data(tempyear,1:(end-1)); % Replace with your variable names
        Y = data(tempyear,end); % Assuming 'DIC' is your dependent variable
        Xref = data(refyear,1:(end-1)); % Replace with your variable names
        Yref = data(refyear,end); % Assuming 'DIC' is your dependent variable
       
        % Save the input tables
        InputXtable{z} = X.Properties.VariableNames;

        % Extract the names of the predzictor variables dynamically
        predictorNames = X.Properties.VariableNames; % Exclude the response variable name if it's the first column
        responsorNames = Y.Properties.VariableNames{1};
   


        for m = 1: length(allCombinationsFormulas)
            models2 = fitlm(horzcat(X,Y), allCombinationsFormulas{m},'RobustOpts','on');
            models2ref = fitlm(horzcat(Xref,Yref), allCombinationsFormulas{m},'RobustOpts','on');
            % eval([MLRModelName ' = model;']);
           MLRModelS2AIC{z,m} = models2.ModelCriterion.AICc;
           MLRModelS2AICref{z,m} = models2ref.ModelCriterion.AICc;
           MLRModelS2R2{z,m} = models2.Rsquared.Adjusted;
           MLRModelS2RMSE{z,m} = models2.RMSE;
        
           % Assuming you have models named model1, model2, model3, etc.
           %MLRmodels.USWC_Subsample = model; % Add your models to this array
            MLRModelS2{z,m}     =   models2;
            MLRModelS2ref{z,m}  =   models2ref;
        end
   end
        MLRModelS2R = cell2mat(MLRModelS2R2); % Assumes all elements are numeric, adjust if necessary
        colsToDelete = any(MLRModelS2R < 0.8, 1);
        MLRModelS2(:, colsToDelete) = [];
        MLRModelS2ref(:, colsToDelete) = [];
        MLRModelS2R2(:, colsToDelete) = [];


        MLRmodelsS2.USWC_Subsample = MLRModelS2;
        MLRmodelsS2ref.USWC_Subsample = MLRModelS2ref;

        fieldName = sprintf('Subregion%d', n);
        fieldNameref = sprintf('Subregion%d_ref', n);

        MLRRECAmodels.(fieldName) = MLRModelS2';
        MLRRECAmodels.(fieldNameref) = MLRModelS2ref';
        
      clear modelComparisonSum 


    for z= 1:numel(Year)

        for m = 1: size(MLRModelS2,2)
 
        R_squaredS2(z,m) = MLRModelS2{z,m}.Rsquared.Ordinary;
        R_squaredS2ref(z,m) = MLRModelS2ref{z,m}.Rsquared.Ordinary;

        Adj_R_squaredS2(z,m) = MLRModelS2{z,m}.Rsquared.Adjusted;
        Adj_R_squaredS2ref(z,m) = MLRModelS2ref{z,m}.Rsquared.Adjusted;

        RMSE_S2(z,m) = MLRModelS2{z,m}.RMSE;
        RMSE_S2ref(z,m) = MLRModelS2ref{z,m}.RMSE;

        AIC_S2(z,m) = MLRModelS2{z,m}.ModelCriterion.AIC;

        [h, p] = lillietest(MLRModelS2{z,m}.Residuals.Raw);
        MLRModelS2lillietesth(z,m) = h; %h = 0, it's normal distribution
        MLRModelS2lillietestp(z,m) = p; %p >0.05; it's normal distribution


        % Concatenate all TermNames into a single string, separated by commas
        MLRModelNamesS2_Formula{z,m} = join(string(MLRmodelsS2.USWC_Subsample{z,m}.Formula.TermNames), ', ');
        MLRModelNamesS2ref_Formula{z,m} = join(string(MLRmodelsS2ref.USWC_Subsample{z,m}.Formula.TermNames), ', ');
        MLRModelName = sprintf('%s%s%s_%s','USEC',string(floor(Year(z))),'MLR');
        MLRModelNamesS2{z} =   MLRModelName;
       
         %calculate an ensemble mean
         enR_squaredS2     = mean(R_squaredS2, 2);
         enAdj_R_squaredS2 = mean(Adj_R_squaredS2, 2);
         enRMSE_S2         = mean(RMSE_S2,2);

         enR_squaredS2ref     = mean(R_squaredS2ref, 2);
         enAdj_R_squaredS2ref = mean(Adj_R_squaredS2ref, 2);
         enRMSE_S2ref         = mean(RMSE_S2ref,2);

         % Create a summary table
         modelComparisonSum.Subsample = table(enR_squaredS2, enAdj_R_squaredS2, enRMSE_S2,MLRModelNamesS2_Formula(:,1:m), ...
                                                   enR_squaredS2ref, enAdj_R_squaredS2ref, enRMSE_S2ref,MLRModelNamesS2ref_Formula(:,1:m),'RowNames', MLRModelNamesS2); % Add more row names if you have more models
         end
    
    end

   for z= 1:numel(Year)

        %[VIFtable, FData] = RECA_DataVIF(data);
        data = [rdata(:,["thetao_A","so_A"]),...
             rdata(:,["o2_A","po4_A","si_A","talk_A","dissic_A"])]; % Replace with your variable names

        %[VIFtable, FData] = RECA_DataVIF(data);
         tempyear = rdata.time==Year(z);
         FData_temp = data(tempyear,:); % Replace with your variable name
         VData_temp = rdata(tempyear,:); % Replace with your variable names
         % prepare for prediction
         %FData_temp = table2array(FData_temp);
         %VData_temp = table2array(VData_temp)z;

          for m = 1: size(MLRModelS2,2)
                 %Define the variables must be included in the model
                TempPredicts(:,m) = predict(MLRModelS2{z,m},FData_temp);
                TempPredict1(:,m) = predict(MLRModelS2ref{z,m},FData_temp);
                dCanthro_RECAtemp(:,m)   =  TempPredicts(:,m) - TempPredict1(:,m) ; 
          end
             
        dCanthro_RECA{:,z}  = mean(dCanthro_RECAtemp,2);
        dCanthro_stdRECA{:,z}  = std(dCanthro_RECAtemp,0,2);

        %Define the variables must be included in the model
        %TempPredicts = predict(MLRmodels.USWC_Subsample{z},VData_temp);
        %TempPredict1 = predict(MLRmodels.USWC_Subsample{1},VData_temp);
        %dCanthro_RECAtemp(:,z)   =  TempPredicts - TempPredict1 ; 
        %check
        % scatter3(VData_temp.lon,VData_temp.lat,-VData_temp.depth,30,TempPredicts(:,1)-VData_temp.disinc_A,'filled')
        % clim([-10 10])
        % colorbar
        % % % % 
        % % % % 
        % scatter3(VData_temp.lon,VData_temp.lat,-VData_temp.depth,30,dCanthro_RECAtemp(:,1),'filled')
        % colorbar
        % clim([0 20])

        clear FData_temp dCanthro_RECAtemp termsTable TempPredicts TempPredict1
  end 

  clear data
toc

