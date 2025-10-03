%This is for another strategy of choosing RECA models.
%In this strategy, it simply tested all the possible combinations and
%choose the combinations with the least VIF as the results


function [VIFtable, VData]= RECA_VIF(data)
tic
        VData = table();

% Define VIF threshold
    vifThreshold = 10; % You can change this to 10 or another value based on your criteria
    %Get the inputs parameters
        X = [data(:,["thetao_A","so_A"]),...
             data(:,["o2_A","no3_A","po4_A","si_A","talk_A"]).*1000]; % Replace with your variable names
        % Get all binary combinations of the parameters
       

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
        VData    = X;
        VIFtable = join(string(X.Properties.VariableNames), ', ');

        parameters = X.Properties.VariableNames; % Get variable names from the table
        combos = nchoosek(parameters, 2);
        % Calculate and add each interaction term to the table
        for i = 1:size(combos, 1)
            % Generate the name for the new interaction variable
            interactionVarName = strjoin(combos(i, :), ':');
            % Calculate the interaction term by multiplying the two columns
            interactionColumn = X.(combos{i, 1}) .* X.(combos{i, 2});  
            % Add the new interaction column to the table
            X.(interactionVarName) = interactionColumn;
        end

           % Start the second iterative process
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
                    highVIF = find(VIF(5:end) > vifThreshold) + 4;        
                  if ~isempty(highVIF)
                    % Exclude the variable with the highest VIF
                     [~, maxIdx] = max(VIF(highVIF));
                     idxToRemove = highVIF(maxIdx);
                     X(:, idxToRemove) = [];
                 else
                    % If no variables exceed the threshold, exit the loop
                    continueRemoving = false;
                 end
            end

        % Save the input tables
        VData    = X;
        VIFtable = join(string(X.Properties.VariableNames), ', ');

   end    
