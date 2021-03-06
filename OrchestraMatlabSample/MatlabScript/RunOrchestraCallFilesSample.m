function RunOrchestraCallFilesSample()

    % First we have to add the directory containing OrchestraInterface.jar
    % file and other *.txt files with configuration to our MATLAB paths
    addpath('../Orchestra/');
    
    % Specify a path to an ORCHESTRA model files
    modelPath = '../Data/';

    % Read model parameters from a model files in the directory specified
    [chemistryFilePath, varDefinitionTable, ioVariableList] = ParseOrchestraModel(modelPath);
    
    % Then define an array of the names of IO variables. This is used in
    % order to set indices which will be used later on
    loopVars = {
        'H[Acetate].con'
        };
    
    % Initialize ORCHESTRA object:
    [orchestraInstance, loopVars] = InitializeOrchestraInterface(chemistryFilePath, varDefinitionTable, ioVariableList, loopVars);
    
    % Now you have an instance of the ORCHESTRA object initialized and may
    % call its method Calculate whenever needed.
    % This is one of the possibilities:
    
    % Create a vector of values for one of our loop vars
    hAcetateCon = 0.01:0.01:0.1;
    
    % Initialize results
    result = zeros(numel(hAcetateCon), numel(ioVariableList));
    
    for idx = 1:numel(hAcetateCon)
        % Run calculation method. it's inputs are:
        % - vector of indices of IO variables to be set (as in variable list);
        % - vector of values to be set for such variables;
        result(idx, :) = orchestraInstance.Calculate(loopVars(1).index, hAcetateCon(idx));
    end
    
    disp(ioVariableList');
    disp(result);
end