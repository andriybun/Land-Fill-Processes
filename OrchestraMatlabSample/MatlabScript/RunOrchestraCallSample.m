function RunOrchestraCallSample()

    % First we have to add the directory containing OrchestraInterface.jar
    % file and other *.txt files with configuration to our MATLAB paths
    addpath('../Orchestra/');
    
    % Then specify the path to chemistry.inp
    chemistryFilePath = '/../Data/chemistry.inp';

    % We define the table of variables specifying a variable name and its
    % default value
    varDefinitionTable = {
        'Ca+2.liter', 1
        'H+.liter', -2.5
    	'CO2[g].tot', 11
        'I', 0.814
        };
    
    % Specify an array with the names of IO variables
    ioVariableList = {
        'H[Acetate].liter'
        'pH'
        };
    
    % Initialize Orchestra interface
    [orchestraInstance, loopVars] = InitializeOrchestraInterface(chemistryFilePath, varDefinitionTable, ioVariableList, {'H[Acetate].liter'});
    
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