function RunOrchestraCallFilesSample()

    % First we have to add the directory containing OrchestraInterface.jar
    % file and other *.txt files with configuration to our MATLAB paths
    addpath('../Orchestra/');
    
    % Specify a path to an ORCHESTRA model files
    modelPath = '../Data/';

    [chemistryFilePath, varDefinitionTable, ioVariableList] = ParseOrchestraModel(modelPath);
    
    % Then specify the name or index of IO variable we want to have a loop over
    loopVar = 'H[Acetate].liter';
    
    % And finally define the range of values for this variable
    loopVarVal = 0:0.02:0.1;
    
    % Call Orchestra interface
    result = RunOrchestraInterface(chemistryFilePath, varDefinitionTable, ioVariableList, loopVar, loopVarVal);
    
    disp(ioVariableList');
    disp(result);
end