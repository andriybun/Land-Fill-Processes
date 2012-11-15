function ORI = initialize_MSWS2(Comp,x)

% Specify the path to a compiled JAR file with the ORCHESTRA interface.
% And import corresponding libraries
javaclasspath([cd '/Orchestra_MSWS2/OrchestraInterface.jar']);
import OrchestraInterface.*;

% We initialize a list of Master species and I/O variables
% We give each variable a name, default value, indicate if it is a static variable
% and indicate where this variable was defined.

variableList = NodeVariableInfoList(length(Comp.all));

for i = 1:length(Comp.all)
    variableList.SetVariable(i, Comp.all(i), Comp.alli(i), false, 'defined by example');
    ioVariableList(i) = Comp.all(i);
end

% Initialize ORCHESTRA object
% To initialize the total amount of H2CO3 at initial pCO2
if x ~= 0
    k1 = find(strcmp('H2CO3.logact',Comp.all));
    ORI = OrchestraModule([cd '/Orchestra_MSWS2/chemistry1.inp'], variableList, ioVariableList);
    ORI = ORI.Calculate([k1], [Comp.alli(k1)]); %(previous method)
% To initialize ORCHESTRA including total amount of H2CO3
else 
    ORI = OrchestraModule([cd '/Orchestra_MSWS2/chemistry.inp'], variableList, ioVariableList);
end

end