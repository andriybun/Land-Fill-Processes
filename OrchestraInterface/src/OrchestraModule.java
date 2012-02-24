//import java.util.*;
import orchestra2.kernel.*;
import orchestra2.exception.*;

public class OrchestraModule {

	Calculator calculator;  // the complete chemical equilibrium module
	NodeType nodeType;      // defines the set of variables (names and default values) that are stored in a cell
	Node node;              // contains the actual values of the variables in this node
	int [] varHandles;      //indices to the different IO variables, for illustration purposes not in array 
	
	public OrchestraModule(String ChemistryFileName, NodeVariableInfoList variableList, String [] ioVariableList) {
		try {
			Initialize(ChemistryFileName, variableList, ioVariableList);
		} catch (ReadException e) {
			e.printStackTrace();
		}
	}
	
	private void Initialize(String ChemistryFileName, NodeVariableInfoList variableList, String [] ioVariableList) throws ReadException {
		
		// First we define a new chemical equilibrium calculator from its own input file
		// we could also 
		calculator = new Calculator(new FileID(null, ChemistryFileName));
		 
		// First we create a NodeType and defining the set of variables that is stored in each Node of this type
		// We give each variable a name, default value, indicate if it is a static variable (all nodes share same variable)
		// and indicate where this variable was defined. (Not really relevant here). 
		nodeType = new NodeType();
		
		for (int i = 0; i < variableList.GetList().length; i++) {
			NodeVariableInfo nd = variableList.GetList()[i]; 
			nodeType.addVariable(nd.varName, nd.defaultValue, nd.isStatic, nd.whereDefined);
		}
		
		// we ask the calculator which variables it wants to store in a cell between two calculations
		// to improve calculation efficiency
		nodeType.useGlobalVariablesFromCalculator(calculator);
		
		// After we have added all variables to the nodeType we can define an (or many) node of this type. 
		// this creates a node with default variables as given above.
		node = new Node(nodeType);
		
		// To write variables to a cell and read output from it we need to create an index to each of the IO variables.
		// We can store and reuse this index for fast access later  
		
		// Create the indices for the IO variables
		varHandles = new int[ioVariableList.length];
		for (int i = 0; i < ioVariableList.length; i++) {
			varHandles[i] = node.nodeType.index(ioVariableList[i]);
		}
	}
	
	public double [] Calculate(int [] varIndices, double [] varValues) throws ReadException {
		for (int i = 0; i < varIndices.length; i++) {
			node.setValue(varHandles[varIndices[i]-1], varValues[i]);
		}

		// Do a calculation
		try {
			calculator.calculate(node);
		} catch (ParserException e) {
			e.printStackTrace();
		} catch (ExitException e) {
			e.printStackTrace();
		}
		
		double [] result = new double[varHandles.length];
		for (int i = 0; i < varHandles.length; i++) {
			result[i] = node.getvalue(varHandles[i]);
		}
		return result;
			
	}

}
