
public class NodeVariableInfo {

	String varName;
	double defaultValue;
	boolean isStatic;
	String whereDefined;
	
	public NodeVariableInfo(String varNameIn, double defaultValueIn, boolean isStaticIn, String whereDefinedIn) {
		SetVariable(varNameIn, defaultValueIn, isStaticIn, whereDefinedIn);
	}
	
	public NodeVariableInfo() {
		varName = "";
		defaultValue = 0.0;
		isStatic = false;
		whereDefined = "";
	}
	
	public void SetVariable(String varNameIn, double defaultValueIn, boolean isStaticIn, String whereDefinedIn) {
		varName = varNameIn;
		defaultValue = defaultValueIn;
		isStatic = isStaticIn;
		whereDefined = whereDefinedIn;
	}
	
}
