
public class NodeVariableInfoList {

	private NodeVariableInfo [] list;
	
	public NodeVariableInfoList() {
		
	}
	
	public NodeVariableInfoList(int num) {
		SetList(num);
	}
	
	public void SetList(int num) {
		list = new NodeVariableInfo[num];
		for (int i = 0; i < num; i++) {
			list[i] = new NodeVariableInfo();
		}
	}

	public NodeVariableInfo [] GetList() {
		return list;
	}
	
	public void SetVariable(int idx, String varNameIn, double defaultValueIn, boolean isStaticIn, String whereDefinedIn) {
		list[idx-1].SetVariable(varNameIn, defaultValueIn, isStaticIn, whereDefinedIn);
	}
}
