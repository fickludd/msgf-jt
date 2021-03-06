package edu.ucsd.msjava.msutil;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import edu.ucsd.msjava.params.ParamObject;
import edu.ucsd.msjava.params.UserParam;


public class Protocol implements ParamObject {
	private String name;
	private String description;
	
	private Protocol(String name, String description)
	{
		this.name = name;
		this.description = description;
	}
	
	public String getName()		{ return name; }
	public String getDescription()	{ return description; }
	
	public String getParamDescription() {
		return name;
	}
	
	// static members
	public static Protocol get(String name)
	{
		return table.get(name);
	}
	
	public static final Protocol NOPROTOCOL;
	public static final Protocol PHOSPHORYLATION;
	public static final Protocol ITRAQ;
	public static final Protocol ITRAQPHOSPHO;
	public static final Protocol TMT;
	
	public static Protocol[] getAllRegisteredProtocols()
	{
		return protocolList.toArray(new Protocol[0]);
	}
	
	private static HashMap<String, Protocol> table;
	private static ArrayList<Protocol> protocolList;
	private static void add(Protocol prot)
	{
		if(table.put(prot.name, prot) == null)
			protocolList.add(prot);
	}
	
	static {
		NOPROTOCOL = new Protocol("NoProtocol", "No protocol");
		PHOSPHORYLATION = new Protocol("Phosphorylation", "Phospho-enriched");
		ITRAQ = new Protocol("iTRAQ", "iTRAQ");
		ITRAQPHOSPHO = new Protocol("iTRAQPhospho", "iTRAQPhospho");
		TMT = new Protocol("TMT", "TMT");
		
		table = new HashMap<String, Protocol>();
		protocolList = new ArrayList<Protocol>();
		
		protocolList.add(NOPROTOCOL);
		add(PHOSPHORYLATION);
		add(ITRAQ);
		add(ITRAQPHOSPHO);
		add(TMT);
		
		// Parse activation methods defined by a user
		File protocolFile = new File("params/protocols.txt");
		if(protocolFile.exists())
		{
			ArrayList<String> paramStrs = UserParam.parseFromFile(protocolFile.getPath(), 2);
			for(String paramStr : paramStrs)
			{
				String[] token = paramStr.split(",");
				String shortName = token[0];
				String description = token[1];
				Protocol newProt = new Protocol(shortName, description);
				add(newProt);
			}
		}
	}
	
}
