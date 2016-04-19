package uninovo.util;

import java.util.Hashtable;


/**
 * A class representing a modification.
 * @author sangtaekim
 *
 */
public class Modification {
	/**
	 * A class representing the modification instance.
	 * @author sangtaekim
	 *
	 */
	public static class Instance {
		private boolean isFixedModification = false;
		private Location location;	// N_Term, C_Term, Anywhere
		private final Modification mod;
		private final char residue;	// if null, no amino acid specificity
		public Instance(Modification mod, char residue)
		{
			this(mod, residue, Location.Anywhere);
		}
		public Instance(Modification mod, char residue, Location location)
		{
			this.mod = mod;
			this.residue = residue;
			this.location = location;
		}
		public Instance fixedModification()	{ isFixedModification = true; return this; }
		
		public Location getLocation()	{ return location; }
		public Modification getModification()	{ return mod; }
		public char getResidue()	{ return residue; }
		public boolean isFixedModification() { return isFixedModification; }
	}
	public static enum Location {
		Anywhere,
		C_Term,
		N_Term,
		Protein_C_Term,
		Protein_N_Term
	}
	// static member
	private static Modification[] modList = 
	{
		new Modification("Carbamidomethylation", new Composition(2,3,1,1,0)),
		new Modification("Carboxymethylation", new Composition(2,2,2,0,0)),
		new Modification("Oxidation", new Composition(0,0,0,1,0)),
		new Modification("Phosphorylation", new Composition(0,1,0,3,0)),
	};
	
	private static Hashtable<String,Modification> modTable;
	
	static {
		modTable = new Hashtable<String, Modification>();
		for(Modification mod : modList)
			modTable.put(mod.getName(), mod);
	}
	
	public static Modification get(String name) { return modTable.get(name); }
	public static Modification register(String name, Composition composition)
	{
		Modification mod = new Modification(name, composition);
		modTable.put(name, mod);
		return mod;
	}
	public static Modification register(String name, double mass)
	{
		Modification mod = new Modification(name, mass);
		modTable.put(name, mod);
		return mod;
	}
	private Composition composition;
	
	private double mass;

	private final String name;
	
	private Modification(String name, Composition composition)
	{
		this.name = name;
		this.composition = composition;
	}

	private Modification(String name, double mass)
	{
		this.name = name;
		this.composition = null;
		this.mass = mass;
	}
	
	public double getAccurateMass() { return mass; } 
	public Composition	getComposition()	{ return composition; }
	
	public float getMass() { return (float)mass; }
	
	public String getName()	{ return name; }
}
