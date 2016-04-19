package suffixtree.actions;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.TreeMap;

import msgap.Parameters;
import msgap.ScoringParameter;
import msgap.ScoringParameterIterator;
import msgap.results.GappedPeptideResults;
import msgap.results.SpectrumMatches;
import msgf.AminoAcidGraph;
import msgf.GeneratingFunction;
import msgf.NominalMass;
import msgf.NominalMassFactory;
import msutil.AminoAcid;
import msutil.AminoAcidSet;
import msutil.Annotation;
import msutil.Composition;
import msutil.Peptide;
import sequences.MassSequence;
import suffixtree.Constants;
import suffixtree.matches.ExactMatchObject;
import suffixtree.matches.MatchObject;
import suffixtree.matches.ModMatchObject;
import suffixtree.matches.MutMatchObject;
import suffixtree.misc.ProgressMeter;
import suffixtree.trees.FRKeywordTreeCompact;
import suffixtree.trees.KeywordTreeCompact;
import uninovoOld.DenovoReconstruction;
import uninovoOld.ScoredSpectrumForMSGF;
import uninovoOld.ScoringData;



public class Scoring {
  
  private static int readExactMatches(
      GappedPeptideResults gpr, 
      String matchFile, 
      TreeMap<Integer,ArrayList<MatchObject>> matches
      ) {
    
    String matchObjectText;
    int resultsCount = 0;
    try {
      BufferedReader br = new BufferedReader(new FileReader(matchFile));
      while ((matchObjectText = br.readLine()) != null) { // while loop begins here
        ExactMatchObject mo = new ExactMatchObject(gpr.getSequences(), matchObjectText);
        resultsCount++;
        int specId = gpr.getSpecId(mo.getQueryIndex());
        if (!matches.containsKey(specId)) matches.put(specId, new ArrayList<MatchObject>());
        matches.get(specId).add(mo);
      } // end while 
    } // end try
    catch (IOException e) {
      System.err.println("Error: " + e);
    }
    return resultsCount;
  }
  
  
  /**
   * Scores and prints out the modified matches. This implementation prints out
   * the best match(es) with a valid probability
   * @param params the parameter object
   * @param sIt the score parameter iterator
   * @param gpr the gap result object
   * @param out the file to write the results
   * @param totalMatches total number of matches stored in the match directory
   * @param matches the specId to matches map
   */
  public static void score(Parameters params,
                           Iterator<ScoringParameter> sIt,
                           GappedPeptideResults gpr,
                           String matchDir,
                           int totalMatches,
                           PrintWriter out) {

    ProgressMeter pm = new ProgressMeter("\nScoring spectrum-matches", totalMatches, System.out);
    
    int matchFileIndex = 0;
    TreeMap<Integer,ArrayList<MatchObject>> matches = new TreeMap<Integer,ArrayList<MatchObject>>();
    File matchFile = new File(matchDir, String.format(KeywordTreeCompact.MATCH_FILE_PREFIX+"%04d.txt", matchFileIndex++));
    readExactMatches(gpr, matchFile.getAbsolutePath(), matches);
    
    //int specCount = 0;
    long time = System.currentTimeMillis();
    int scoredMatches = 0;

    while (sIt.hasNext()) {
      ScoringParameter sp = sIt.next();
      int specId = sp.getSpecID();
      
      if(matches.isEmpty() || specId > matches.lastKey()) {
        matches.clear();
        matchFile = new File(matchDir, String.format(KeywordTreeCompact.MATCH_FILE_PREFIX+"%04d.txt", matchFileIndex++));
        if (!matchFile.exists()) break;
        //System.out.println(matchFile.getAbsolutePath());
        readExactMatches(gpr, matchFile.getAbsolutePath(), matches);
      }
        
      //boolean hasMatch = false;
      if (matches.containsKey(specId)) {
        // score all the matches
        SpectrumMatches sm = new SpectrumMatches(sp, params);
        
        float bestProb = sm.addMatches(matches.get(specId));
        scoredMatches += matches.get(specId).size();
        pm.update(scoredMatches);
        
        if (bestProb < 1.0f) {
          
          for (MatchObject mo : sm.getMatches()) {
            
            if (mo.getProb()<=bestProb) {
              // offset is the ParentMass - Theoretical Mass
              float offset = sp.getOriginalParentMass()-mo.getPeptide().getMass()-(float)Composition.H2O;
              String line = mo.getSummaryLine(gpr.getFileName(specId), gpr.getScanNumber(specId), gpr.getActmethod(specId), gpr.getPrecursorMass(specId), gpr.getCharge(specId), offset);
              out.println(line);
            }
          }
        }
       
      }
      if (specId >= gpr.getMaxId()) {
        if (specId > gpr.getMaxId()) {
          System.out.println("WARNING: Potential spectrum that was not scored");
        }
        break; // we are done iterating the iterator
      }
    }
    System.out.println();
    
    time = System.currentTimeMillis() - time;
    System.out.printf("Scored %d matches in %d seconds.\n", scoredMatches, time/1000);
  }
  
  /**
   * Scores and prints out the modified matches. This implementation prints out
   * the best match(es) with a valid probability
   * @param gpr the gap result object
   * @param out the file to write the results
   * @param totalMatches total number of matches stored in the match directory
   * @param matches the specId to matches map
   * added by Kyowon for adaNovo
   */
  public static void score(GappedPeptideResults gpr,
          String matchDir,
          Iterator<ScoringData> sd,
          int totalMatches,
          AminoAcidSet aaSet,
          PrintWriter out) {

    ProgressMeter pm = new ProgressMeter("\nScoring spectrum-matches", totalMatches, System.out);
    NominalMassFactory factory  = null;
    
    int matchFileIndex = 0;
    TreeMap<Integer,ArrayList<MatchObject>> matches = new TreeMap<Integer,ArrayList<MatchObject>>();
    File matchFile = new File(matchDir, String.format(KeywordTreeCompact.MATCH_FILE_PREFIX+"%04d.txt", matchFileIndex++));
    readExactMatches(gpr, matchFile.getAbsolutePath(), matches);
    
    //int specCount = 0;
    long time = System.currentTimeMillis();
    int scoredMatches = 0;
    //matches.firstKey()
    while (sd.hasNext()) {
    	ScoringData s = sd.next();
        int specId = s.getSpecID();
        
        if(matches.isEmpty() || specId > matches.lastKey()) {
          matches.clear();
          matchFile = new File(matchDir, String.format(FRKeywordTreeCompact.MATCH_FILE_PREFIX+"%04d.txt", matchFileIndex++));
          if(!matchFile.exists()) break;
          
          //System.out.println(matchFile.getAbsolutePath());
          readExactMatches(gpr, matchFile.getAbsolutePath(), matches);
        }
        
        ScoredSpectrumForMSGF<NominalMass> sg = new ScoredSpectrumForMSGF<NominalMass>(s.getGraph());
		if(factory == null) factory = new NominalMassFactory(aaSet, s.getGraph().getEnzyme(), 70);
     	 AminoAcidGraph ag =  new AminoAcidGraph(factory, s.getGraph().getParentMass(), sg);
			
		 GeneratingFunction<NominalMass> gf = new GeneratingFunction<NominalMass>(ag).enzyme(s.getGraph().getEnzyme());
		 gf.doNotBacktrack().doNotCalcNumber();
		 gf.computeGeneratingFunction();
          
        //boolean hasMatch = false;
        if (matches.containsKey(specId)) {
        // score all the matches
         
	        ArrayList<MatchObject> mos = matches.get(specId);
	        scoredMatches += matches.get(specId).size();
	        pm.update(scoredMatches);
	        
	        float minProb = 1;
	        float massError = 0;
	        
	        MatchObject bestMo = null;
	        
		      for (MatchObject mo : mos) {  
		    	  if(!mo.isCleaved(s.getGraph().getEnzyme())) continue;
		    	  Peptide pep = new Peptide(mo.getPeptide().toString(), aaSet);
		    	  
		    	  for(Peptide mpep : pep.getModifiedPeptides(aaSet)){
		    		
			     	 if(!s.getGraph().getSinkNode().isCorrect(mpep)) continue;
			     	  
					 Annotation annotation = new Annotation(aaSet.getAminoAcid(mo.getLeftFlankingAA()), mpep,  aaSet.getAminoAcid(mo.getRightFlankingAA()));
					
					 int score = gf.getScore(annotation);// s.getGraph().getScore(pep);
					 float prob = gf.getSpectralProbability(Math.min(score, gf.getMaxScore()-1));
				//	 System.out.println( ag.getScore(pep));
					// System.out.println(annotation.getPeptide().getMass()+"\t" + gf.getMaxScore() + "\t" + s.getGraph().getSinkNode().getMass());
						
					 if(minProb > prob){
						 boolean  isCorrect = false;
				     	 for(DenovoReconstruction g : s.getRecs()){
				     		 if(g.isCorrect(pep)){
				     			 isCorrect = true;
				     			 
				     			 break;
				     		 }
				     	 }
				     	 
				     	 if(!isCorrect) continue;
				     	 
				     	
						 mo.setProb(prob);
						 mo.setScore(score);
						 mo.setPepideWithExactModifications(mpep);
						 // System.out.println(mpep+"\t"+annotation + "\t" + score);
					    	
						 bestMo = mo;
						 minProb = prob;
						 massError = s.getGraph().getSinkNode().getMass() - mpep.getMass();
					 }
		    	  }
				 /*
		     	 
		     	 if(maxScore < score){
		     		 mo.setProb(score);
		     		 bestMo = mo;
		     		 maxScore = score;
		     		massError = s.getGraph().getSinkNode().getMass() - pep.getMass();
		     	 }
		     	 */
		          
		      }
		      if(bestMo != null){
		    	  String line = bestMo.getSummaryLine(gpr.getFileName(specId), gpr.getScanNumber(specId), gpr.getActmethod(specId), gpr.getPrecursorMass(specId), gpr.getCharge(specId), massError);
			    	 out.println(line);
	          }
       
      }
      if (specId >= gpr.getMaxId()) {
        if (specId > gpr.getMaxId()) {
          System.out.println("WARNING: Potential spectrum that was not scored");
        }
        break; // we are done iterating the iterator
      }
    }
    System.out.println();
    
    time = System.currentTimeMillis() - time;
    System.out.printf("Scored %d matches in %d seconds.\n", scoredMatches, time/1000);
  }
  

  private static int readModMatches(MassSequence db, 
                                    GappedPeptideResults gpr, 
                                    String matchFile, 
                                    TreeMap<Integer,ArrayList<MatchObject>> matches) {
    String matchObjectText;
    int resultsCount = 0;
    try {
      BufferedReader br = new BufferedReader(new FileReader(matchFile));
      while ((matchObjectText = br.readLine()) != null) { // while loop begins here
        ModMatchObject mo = new ModMatchObject(db, gpr.getSequences(), matchObjectText);
        resultsCount++;
        int specId = gpr.getSpecId(mo.getQueryIndex());
        if (!matches.containsKey(specId)) matches.put(specId, new ArrayList<MatchObject>());
        matches.get(specId).add(mo);
      } // end while 
    } // end try
    catch (IOException e) {
      System.err.println("Error: " + e);
    }
    return resultsCount;
  }
  
  
  private static int readMutMatches (
      MassSequence db, 
      GappedPeptideResults gpr, 
      String matchFile, 
      TreeMap<Integer,ArrayList<MatchObject>> matches
      ) {

    String matchObjectText;
    int resultsCount = 0;
    try {
      BufferedReader br = new BufferedReader(new FileReader(matchFile));
      while ((matchObjectText = br.readLine()) != null) { // while loop begins here
        MutMatchObject mo = new MutMatchObject(db, gpr.getSequences(), matchObjectText);
        resultsCount++;
        int specId = gpr.getSpecId(mo.getQueryIndex());
        if (!matches.containsKey(specId)) matches.put(specId, new ArrayList<MatchObject>());
        matches.get(specId).add(mo);
      } // end while 
    } // end try
    catch (IOException e) {
      System.err.println("Error: " + e);
    }
    return resultsCount;
  }
  
  
  
  /**
   * Scores and prints out the mutated matches. This implementation prints out
   * the best match(es) with a valid probability
   * @param params the parameter object
   * @param sIt the score parameter iterator
   * @param gpr the gap result object
   * @param out the file to write the results
   * @param totalMatches total number of matches stored in the match directory
   * @param matches the specId to matches map
   */
  public static void scoreMutatedMatches(Parameters params,
                                         MassSequence db,
                                         Iterator<ScoringParameter> sIt,
                                         GappedPeptideResults gpr,
                                         String matchDir,
                                         int totalMatches,
                                         PrintWriter out) {

    ProgressMeter pm = new ProgressMeter("\nScoring spectrum-matches", totalMatches, System.out);
    
    int matchFileIndex = 0;
    TreeMap<Integer,ArrayList<MatchObject>> matches = new TreeMap<Integer,ArrayList<MatchObject>>();
    File matchFile = new File(matchDir, String.format(FRKeywordTreeCompact.MATCH_FILE_PREFIX+"%04d.txt", matchFileIndex++));
    readMutMatches(db, gpr, matchFile.getAbsolutePath(), matches);
    
    //int specCount = 0;
    long time = System.currentTimeMillis();
    int scoredMatches = 0;

    while (sIt.hasNext()) {
      ScoringParameter sp = sIt.next();
      int specId = sp.getSpecID();
      
      if(matches.isEmpty() || specId > matches.lastKey()) {
        matches.clear();
        matchFile = new File(matchDir, String.format(FRKeywordTreeCompact.MATCH_FILE_PREFIX+"%04d.txt", matchFileIndex++));
        //System.out.println(matchFile.getAbsolutePath());
        readMutMatches(db, gpr, matchFile.getAbsolutePath(), matches);
      }
        
      //boolean hasMatch = false;
      if (matches.containsKey(specId)) {
        // score all the matches
        SpectrumMatches sm = new SpectrumMatches(sp, params);
        
        float bestProb = sm.addMatches(matches.get(specId));
        scoredMatches += matches.get(specId).size();
        pm.update(scoredMatches);
        
        if (bestProb < 1.0f) {
          
          for (MatchObject mo : sm.getMatches()) {
            
            if (mo.getProb()<=bestProb) {
              // offset is the ParentMass - Theoretical Mass
              float offset = sp.getOriginalParentMass()-mo.getPeptide().getMass()-(float)Composition.H2O;
              String line = mo.getSummaryLine(gpr.getFileName(specId), gpr.getScanNumber(specId), gpr.getActmethod(specId), gpr.getPrecursorMass(specId), gpr.getCharge(specId), offset);
              out.println(line);
            }
          }
        }
       
      }
      if (specId >= gpr.getMaxId()) {
        if (specId > gpr.getMaxId()) {
          System.out.println("WARNING: Potential spectrum that was not scored");
        }
        break; // we are done iterating the iterator
      }
    }
    System.out.println();
    
    time = System.currentTimeMillis() - time;
    System.out.printf("Scored %d matches in %d seconds.\n", scoredMatches, time/1000);
  }
  
  /**
   * Scores and prints out the mutated matches. This implementation prints out
   * the best match(es) with a valid probability
   * @param gpr the gap result object
   * @param out the file to write the results
   * @param totalMatches total number of matches stored in the match directory
   * @param matches the specId to matches map
   * added by Kyowon
   */
  public static void scoreMutatedMatches(MassSequence db,
          GappedPeptideResults gpr,
          String matchDir,
          Iterator<ScoringData> sd,
          int totalMatches,
          AminoAcidSet aaSet,
          PrintWriter out) {

    ProgressMeter pm = new ProgressMeter("\nScoring spectrum-matches", totalMatches, System.out);
    //NominalMassFactory factory  = null;
    
    int matchFileIndex = 0;
    TreeMap<Integer,ArrayList<MatchObject>> matches = new TreeMap<Integer,ArrayList<MatchObject>>();
    File matchFile = new File(matchDir, String.format(FRKeywordTreeCompact.MATCH_FILE_PREFIX+"%04d.txt", matchFileIndex++));
    readMutMatches(db, gpr, matchFile.getAbsolutePath(), matches);
    
    //int specCount = 0;
    long time = System.currentTimeMillis();
    int scoredMatches = 0;

    while (sd.hasNext()) {
    	ScoringData s = sd.next();
        int specId = s.getSpecID();
        
        if(matches.isEmpty() || specId > matches.lastKey()) {
          matches.clear();
          matchFile = new File(matchDir, String.format(FRKeywordTreeCompact.MATCH_FILE_PREFIX+"%04d.txt", matchFileIndex++));
          if(!matchFile.exists()) break;
          
         // System.out.println(matchFile.getAbsolutePath());
          readMutMatches(db, gpr, matchFile.getAbsolutePath(), matches);
        }
          
      //  ScoredSpectrumForMSGF<NominalMass> sg = new ScoredSpectrumForMSGF<NominalMass>(s.getGraph());
	//	if(factory == null) factory = new NominalMassFactory(aaSet, s.getGraph().getEnzyme(), 70);
     //	 AminoAcidGraph ag =  new AminoAcidGraph(factory, s.getGraph().getParentMass(), sg);
			
	//	 GeneratingFunction<NominalMass> gf = new GeneratingFunction<NominalMass>(ag).enzyme(s.getGraph().getEnzyme());
	//	 gf.doNotBacktrack().doNotCalcNumber();
	//	 gf.computeGeneratingFunction();
          
       
        //boolean hasMatch = false;
        if (matches.containsKey(specId)) {
        
	    	ArrayList<MatchObject> mos = matches.get(specId);
	        scoredMatches += mos.size();
	        pm.update(scoredMatches);
	          
	        //System.out.println(s.getGraph().getAllNodes());
	      float massError = 0;
	 //     float minProb = 1;
	 //     float delta = 1000f;
	      
	      String line = null;
	      HashSet<String> matchedPepStrings = new HashSet<String>();
	      for (MatchObject mo : mos) {
	    	  if(!mo.isCleaved(s.getGraph().getEnzyme())) continue;
	     	 Peptide pep = new Peptide(mo.getPeptide().toString(), aaSet);
	     	 
	     	 for(Peptide mpep : pep.getModifiedPeptides(aaSet)){
		    	 if(!s.getGraph().getSinkNode().isCorrect(mpep)) continue;
		     	 
		 //    	 Annotation annotation = new Annotation(aaSet.getAminoAcid(mo.getLeftFlankingAA()), mpep,  aaSet.getAminoAcid(mo.getRightFlankingAA()));
			//	 int score =  gf.getScore(annotation);//s.getGraph().getScore(pep);
			//	 float prob = gf.getSpectralProbability(Math.min(score, gf.getMaxScore()-1));
				 massError = s.getGraph().getSinkNode().getMass() - mpep.getMass();
				 
			//	 if(minProb >= prob){
				//	 if(minProb == prob && massError > delta) continue;
					 
				/*	 boolean  isCorrect = false;
			    	 for(DenovoReconstruction g : s.getRecs()){
			    		 if(g.isCorrect(mpep)){
			    			 isCorrect = true;
			    			 
			    			 break;
			    		 }
			    	 }
			    	 
			    	 if(!isCorrect) continue;
			    	*/ 
				 
		
				 if(matchedPepStrings.contains(mpep.toString())) continue;
		    		
				//	 mo.setProb(prob);
				//	 mo.setScore(score);
					 mo.setPepideWithExactModifications(mpep);
			//		 minProb = prob;
				//	 delta = massError;
					 line = mo.getSummaryLine(gpr.getFileName(specId), gpr.getScanNumber(specId), gpr.getActmethod(specId), gpr.getPrecursorMass(specId), gpr.getCharge(specId), massError);
						if(line != null){
			    			out.println(line);
			    			matchedPepStrings.add(mpep.toString());
			    		}
		     
				// }
	     	 }
        }
	   // if(line != null) out.println(line);
	    	  
       
      }
      if (specId >= gpr.getMaxId()) {
        if (specId > gpr.getMaxId()) {
          System.out.println("WARNING: Potential spectrum that was not scored");
        }
        break; // we are done iterating the iterator
      }
    }
    System.out.println();
    
    time = System.currentTimeMillis() - time;
    System.out.printf("Scored %d matches in %d seconds.\n", scoredMatches, time/1000);
  }
  
  /**
   * Scores and prints out the modified matches. This implementation prints out
   * the best match(es) with a valid probability
   * @param params the parameter object
   * @param sIt the score parameter iterator
   * @param gpr the gap result object
   * @param out the file to write the results
   * @param totalMatches total number of matches stored in the match directory
   * @param matches the specId to matches map
   */
  public static void scoreModdedMatches(Parameters params,
                                        MassSequence db,
                                        Iterator<ScoringParameter> sIt,
                                        GappedPeptideResults gpr,
                                        String matchDir,
                                        int totalMatches,
                                        PrintWriter out) {

    ProgressMeter pm = new ProgressMeter("\nScoring spectrum-matches", totalMatches, System.out);
    
    int matchFileIndex = 0;
    TreeMap<Integer,ArrayList<MatchObject>> matches = new TreeMap<Integer,ArrayList<MatchObject>>();
    File matchFile = new File(matchDir, String.format(FRKeywordTreeCompact.MATCH_FILE_PREFIX+"%04d.txt", matchFileIndex++));
    readModMatches(db, gpr, matchFile.getAbsolutePath(), matches);
    
    //int specCount = 0;
    long time = System.currentTimeMillis();
    int scoredMatches = 0;

    while (sIt.hasNext()) {
      ScoringParameter sp = sIt.next();
      int specId = sp.getSpecID();
      
      if(matches.isEmpty() || specId > matches.lastKey()) {
        matches.clear();
        matchFile = new File(matchDir, String.format(FRKeywordTreeCompact.MATCH_FILE_PREFIX+"%04d.txt", matchFileIndex++));
        //System.out.println(matchFile.getAbsolutePath());
        readModMatches(db, gpr, matchFile.getAbsolutePath(), matches);
      }
        
      //boolean hasMatch = false;
      if (matches.containsKey(specId)) {
        // score all the matches
        SpectrumMatches sm = new SpectrumMatches(sp, params);
        
        float bestProb = sm.addMatches(matches.get(specId));
        scoredMatches += matches.get(specId).size();
        pm.update(scoredMatches);
        
        if (bestProb < 1.0f) {
          
          for (MatchObject mo : sm.getMatches()) {
            
            if (mo.getProb()<=bestProb) {
              // offset is the ParentMass - Theoretical Mass
              float offset = sp.getOriginalParentMass()-mo.getPeptide().getMass()-(float)Composition.H2O;
              String line = mo.getSummaryLine(gpr.getFileName(specId), gpr.getScanNumber(specId), gpr.getActmethod(specId), gpr.getPrecursorMass(specId), gpr.getCharge(specId), offset);
              out.println(line);
            }
          }
        }
       
      }
      if (specId >= gpr.getMaxId()) {
        if (specId > gpr.getMaxId()) {
          System.out.println("WARNING: Potential spectrum that was not scored");
        }
        break; // we are done iterating the iterator
      }
    }
    System.out.println();
    
    time = System.currentTimeMillis() - time;
    System.out.printf("Scored %d matches in %d seconds.\n", scoredMatches, time/1000);
  }
  
  /**
   * Scores and prints out the modified matches. This implementation prints out
   * the best match(es) with a valid probability
   * @param gpr the gap result object
   * @param out the file to write the results
   * @param totalMatches total number of matches stored in the match directory
   * @param matches the specId to matches map
   * added by Kyowon
   */
  public static void scoreModdedMatches(MassSequence db,
                                        GappedPeptideResults gpr,
                                        String matchDir,
                                        Iterator<ScoringData> sd,
                                        int totalMatches,
                                        AminoAcidSet aaSet,
                                        PrintWriter out) {

    ProgressMeter pm = new ProgressMeter("\nScoring spectrum-matches", totalMatches, System.out);
   // NominalMassFactory factory  = null;
    
    int matchFileIndex = 0;
    TreeMap<Integer,ArrayList<MatchObject>> matches = new TreeMap<Integer,ArrayList<MatchObject>>();
    File matchFile = new File(matchDir, String.format(FRKeywordTreeCompact.MATCH_FILE_PREFIX+"%04d.txt", matchFileIndex++));
    readModMatches(db, gpr, matchFile.getAbsolutePath(), matches);
    
    //int specCount = 0;
    long time = System.currentTimeMillis();
    int scoredMatches = 0;
 
    while (sd.hasNext()) {
      ScoringData s = sd.next();
      int specId = s.getSpecID();
      
      if(matches.isEmpty() || specId > matches.lastKey()) {
    	 // System.out.println(matchFileIndex + "\t" + specId + "\t" + matches.lastKey());
          
        matches.clear();
        matchFile = new File(matchDir, String.format(FRKeywordTreeCompact.MATCH_FILE_PREFIX+"%04d.txt", matchFileIndex++));
        if(!matchFile.exists()) break;
        readModMatches(db, gpr, matchFile.getAbsolutePath(), matches);
      }
        
      //ScoredSpectrumForMSGF<NominalMass> sg = new ScoredSpectrumForMSGF<NominalMass>(s.getGraph());
	//	if(factory == null) factory = new NominalMassFactory(aaSet, s.getGraph().getEnzyme(), 70);
	//	AminoAcidGraph ag =  new AminoAcidGraph(factory, s.getGraph().getParentMass(), sg);
			
	//	 GeneratingFunction<NominalMass> gf = new GeneratingFunction<NominalMass>(ag).enzyme(s.getGraph().getEnzyme());
	//	 gf.doNotBacktrack().doNotCalcNumber();
	//	 gf.computeGeneratingFunction();
       
		 
      //boolean hasMatch = false;
      if (matches.containsKey(specId)) {
        // score all the matches
    	ArrayList<MatchObject> mos = matches.get(specId);
        scoredMatches += mos.size();
        pm.update(scoredMatches);
        
       //System.out.println(s.getGraph().getAllNodes());
   //   float minProb = 1;
     
      String line = null;
   //   float delta = 1000f;
      
      HashSet<String> matchedPepStrings = new HashSet<String>();
      for (MatchObject mo : mos) {
    	  if(!mo.isCleaved(s.getGraph().getEnzyme())) continue;
    	 
    	  Peptide pep = new Peptide(mo.getUnmodifiedPeptide().toString(), aaSet);
     	 
    	  ModMatchObject mmo =  (ModMatchObject)mo;
    	  
    	  
    	 for(Peptide mpep : pep.getModifiedPeptides(aaSet)){
    		 float offset = s.getGraph().getSinkNode().getMass() - mpep.getMass();
    		
    		 if(offset < Constants.MIN_MOD || offset > Constants.MAX_MOD) continue;
    		
    		 for( int modIndex = mmo.getModStart(); modIndex < mmo.getModEnd();modIndex++){
    		     Peptide tmpep = new Peptide(mpep);
    		    
    		     if(!s.getGraph().getSinkNode().isCorrect(mpep)){ // if blind modified
		    		// if(mo.getPeptide().isModified()){
	    			 AminoAcid aa = tmpep.get(modIndex);
	    			// System.out.println(mpep);
	    			 if(!aa.isModified()){ 
		    			 if(aa.getMass() + offset < Constants.MIN_AMINO_ACID_MASS) continue;
		    			 tmpep.remove(modIndex);
		    			 tmpep.add(modIndex, AminoAcid.getCustomAminoAcid(Character.toLowerCase(aa.getResidue()), "", aa.getMass() + offset));
	    			 }
		    		//}
    		     }else continue; // only blind modified peptides survive
    			
					 
	    		 if(!s.getGraph().getSinkNode().isCorrect(tmpep)) continue; 
	    		 	    		 
	    	//	 Annotation annotation = new Annotation(aaSet.getAminoAcid(mo.getLeftFlankingAA()), tmpep,  aaSet.getAminoAcid(mo.getRightFlankingAA()));
 			    
		  //  	 int score = gf.getScore(annotation);// s.getGraph().getScore(pep);
		//		 float prob = gf.getSpectralProbability(Math.min(score, gf.getMaxScore()-1));
				
		//		 if(minProb >= prob){
		//			 if(minProb == prob && Math.abs(offset) > delta) continue;
					 
				/*	 boolean  isCorrect = false;
			    	 for(DenovoReconstruction g : s.getRecs()){
			    		 if(g.isCorrect(tmpep)){
			    			 isCorrect = true;
			    			 
			    			 break;
			    		 }
			    	 }
			    	 
			    	 if(!isCorrect) continue;*/
			    	 
			    	 //mo.setProb(prob);
					// mo.setScore(score);
	    		 if(matchedPepStrings.contains(tmpep.toString())) continue;
	    		 
					 mo.setPepideWithExactModifications(tmpep);
					 
			    	 line = mo.getSummaryLine(gpr.getFileName(specId), gpr.getScanNumber(specId), gpr.getActmethod(specId), gpr.getPrecursorMass(specId), gpr.getCharge(specId), offset);

			    		if(line != null){
			    			out.println(line);
			    			matchedPepStrings.add(tmpep.toString());
			    		}
					
				/*	 minProb = prob;
					 delta = Math.abs(offset);
					 if(modIndex != mmo.getModEnd()+1){
						 mo.setProb(prob);
						 mo.setScore(score);
						 mo.setPepideWithExactModifications(tmpep);
						 line = mo.getSummaryLine(gpr.getFileName(specId), gpr.getScanNumber(specId), gpr.getActmethod(specId), gpr.getPrecursorMass(specId), gpr.getCharge(specId), offset);
					 }
				 }*/
    		 }
    	 }
        }
 	 
      }
      if (specId >= gpr.getMaxId()) {
        if (specId > gpr.getMaxId()) {
          System.out.println("WARNING: Potential spectrum that was not scored");
        }
        break; // we are done iterating the iterator
      }
    }
    System.out.println();
    
    time = System.currentTimeMillis() - time;
    System.out.printf("Scored %d matches in %d seconds.\n", scoredMatches, time/1000);
  }

  /**
   * Scores and prints out the matches. This implementation prints out
   * the best match(es) with a valid probability
   * @param params the parameter object
   * @param sIt the score parameter iterator
   * @param gpr the gap result object
   * @param out the file to write the results
   * @param matches the specId to matches map
   */
  public static void score(Parameters params,
                           Iterator<ScoringParameter> sIt,
                           GappedPeptideResults gpr,
                           PrintWriter out, 
                           HashMap<Integer,ArrayList<MatchObject>> matches) {

    ProgressMeter pm = new ProgressMeter("\nScoring spectrum-matches", matches.size(), System.out);
    
    int specCount = 0;
    long time = System.currentTimeMillis();
    int scoredMatches = 0;

    while (!matches.isEmpty() && sIt.hasNext()) {
      ScoringParameter sp = sIt.next();
      int specId = sp.getSpecID();
      
      //boolean hasMatch = false;
      if (matches.containsKey(specId)) {
        // score all the matches
        SpectrumMatches sm = new SpectrumMatches(sp, params);
        pm.update(++specCount);
        
        float bestProb = sm.addMatches(matches.get(specId));
        scoredMatches += matches.get(specId).size();
        
        if (bestProb < 1.0f) {
          
          for (MatchObject mo : sm.getMatches()) {
            
            if (mo.getProb()<=bestProb) {
              // offset is the ParentMass - Theoretical Mass
              float offset = sp.getOriginalParentMass()-mo.getPeptide().getMass()-(float)Composition.H2O;
              String line = mo.getSummaryLine(gpr.getFileName(specId), gpr.getScanNumber(specId), gpr.getActmethod(specId), gpr.getPrecursorMass(specId), gpr.getCharge(specId), offset);
              out.println(line);
            }
          }
        }
       
      }
      if (specId >= gpr.getMaxId()) {
        if (specId > gpr.getMaxId()) {
          System.out.println("WARNING: Potential spectrum that was not scored");
        }
        break; // we are done iterating the iterator
      }
    }
    System.out.println();
    
    time = System.currentTimeMillis() - time;
    System.out.printf("Scored %d matches in %d seconds.\n", scoredMatches, time/1000);
  }
  
  
  
  // score a particular spectrum
  public static void score(Parameters params, String filename, int scanNum) {
    Iterator<ScoringParameter> i = new ScoringParameterIterator(params);
    
    String peptideF = "R.GLGILDSALNELQGDTLDGETVFK.L";
   // String peptide = "GLGILDSALNELQGDTLDGETVFK";
    while (i.hasNext()) {
      
      ScoringParameter sp = i.next();
      String[] tokens = sp.getSpecFileName().split("/");
      String name = tokens[tokens.length-1].trim();
      //System.out.println(sp.getScanNum());
      //System.out.println(name);
      if (name.equals(filename) && scanNum==sp.getScanNum()) {
        System.out.println(name); 
        SpectrumMatches sm = new SpectrumMatches(sp, params);
        Annotation pep = new Annotation(peptideF, params.aaSet());
        int score = sm.getScore(pep);
        float prob = sm.getProbability(pep, score);
        System.out.printf("%s\t%d\t%.3e\n", name, score, prob);
      }
      
    }
  }
  
  
  public static void main(String[] args) {
    String userHome = System.getProperty("user.home");
    
    String[] in3 = {String.format("OutputFile %s/Data/Spectra/Sone/LTQFT3/output6", userHome), 
                    String.format("DBFile %s/Data/Databases/Sone/pro/SOne_uniprot_plus_contaminants.fasta", userHome),
                    String.format("InputFile %s/Data/Spectra/Sone/LTQFT3", userHome)};
    score(new Parameters(in3), "ShewFed037_LTQFT_1_12Nov04_Pegasus_0804-4_dta.ms2", 10998);  
  }
  
  
  
}
