package edu.ucsd.msjava.msutil;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Iterator;

import uk.ac.ebi.jmzidml.model.mzidml.CvParam;

import edu.ucsd.msjava.mzid.Constants;
import edu.ucsd.msjava.parser.MS2SpectrumParser;
import edu.ucsd.msjava.parser.MgfSpectrumParser;
import edu.ucsd.msjava.parser.MzXMLSpectraIterator;
import edu.ucsd.msjava.parser.MzXMLSpectraMap;
import edu.ucsd.msjava.parser.PNNLSpectraIterator;
import edu.ucsd.msjava.parser.PNNLSpectraMap;
import edu.ucsd.msjava.parser.PklSpectrumParser;
import edu.ucsd.msjava.parser.SpectrumParser;

import se.lth.immun.msgf.MzMLCache;
import se.lth.immun.msgf.mzml.MzMLSpectraIterator;
import se.lth.immun.msgf.mzml.MzMLSpectraMap;

public class SpectraAccessor {
	private final File specFile;
	private final SpecFileFormat specFormat;
	private double intensityThreshold = 0;
	
	//private MzMLAdapter mzmlAdapter = null;
	private MzMLCache _mzmlCache = null;
	public MzMLCache mzmlCache(File f) {
		if (_mzmlCache == null)
			_mzmlCache = new MzMLCache(f, 2, 2, intensityThreshold);
		return _mzmlCache;
	}
	
	SpectrumAccessorBySpecIndex specMap = null;
	Iterator<Spectrum> specItr = null;
	
	public SpectraAccessor(File specFile) {
		this(specFile, SpecFileFormat.getSpecFileFormat(specFile.getName()),0.0);
	}
	
	public SpectraAccessor(File specFile, double intensityThreshold) {
		this(specFile, SpecFileFormat.getSpecFileFormat(specFile.getName()), intensityThreshold);
	}
	
	public SpectraAccessor(File specFile, SpecFileFormat specFormat, double intensityThreshold) {
		this.specFile = specFile;
		this.specFormat = specFormat;
		this.intensityThreshold = intensityThreshold;
	}
	
	public SpectrumAccessorBySpecIndex getSpecMap()
	{
		if(specMap == null)
		{
			if(specFormat == SpecFileFormat.MZXML)
				specMap = new MzXMLSpectraMap(specFile.getPath());
			else if(specFormat == SpecFileFormat.MZML || 
					specFormat == SpecFileFormat.MZML_GZ)
				specMap = new MzMLSpectraMap(mzmlCache(specFile));
			else if(specFormat == SpecFileFormat.DTA_TXT)
				specMap = new PNNLSpectraMap(specFile.getPath());
			else
			{
				SpectrumParser parser = null;
				if(specFormat == SpecFileFormat.MGF)
					parser = new MgfSpectrumParser();
				else if(specFormat == SpecFileFormat.MS2)
					parser = new MS2SpectrumParser();
				else if(specFormat == SpecFileFormat.PKL)
					parser = new PklSpectrumParser();
				else
					return null;
				specMap = new SpectraMap(specFile.getPath(), parser);
			}			
		}
		
		if(specMap == null)
		{
			System.out.println("File: " + specFile.getAbsolutePath());
			System.out.println("Format: " + specFormat.getPSIName());
		}
		return specMap;
	}
	
	public Iterator<Spectrum> getSpecItr()
	{
		if(specItr == null)
		{
			if(specFormat == SpecFileFormat.MZXML)
				specItr = new MzXMLSpectraIterator(specFile.getPath());
			else if(specFormat == SpecFileFormat.MZML || 
					specFormat == SpecFileFormat.MZML_GZ)
				specItr = new MzMLSpectraIterator(mzmlCache(specFile));
			else if(specFormat == SpecFileFormat.DTA_TXT)
				try {
					specItr = new PNNLSpectraIterator(specFile.getPath());
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				}
			else
			{
				SpectrumParser parser = null;
				if(specFormat == SpecFileFormat.MGF)
					parser = new MgfSpectrumParser();
				else if(specFormat == SpecFileFormat.MS2)
					parser = new MS2SpectrumParser();
				else if(specFormat == SpecFileFormat.PKL)
					parser = new PklSpectrumParser();
				else
					return null;
				try {
					specItr = new SpectraIterator(specFile.getPath(), parser);
				} catch (FileNotFoundException e) {
					e.printStackTrace();
				}
			}
		}
		
		return specItr;
	}
	
	public Spectrum getSpectrumBySpecIndex(int specIndex)
	{
		return getSpecMap().getSpectrumBySpecIndex(specIndex);
	}
	
	public Spectrum getSpectrumById(String specId)
	{
		return getSpecMap().getSpectrumById(specId);
	}
	
	public String getID(int specIndex)
	{
		return getSpecMap().getID(specIndex);
	}
	
	public float getPrecursorMz(int specIndex)
	{
		return getSpecMap().getPrecursorMz(specIndex);
	}
	
	public String getTitle(int specIndex)
	{
		return getSpecMap().getTitle(specIndex);
	}
	
	public CvParam getSpectrumIDFormatCvParam()
	{
		CvParam cvParam = null;
		if(specFormat == SpecFileFormat.DTA_TXT 
				|| specFormat == SpecFileFormat.MGF
				|| specFormat == SpecFileFormat.PKL
				|| specFormat == SpecFileFormat.MS2
				)
			cvParam = Constants.makeCvParam("MS:1000774", "multiple peak list nativeID format");
		else if(specFormat == SpecFileFormat.MZXML)
			cvParam = Constants.makeCvParam("MS:1000776", "scan number only nativeID format");
		else if(specFormat == SpecFileFormat.MZDATA)
			cvParam = Constants.makeCvParam("MS:1000777", "spectrum identifier nativeID format");
		else if(specFormat == SpecFileFormat.MZML || 
				specFormat == SpecFileFormat.MZML_GZ) {
			se.lth.immun.mzml.CvParam cvp = mzmlCache(specFile).spectrumIDFormatCvParam();
			cvParam = Constants.makeCvParam(cvp.accession(), cvp.name());
		}
			
		
		return cvParam;
	}
	
}
