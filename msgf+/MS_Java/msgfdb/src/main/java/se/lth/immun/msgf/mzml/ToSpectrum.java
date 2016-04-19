package se.lth.immun.msgf.mzml;

import java.util.Collections;

import se.lth.immun.msgf.MzMLCache.SpecData;
import se.lth.immun.msgf.MzMLCache.SpecPrec;
import edu.ucsd.msjava.msutil.ActivationMethod;
import edu.ucsd.msjava.msutil.Peak;
import edu.ucsd.msjava.msutil.Spectrum;

public class ToSpectrum {
	
	static Spectrum from(int index, SpecData sd) {
		Spectrum s = new Spectrum();
		
		s.setID(sd.id());
		s.setScanNum(sd.scanNum());
		s.setMsLevel(sd.msLevel());
		s.setIsCentroided(sd.isCentroided());
		
		SpecPrec sp = sd.precs()[0];
		s.setIsolationWindowTargetMz((float)sp.isolationWindowTargetMz());
		s.setPrecursor(new Peak((float)sp.precMz(), (float)sp.precInt(), sp.precZ()));
		for (String acc : sp.activationAcc()) {
			ActivationMethod am = ActivationMethod.getByCV(acc);
			if (am != null) {
				s.setActivationMethod(am);
				if (am == ActivationMethod.ETD)
					break;
			}
		}
		
		for (int i=0; i<sd.mzs().length; i++)
			s.add(new Peak((float)sd.mzs()[i], (float)sd.intensities()[i], 1));
		
		s.setSpecIndex(index);
		
		Collections.sort(s);
		s.determineIsCentroided();
		
		return s;
	}
}
