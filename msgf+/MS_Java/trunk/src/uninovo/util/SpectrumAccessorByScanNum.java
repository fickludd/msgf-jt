package uninovo.util;

import java.util.ArrayList;

public interface SpectrumAccessorByScanNum {
	public ArrayList<Integer> getScanNumList();
	public Spectrum getSpectrumByScanNum(int scanNum);
}
