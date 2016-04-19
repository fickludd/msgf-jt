package uninovo.parser;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;

public class BufferedRandomAccessLineReader implements LineReader {
	private static final int DEFAULT_BUFFER_SIZE = 1 << 16;
	public static void main(String argv[]) throws Exception
	{
		long time = System.currentTimeMillis();
		String fileName = "/home/sangtaekim/Research/Data/ABRF/2011/UniProt.Yeast.NFISnr.contamsS48.fasta";
		BufferedRandomAccessLineReader in = new BufferedRandomAccessLineReader(fileName, 1 << 16);
//		BufferedReader in = new BufferedReader(new FileReader(fileName));
//		RandomAccessFile in = new RandomAccessFile(fileName, "r");
		String s;
		int lineNum=0;
		long pos = 0;
		while((s=in.readLine()) != null)
		{
			lineNum++;
			if(lineNum == 48232)
				System.out.println(lineNum+" "+s+" "+(pos=in.getPosition()));
		}
		in.seek(pos);
		System.out.println(in.readLine());
		System.out.println("Time: " + (System.currentTimeMillis()-time));
	}
	private byte[] buffer;
	private int bufLength = -1;
	int bufPointer;
	int bufSize;

	long bufStartingPos;
	private final byte CR = (byte)'\r';
	private long fileSize;
	private final FileChannel in;
	private final byte NL = (byte)'\n';
	private long pointer;
	
	int startIndex;
	
	public BufferedRandomAccessLineReader(String fileName)
	{
		this(fileName, DEFAULT_BUFFER_SIZE);
	}
	
	public BufferedRandomAccessLineReader(String fileName, int bufSize)
	{
		FileInputStream fin = null;
		try {
			fin = new FileInputStream(fileName);
		} catch (FileNotFoundException e1) {
			e1.printStackTrace();
		}

		in = fin.getChannel();
		try {
			fileSize = in.size();
		} catch (IOException e) {
			e.printStackTrace();
		}
	
		this.bufSize = bufSize;
		pointer = 0;
		fillBuffer();
	}
	
	public void close() throws IOException	{ in.close(); }
	
	private int fillBuffer()
	{
		ByteBuffer tempBuffer = null;
		int bytesRead = -1;
		try {
			tempBuffer = ByteBuffer.allocate(bufSize);
			bytesRead = in.read(tempBuffer);
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		
		buffer = tempBuffer.array();
		bufLength = bytesRead;
		startIndex = 0;
		bufPointer = 0;
		bufStartingPos = pointer;
		
		return bytesRead;
	}
	
	public long getPosition()	{ return pointer; }
	
	@Override
	public String readLine()
	{
		if(pointer >= fileSize)
			return null;
		String str = readLineFromBuffer();
		
		if(bufPointer == bufLength && bufLength == bufSize)
		{
			fillBuffer();
			str = str + readLine();
		}
		else if(pointer < fileSize)
		{
			bufPointer++;
			pointer++;
			startIndex = bufPointer;
		}
		return str;
	}
	
	private String readLineFromBuffer()	// line terminating char: \n or \r\n
	{
		if(pointer >= fileSize)
			return null;
		while(pointer < fileSize && bufPointer < bufLength)
		{
			if(buffer[bufPointer] != NL)
			{
				bufPointer++;
				pointer++;
			}
			else
				break;
		}
		
		String str;
		if(bufPointer > 0 && buffer[bufPointer-1] == CR)
			str = new String(buffer, startIndex, (bufPointer - startIndex - 1));
		else
			str = new String(buffer, startIndex, (bufPointer - startIndex));
		return str;
	}
	public void reset()		{ pointer = 0; startIndex = 0;}
	public void seek(long position)	
	{ 
		pointer = position;
		if(position >= bufStartingPos && position < bufStartingPos+bufSize)
		{
			startIndex = bufPointer = (int)(position - bufStartingPos);
		}
		else
		{
			try {
				in.position(pointer);
			} catch (IOException e) {
				e.printStackTrace();
			}
			fillBuffer();
		}
	}
	
	public int size()	{ return buffer.length; }
	
}
