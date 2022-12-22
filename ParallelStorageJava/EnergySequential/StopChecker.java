package EnergySequential;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

/**
 * <p>
 * Title: MPJTest2
 * </p>
 * <p>
 * Description: UPS Load Planning
 * </p>
 * <p>
 * Copyright: Copyright (c) 2012
 * </p>
 * <p>
 * Company: Princeton University
 * </p>
 * 
 * @author Belgacem Bouzaiene-Ayari
 * @version 5.40
 */
public class StopChecker extends Thread {
	private volatile File file;
	private volatile boolean done;

	/**
	 * 
	 */
	public StopChecker(File file_) {
		done = false;
		file = file_;
	}

	/**
	 * 
	 */
	public void end() {
		done = true;
	}

	/**
	 * 
	 */
	public void run() {
		while (!done) {
			System.out.println(" sleeping for 10 secs. . .");
			System.out.flush();
			try {
				Thread.sleep(20000);
			} catch (InterruptedException e) {
				e.printStackTrace();
			}
			BufferedReader reader = null;
			try {
				reader = new BufferedReader(new FileReader(file));
			} catch (FileNotFoundException e1) {
				e1.printStackTrace();
			}
			if (reader != null) {
				System.out.println(" checking stop file. . .");
				System.out.flush();
				try {
					String line = reader.readLine();
					System.out.println("  line: "+line);
					System.out.flush();
					if (line != null && line.trim().startsWith("1")) {
						// stopping
						System.out.println("  will end!");
						System.out.flush();
						break;
					}
				} catch (IOException e) {
					e.printStackTrace();
				}
				try {
					reader.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
			}
		}

		if (!done) {
			// exiting!
			System.exit(-1);
		}
	}
}
