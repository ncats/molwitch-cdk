/*
 * NCATS-MOLWITCH-CDK
 *
 * Copyright (c) 2024.
 *
 * This work is free software; you can redistribute it and/or modify it under the terms of the
 * GNU Lesser General Public License as published by the Free Software Foundation;
 * either version 2.1 of the License, or (at your option) any later version.
 *
 * This work is distributed in the hope that it will be useful, but without any warranty;
 * without even the implied warranty of merchantability or fitness for a particular purpose.
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with this library;
 *  if not, write to:
 *
 *  the Free Software Foundation, Inc.
 *  59 Temple Place, Suite 330
 *  Boston, MA 02111-1307 USA
 */

package gov.nih.ncats.molwitch.cdk;

import java.io.BufferedReader;
import java.io.ByteArrayInputStream;
import java.io.CharArrayReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;
import java.util.Scanner;

public class ProgramClearingMol2000Wrapper extends BufferedReader{

	private static String  M_END                = "M  END";
	private int currentLineInMol=0;
	private Integer lastMarkedLineInMol;
	public ProgramClearingMol2000Wrapper(Reader reader) throws IOException{
		super(reader);
	}
	
	

	@Override
	public void mark(int readAheadLimit) throws IOException {
		lastMarkedLineInMol = currentLineInMol;
		super.mark(readAheadLimit);
	}

	@Override
	public void reset() throws IOException {
		currentLineInMol = lastMarkedLineInMol;
		super.reset();
	}

	@Override
	public boolean markSupported() {
		return super.markSupported();
	}



	@Override
	public String readLine() throws IOException {
		
		String line= super.readLine();
		if(line !=null){
			currentLineInMol++;
		}
		if(currentLineInMol==2){
			//program line blank it out
			line = "";
		}
		if(M_END.equals(line)){
			currentLineInMol=-1;
		}
//		System.out.printf("line ='%s' count = %d%n", line, currentLineInMol);
		return line;
		
	}


	
	
	

}
