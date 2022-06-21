import java.util.ArrayList;
import java.util.Arrays;




public class convertCoordinates {

	private static String convertLine(String line) {
		String[] fields = line.split("\t");
		String[] chrom = fields[2].split("[-]");
		if(line.charAt(0) == '@'){
			return line;
		}
		else if ( fields[5].equals("*") && ( fields[2].length() != 1 || Integer.parseInt(fields[3]) != 0 ) ){
				fields[2] = "*";
				fields[3] = "0";
				String newline = Arrays.toString(fields).replace(", ", "\t").replaceAll("[\\[\\]]", "");
				return newline;
				
			}
		else if (chrom.length == 1) {
			return line;
		}
		else if(chrom.length == 100){
			
			int start1 = Integer.parseInt(chrom[1]);
			int stop1 = Integer.parseInt(chrom[2]);
			int start2 = Integer.parseInt(chrom[3]);
			int beginning = Integer.parseInt(fields[3]);
			int read1start = start1 + beginning - 1;
            int read1length = stop1 - read1start + 1;
            int Nlength = start2-stop1-1;
            String cigar = fields[5];
            String newcigar = "";
            
            
            
            if(cigar.equals("*")){
            	newcigar = "*";
            }
            else if (read1length <= 0) {
				newcigar = cigar;
				read1start = start2 + beginning - stop1 + start1 - 2;
			}
            else{
            	String[] cigarnumsTmp = cigar.split("[MIDNSHP]");
                int[] cigarnums = new int[cigarnumsTmp.length];
                for(int i = 0; i < cigarnumsTmp.length; i++){
                	cigarnums[i] = Integer.parseInt(cigarnumsTmp[i]);
                }
                String cigarlettersTmp[] = cigar.split("[0-9]+");
                String cigarletters[] = Arrays.copyOfRange(cigarlettersTmp, 1, cigarlettersTmp.length);
            	int currentpos = 0;
            	boolean found = false;
            	
            	for(int i = 0; i < cigarnums.length; i++){
            		if(!cigarletters[i].equals("I") && !cigarletters[i].equals("S") && !cigarletters[i].equals("H")){
            			currentpos += cigarnums[i];
            		}
            		
            		if(currentpos > read1length && !found){
            			found = true;
            			int len1 = cigarnums[i] - (currentpos - read1length);
            			int len2 = cigarnums[i] - len1;
            			newcigar = newcigar + String.valueOf(len1) + cigarletters[i] + String.valueOf(Nlength) + "N" + String.valueOf(len2 + cigarletters[i]);
            		}
            		else{
            			newcigar = newcigar + String.valueOf(cigarnums[i]) + cigarletters[i];
            		}
            	}
            }                
            String sep = "\t";
            String newLine = fields[0] + sep + fields[1] + sep + chrom[0] + sep + read1start + sep + fields[4] + sep + newcigar + sep;
            for(int i = 6; i < fields.length-1; i++){
            	newLine += fields[i]+sep;
            }
            newLine += fields[fields.length-1];
            
            return newLine;

		}
		else{
			
			String newline = parseMultiexon(fields,chrom);
			return newline;
			
			
		}
	}

	
	
	private static String parseMultiexon(String[] fields, String[] chrom) {
		ArrayList<Integer> exonStarts = new ArrayList<Integer>();
		ArrayList<Integer> exonEnds = new ArrayList<Integer>();
		ArrayList<Integer> cmLength = new ArrayList<Integer>();
		
		int readStartRel = Integer.parseInt(fields[3]);
		String newcigar = "";
		
		/*
		 * get numbers and characters from cigar string
		 */
		String cigar = fields[5];
		String[] cigarnumsTmp = cigar.split("[MIDNSHP]");
        int[] cigarnums = new int[cigarnumsTmp.length];
        int cigarsum = 0;
        for(int i = 0; i < cigarnumsTmp.length; i++){
        	cigarnums[i] = Integer.parseInt(cigarnumsTmp[i]);
        	cigarsum += cigarnums[i]; 
        }
        String cigarlettersTmp[] = cigar.split("[0-9]+");
        String cigarletters[] = Arrays.copyOfRange(cigarlettersTmp, 1, cigarlettersTmp.length);
		
        /*
		 * get exon starts and ends as well as individual and cumulative exon lengths
		 */
		int len = 0;
		for(int i = 1; i < chrom.length; i++){
			if( i % 2 == 0){
				exonEnds.add(Integer.parseInt(chrom[i]));
				len += exonEnds.get(exonEnds.size()-1) - exonStarts.get(exonStarts.size()-1) +1;
				cmLength.add(new Integer(len));
			}
			else{
				exonStarts.add(Integer.parseInt(chrom[i]));
			}
		}
		/*
		 * find the exon that contains the read start
		 */
		int readStartIdx = -1;
		int readStartAbs = -1;
		if(readStartRel <= cmLength.get(0)){
			readStartIdx = 0;
			readStartAbs = exonStarts.get(0) + readStartRel - 1;
		}
		for(int i = 1; i < cmLength.size(); i++){
			if(readStartRel <= cmLength.get(i) && readStartRel > cmLength.get(i-1)){
				readStartIdx = i;
				readStartAbs = exonStarts.get(i) + (readStartRel-cmLength.get(i-1))-1;
			}
		}

        
        if(cigar.equals("*")){
        	newcigar = "*";
        }
        else{
        	int currCigarPos = 0;
        	int currReadPos = exonEnds.get(readStartIdx) - readStartAbs + 1;
            int readRemain = cigarsum;
        	boolean found = false;

        	for(int i = 0; i < cigarnums.length; i++){
        		
        		if(!cigarletters[i].equals("I") && !cigarletters[i].equals("S") && !cigarletters[i].equals("H")){
        			currCigarPos += cigarnums[i];

        		}
        		
        		while(currCigarPos > currReadPos && !found && readStartIdx < exonEnds.size()-1){

 
        			int len1 = cigarnums[i] - (currCigarPos - currReadPos);
        			cigarnums[i] -= len1;
        			int Nlength = exonStarts.get(readStartIdx+1) - exonEnds.get(readStartIdx) - 1 ;
        			newcigar = newcigar + String.valueOf(len1) + cigarletters[i] + String.valueOf(Nlength) + "N";
        			readRemain -= len1;
        			readStartIdx++;
        			if(exonEnds.get(readStartIdx) - exonStarts.get(readStartIdx) + 1 < readRemain){
        				currReadPos += (exonEnds.get(readStartIdx) - exonStarts.get(readStartIdx) + 1);
        			}
        			else{
        				currReadPos += readRemain;
        			}
        			
        		}

        		newcigar = newcigar + String.valueOf(cigarnums[i]) + cigarletters[i];
        	}
        }

		
		String sep = "\t";
		String newLine = fields[0] + sep + fields[1] + sep + chrom[0] + sep + readStartAbs + sep + fields[4] + sep + newcigar + sep;
        for(int i = 6; i < fields.length-1; i++){
        	newLine += fields[i]+sep;
        }
        newLine += fields[fields.length-1];
        
		return newLine;
	}



	/**
	 * @param args
	 */
	public static void main(String[] args) {
		try {
	       java.io.BufferedReader stdin = new java.io.BufferedReader(new java.io.InputStreamReader(System.in));
	       String line;
	       
	       
			while((line = stdin.readLine()) !=null){
				 String convertedLine = convertLine(line);
				 
				 System.out.println(convertedLine);
			}
		}	  
		catch (java.io.IOException e) { System.out.println(e); }
		catch (NumberFormatException e) { System.out.println(e); }


	}


}
