import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.concurrent.atomic.AtomicInteger;

import org.junit.Test;

import gov.nih.ncats.molwitch.Chemical;

public class TestExoticAtoms {
	   @Test
       public void ensureWeightForStrangeAtoms() throws IOException {
           String mol= "[Am]";
                    
           
           Chemical c = Chemical.parse(mol);
           for(int i=1;i<=118;i++){
        	   c.getAtom(0).setAtomicNumber(i);
        	   double  mass = c.getMass();
        	   assertTrue(mass>0.01);
           }
       }
       
       @Test
       public void ensureWeightForStrangeAtomsTogether() throws IOException {
             
           
           Chemical c = Chemical.parse("\n" + 
           		"   JSDraw209302014022D\n" + 
           		"\n" + 
           		"  7  6  0  0  0  0            999 V2000\n" + 
           		"   17.5510   -7.6960    0.0000 Es  0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   23.0360   -5.6255    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   21.6850   -4.8455    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   21.6850   -3.2855    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   24.3870   -4.8455    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   24.3870   -3.2855    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   23.0360   -2.5055    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"  2  3  2  0  0  0  0\n" + 
           		"  3  4  1  0  0  0  0\n" + 
           		"  2  5  1  0  0  0  0\n" + 
           		"  5  6  2  0  0  0  0\n" + 
           		"  6  7  1  0  0  0  0\n" + 
           		"  7  4  2  0  0  0  0\n" + 
           		"M  END");

           assertEquals(330.112, c.getMass(),0.001);
           assertEquals("C6H6Es", c.getFormula());
           assertEquals(330.112, c.getMass(),0.001);
          
       }
       

       @Test
       public void ensureMolecularFormulaDoesntHaveCharges() throws IOException {
           String mol= "\n" + 
           		"   JSDraw209302010192D\n" + 
           		"\n" + 
           		" 16 16  0  0  0  0            999 V2000\n" + 
           		"   19.7870   -5.7460    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   18.4350   -4.9660    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   17.0840   -5.7460    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   18.4350   -3.4060    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   19.7870   -7.3060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   21.1370   -8.0860    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   22.4890   -7.3060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   23.8400   -8.0860    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   22.4890   -5.7460    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   21.1370   -9.6460    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   19.7870  -10.4260    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   18.4350   -9.6460    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   18.4350   -8.0860    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   19.7870  -11.9860    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   17.0840   -7.3060    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"   15.7330   -8.0860    0.0000 *   0  0  0  0  0  0  0  0  0  0  0  0\n" + 
           		"  1  2  1  0  0  0  0\n" + 
           		"  2  3  1  0  0  0  0\n" + 
           		"  2  4  2  0  0  0  0\n" + 
           		"  1  5  1  0  0  0  0\n" + 
           		"  5  6  2  0  0  0  0\n" + 
           		"  6  7  1  0  0  0  0\n" + 
           		"  7  8  1  0  0  0  0\n" + 
           		"  7  9  2  0  0  0  0\n" + 
           		"  6 10  1  0  0  0  0\n" + 
           		" 10 11  2  0  0  0  0\n" + 
           		" 12 11  1  0  0  0  0\n" + 
           		" 13 12  2  0  0  0  0\n" + 
           		"  5 13  1  0  0  0  0\n" + 
           		" 11 14  1  0  0  0  0\n" + 
           		"  3 15  2  0  0  0  0\n" + 
           		" 15 16  1  0  0  0  0\n" + 
           		"M  END";
                    
           Chemical c = Chemical.parse(mol);
           double  mass = c.getMass();
           
           assertEquals(270.056, mass,0.001);
           assertEquals("C10H6BrO4", c.getFormula());
           Iterator<Chemical> iter = c.copy().connectedComponents();
           final Map<String, AtomicInteger> formula = new HashMap<>();
           while(iter.hasNext()){
        	   Chemical c2 = iter.next();
               double  mass2 = c2.getMass();
               
               assertEquals(270.056, mass2,0.001);
               assertEquals("C10H6BrO4", c2.getFormula());
           }
       }
    
}
