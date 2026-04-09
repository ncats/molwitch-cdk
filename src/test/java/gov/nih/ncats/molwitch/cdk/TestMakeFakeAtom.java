package gov.nih.ncats.molwitch.cdk;

import gov.nih.ncats.molwitch.Atom;
import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.ChemicalBuilder;
import org.junit.Ignore;
import org.junit.Test;

import static org.junit.Assert.*;

public class TestMakeFakeAtom {
    @Test
    public void makeAtomWithASymbol(){
        //FDA often makes their pseudo atoms called A

        Chemical c = new Chemical();
        Atom atom = c.addAtom("A");
        assertEquals("A", atom.getSymbol());
    }

    @Test
    public void makeAtomWithASymbolBuilder(){
        //FDA often makes their pseudo atoms called A

        ChemicalBuilder builder = new ChemicalBuilder();

        builder.addAtom("A");
        Chemical c = builder.build();
        assertEquals("A", c.getAtom(0).getSymbol());
    }
}
