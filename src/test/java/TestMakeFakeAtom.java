/*
 * NCATS-MOLWITCH-CDK
 *
 * Copyright (c) 2025.
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

import gov.nih.ncats.molwitch.Atom;
import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.ChemicalBuilder;
import org.junit.Ignore;
import org.junit.Test;

import static org.junit.Assert.*;
@Ignore
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
