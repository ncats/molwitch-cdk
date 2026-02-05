package gov.nih.ncats.molwitch.cdk;

import org.junit.Assert;
import org.junit.Test;

public class CdkChemical2FactoryImplTest {

    @Test
    public void testIsSmarts() {
        String smiles = "[$([NX3](=O)=O),$([NX3+](=O)[O-])][!#8]";
        boolean result =  CdkChemical2FactoryImpl.looksLikeSmarts(smiles);
        Assert.assertTrue(result);
    }

    @Test
    public void testNotSmarts() {
        String smiles = "[C-]#[N+]C(F)(F)F";
        boolean result = CdkChemical2FactoryImpl.looksLikeSmarts(smiles);
        Assert.assertFalse(result);
    }
}
