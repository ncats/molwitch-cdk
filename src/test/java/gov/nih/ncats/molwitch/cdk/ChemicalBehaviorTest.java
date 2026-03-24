package gov.nih.ncats.molwitch.cdk;

import gov.nih.ncats.molwitch.Atom;
import gov.nih.ncats.molwitch.Bond;
import gov.nih.ncats.molwitch.Chemical;
import org.junit.Test;

import java.io.IOException;
import java.io.InputStream;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

public class ChemicalBehaviorTest {

    @Test
    public void setAtomMapToPositionAssignsOneBasedMaps() throws Exception {
        Chemical chemical = Chemical.createFromSmiles("CCO");

        chemical.setAtomMapToPosition();

        assertTrue(chemical.hasAtomToAtomMappings());
        assertEquals(1, chemical.getAtom(0).getAtomToAtomMap().getAsInt());
        assertEquals(2, chemical.getAtom(1).getAtomToAtomMap().getAsInt());
        assertEquals(3, chemical.getAtom(2).getAtomToAtomMap().getAsInt());
    }

    @Test
    public void clearAtomMapsRemovesAllMaps() throws Exception {
        Chemical chemical = Chemical.createFromSmiles("CCO");
        chemical.setAtomMapToPosition();

        chemical.clearAtomMaps();

        assertFalse(chemical.hasAtomToAtomMappings());
        assertFalse(chemical.getAtom(0).getAtomToAtomMap().isPresent());
        assertFalse(chemical.getAtom(1).getAtomToAtomMap().isPresent());
        assertFalse(chemical.getAtom(2).getAtomToAtomMap().isPresent());
    }

    @Test
    public void computeShortestPathUsesMinimalNumberOfBonds() throws Exception {
        Chemical chemical = Chemical.createFromSmiles("CCCO");

        List<Bond> path = chemical.computeShortestPath(chemical.getAtom(0), chemical.getAtom(3));

        assertEquals(3, path.size());
    }

    @Test
    public void computeAllPathsReturnsBothRingRoutes() throws Exception {
        Chemical chemical = Chemical.createFromSmiles("C1CCC1");

        List<List<Bond>> paths = chemical.computeAllPaths(chemical.getAtom(0), chemical.getAtom(2));

        assertEquals(2, paths.size());
        assertEquals(2, paths.get(0).size());
        assertEquals(2, paths.get(1).size());
    }

    @Test
    public void computeStereochemistryTypeIdentifiesAchiralMolecule() throws Exception {
        Chemical achiral = Chemical.createFromSmiles("CCO");

        assertEquals(Chemical.StereochemistryType.ACHIRAL, achiral.computeStereochemistryType().get());
        assertEquals(Chemical.OpticalActivity.NONE, achiral.computeOpticalActivity().get());
    }

    @Test
    public void computeStereochemistryTypeIdentifiesAbsoluteFixture() throws Exception {
        Chemical chiral = loadMol("/mols/1-Phenylethanol-R.mol");

        assertEquals(Chemical.StereochemistryType.ABSOLUTE, chiral.computeStereochemistryType().get());
        assertEquals(Chemical.OpticalActivity.UNSPECIFIED, chiral.computeOpticalActivity().get());
    }

    @Test
    public void computeStereochemistryTypeIdentifiesRacemicFixture() throws Exception {
        Chemical racemic = loadMol("/mols/1-Phenylethanol-racemic.mol");

        assertEquals(Chemical.StereochemistryType.RACEMIC, racemic.computeStereochemistryType().get());
        assertEquals(Chemical.OpticalActivity.PLUS_MINUS, racemic.computeOpticalActivity().get());
    }

    private Chemical loadMol(String resourcePath) throws IOException {
        try (InputStream in = getClass().getResourceAsStream(resourcePath)) {
            return Chemical.parseMol(in);
        }
    }
}
