package gov.nih.ncats.molwitch.cdk;

import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.cdk.search.CdkMolSearcher;
import org.apache.commons.io.IOUtils;
import org.junit.Assert;
import org.junit.Test;

import java.util.Optional;

public class CdkMolSearcherTest {

    private static final String ASSERTION_TITLE_PRESENT = "Search must return an array of int";
    private static final String ASSERTION_TITLE_MAPPING = "Search must return an array of atom mappings";

    @Test
    public void smartsQueryReturnsExpectedMapping() throws Exception {
        CdkMolSearcher searcher = new CdkMolSearcher("[#7,#8]c1ccccc1");

        Optional<int[]> hit1 = searcher.search(Chemical.parse("Nc1ccccc1"));
        Assert.assertTrue( hit1.isPresent());
        Assert.assertArrayEquals(ASSERTION_TITLE_PRESENT, new int[]{0, 1, 6, 5, 4, 3, 2}, hit1.get());

        Optional<int[]> hit2 = searcher.search(Chemical.parse("CCCCCNc1ccccc1"));
        Assert.assertTrue(ASSERTION_TITLE_PRESENT, hit2.isPresent());
        Assert.assertArrayEquals("Search must return an array of atom mappings", new int[]{5, 6, 11, 10, 9, 8, 7}, hit2.get());
    }

    @Test
    public void atomListQuery() throws  Exception {
        String molfileName = "5-membered-ring-atom-list.mol";
        String mol = IOUtils.toString(
                this.getClass().getResourceAsStream("/mols/" + molfileName),
                "UTF-8"
        );
        Chemical queryMol = Chemical.parseMol(mol);
        String target1 = IOUtils.toString(
                this.getClass().getResourceAsStream("/mols/pyrrole.mol" ),
                "UTF-8"
        );
        Chemical pyrrole = Chemical.parseMol(target1);

        String target2 = IOUtils.toString(
                this.getClass().getResourceAsStream("/mols/thf.mol" ),
                "UTF-8"
        );
        Chemical thf = Chemical.parseMol(target1);
        CdkMolSearcher searcher = new CdkMolSearcher(queryMol);
        Optional<int[]> hit1 = searcher.search(pyrrole);
        Assert.assertTrue("", hit1.isPresent());
        Assert.assertTrue(hit1.get().length > 0);
        Optional<int[]> hit2 = searcher.search(thf);
        Assert.assertTrue(hit2.get().length > 0);
    }

    @Test
    public void smartsQueryReturnsEmptyWhenNoMatch() throws Exception {
        CdkMolSearcher searcher = new CdkMolSearcher("[#7,#8]c1ccccc1");
        Assert.assertFalse(searcher.search(Chemical.parse("c1ccccc1")).isPresent());
    }

    @Test
    public void chemQueryFindsSimpleMatch() throws Exception {
        String mol = IOUtils.toString(
                this.getClass().getResourceAsStream("/mols/benzene.mol"),
                "UTF-8"
        );

        Chemical queryChemical = Chemical.parseMol(mol);
        CdkMolSearcher searcher = new CdkMolSearcher(queryChemical);

        mol = IOUtils.toString(
                this.getClass().getResourceAsStream("/mols/benzoic.acid.mol"),
                "UTF-8"
        );
        Chemical benzoicAcidHit = Chemical.parseMol(mol);
        Optional<int[]> hitMapping = searcher.search(benzoicAcidHit);
        Assert.assertEquals(6, hitMapping.get().length);
    }

    @Test
    public void chemQueryFindsNoMatch() throws Exception {
        String mol = IOUtils.toString(
                this.getClass().getResourceAsStream("/mols/coronene.mol"),
                "UTF-8"
        );

        Chemical queryChemical = Chemical.parseMol(mol);
        CdkMolSearcher searcher = new CdkMolSearcher(queryChemical);

        mol = IOUtils.toString(
                this.getClass().getResourceAsStream("/mols/benzoic.acid.mol"),
                "UTF-8"
        );
        Chemical benzoicAcidHit = Chemical.parseMol(mol);
        Optional<int[]> hitMapping = searcher.search(benzoicAcidHit);
        Assert.assertTrue(hitMapping.isEmpty());
    }

    @Test
    public void invalidSmartsThrowsOnParse(){
        String invalidSmarts = "-N(=O)=O";
        Exception result = Assert.assertThrows(RuntimeException.class, () -> {
            new CdkMolSearcher(invalidSmarts);
        } );
        Assert.assertNotNull(result.getMessage());
    }
}
