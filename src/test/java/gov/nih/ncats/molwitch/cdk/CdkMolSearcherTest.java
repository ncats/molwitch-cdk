package gov.nih.ncats.molwitch.cdk;

import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.cdk.search.CdkMolSearcher;
import org.junit.Assert;
import org.junit.Test;

import java.util.Optional;
import java.util.logging.Logger;

public class CdkMolSearcherTest {

    private static Logger logger = Logger.getLogger("CdkMolSearcherTest");

    @Test
    public void smartsQueryReturnsExpectedMapping() throws Exception {
        CdkMolSearcher searcher = new CdkMolSearcher("[#7,#8]c1ccccc1");

        Optional<int[]> hit1 = searcher.search(Chemical.parse("Nc1ccccc1"));
        Assert.assertTrue(hit1.isPresent());
        Assert.assertArrayEquals(new int[]{0, 1, 6, 5, 4, 3, 2}, hit1.get());

        Optional<int[]> hit2 = searcher.search(Chemical.parse("CCCCCNc1ccccc1"));
        Assert.assertTrue(hit2.isPresent());
        Assert.assertArrayEquals(new int[]{5, 6, 11, 10, 9, 8, 7}, hit2.get());
    }

    @Test
    public void smartsQueryReturnsEmptyWhenNoMatch() throws Exception {
        CdkMolSearcher searcher = new CdkMolSearcher("[#7,#8]c1ccccc1");
        Assert.assertFalse(searcher.search(Chemical.parse("c1ccccc1")).isPresent());
    }
}
