import gov.nih.ncats.molwitch.Bond;
import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.cdk.CdkChemicalImpl;
import org.apache.commons.io.IOUtils;
import org.junit.Assert;
import org.junit.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.List;

public class TestStereoManipulation {

    private class TestMol{
        private String molfileName;

        public TestMol(String molfileName, int expectedUp, int expectedDown) {
            this.molfileName = molfileName;
            this.expectedUp = expectedUp;
            this.expectedDown = expectedDown;
        }

        public String getMolfileName() {
            return molfileName;
        }

        public void setMolfileName(String molfileName) {
            this.molfileName = molfileName;
        }

        public int getExpectedUp() {
            return expectedUp;
        }

        public void setExpectedUp(int expectedUp) {
            this.expectedUp = expectedUp;
        }

        public int getExpectedDown() {
            return expectedDown;
        }

        public void setExpectedDown(int expectedDown) {
            this.expectedDown = expectedDown;
        }

        private int expectedUp;
        private int expectedDown;
    }

    @Test
    public void testSet() throws IOException {
        List<TestMol> testMols = Arrays.asList(
                new TestMol("1-Phenylethanol-R", 1, 0),
                new TestMol("1-Phenylethanol-S", 0, 1),
                new TestMol("1-Phenylethanol-racemic", 0, 0),
                new TestMol("33W7SJ9TBX", 0, 2),
                new TestMol("VG7S7JRA56", 6, 3),
                new TestMol("DRR5D9W4K6", 2, 4),
                new TestMol("N6WK7SF4JA", 5, 3),
                new TestMol("G783UGT4GL", 1, 1),
                new TestMol("L6N99T5XF6", 0, 0)
        );
        testMols.forEach(m->{
            try {
                String molfileText = IOUtils.toString(
                    this.getClass().getResourceAsStream("mols/" + m.molfileName + ".mol"),
                    "UTF-8"
                );
                Chemical before= null;
                before = Chemical.parse(molfileText);
                Chemical flipped = CdkChemicalImpl.flipStereocenters(before);
                Assert.assertEquals(m.expectedDown, flipped.bonds().filter(b->b.getStereo() == Bond.Stereo.DOWN ).count());
                Assert.assertEquals(m.expectedUp, flipped.bonds().filter(b->b.getStereo() == Bond.Stereo.UP).count());
            } catch (IOException e) {
                System.err.println("Error processing molfile "+ e.getMessage());
                Assert.fail("Error processing molfile " + m.molfileName);
            }
        });

    }
}
