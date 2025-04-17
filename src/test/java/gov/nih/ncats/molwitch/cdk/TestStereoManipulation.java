package gov.nih.ncats.molwitch.cdk;

import gov.nih.ncats.molwitch.Bond;
import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.MolwitchException;
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

    private class MolEquivalenceTestData {
        public MolEquivalenceTestData(Chemical chemical1, Chemical chemical2, boolean expectedToBeEquivalent) {
            this.chemical1 = chemical1;
            this.chemical2 = chemical2;
            this.expectedToBeEquivalent = expectedToBeEquivalent;
        }

        public Chemical getChemical1() {
            return chemical1;
        }

        public void setChemical1(Chemical chemical1) {
            this.chemical1 = chemical1;
        }

        public Chemical getChemical2() {
            return chemical2;
        }

        public void setChemical2(Chemical chemical2) {
            this.chemical2 = chemical2;
        }

        public boolean isExpectedToBeEquivalent() {
            return expectedToBeEquivalent;
        }

        public void setExpectedToBeEquivalent(boolean expectedToBeEquivalent) {
            this.expectedToBeEquivalent = expectedToBeEquivalent;
        }

        private Chemical chemical1;
        private Chemical chemical2;
        private boolean expectedToBeEquivalent;

    }

    @Test
    public void testSet() throws IOException {
        List<TestMol> testMols = generateTestMolSet();
        testMols.forEach(m->{
            try {
                String molfileText = IOUtils.toString(
                    this.getClass().getResourceAsStream("/mols/" + m.molfileName + ".mol"),
                    "UTF-8"
                );
                Chemical before= Chemical.parse(molfileText);
                Chemical flipped = CdkChemicalImpl.flipStereocenters(before);
                Assert.assertEquals(m.expectedDown, flipped.bonds().filter(b->b.getStereo() == Bond.Stereo.DOWN ).count());
                Assert.assertEquals(m.expectedUp, flipped.bonds().filter(b->b.getStereo() == Bond.Stereo.UP).count());
            } catch (IOException e) {
                System.err.println("Error processing molfile "+ e.getMessage());
                Assert.fail("Error processing molfile " + m.molfileName);
            }
        });
    }

    @Test
    public void testpermuteEpimers() throws IOException {
        List<TestMol> testMols = generateEpimericTestMolSet();
        testMols.forEach(m-> {
            System.out.printf("Starting test for %s\n", m.getMolfileName());
            try {
                String molfileText = IOUtils.toString(
                        this.getClass().getResourceAsStream("/mols/" + m.molfileName + ".mol"),
                        "UTF-8"
                );
                Chemical before = Chemical.parse(molfileText);
                CdkChemicalImpl b= (CdkChemicalImpl) before.getImpl();
                List<Chemical> epimers= b.permuteEpimers();
                if( epimers.size()== 0){
                    System.out.printf("no epimers found for %s\n", m.getMolfileName());
                }
                Assert.assertTrue(epimers.size()>0);
            } catch (IOException e) {
                System.err.println("Error processing molfile " + e.getMessage());
                Assert.fail("Error processing molfile " + m.molfileName);
            }
        });
    }

    @Test
    public void testSet2() throws IOException {
        List<TestMol> testMols = generateTestMolSet();
        testMols.forEach(m->{
            try {
                System.out.printf("testing mol %s\n", m.molfileName);
                String molfileText = IOUtils.toString(
                        this.getClass().getResourceAsStream("/mols/" + m.molfileName + ".mol"),
                        "UTF-8"
                );
                Chemical chemical = Chemical.parse(molfileText);
                Chemical flippedChemical= chemical.getImpl().flipAllChiralCenters();
                if(m.expectedDown !=  getStereoBondCount(flippedChemical, Bond.Stereo.DOWN)
                    || m.expectedUp != getStereoBondCount(flippedChemical, Bond.Stereo.UP)){
                    File fileOutputFile = new File(m.molfileName + "_mod.mol");
                    Files.writeString(fileOutputFile.toPath(), chemical.toMol());
                }
                Assert.assertEquals(m.expectedDown, getStereoBondCount(flippedChemical, Bond.Stereo.DOWN)); //chemical.bonds().filter(b->b.getStereo() == Bond.Stereo.DOWN ).count()
                Assert.assertEquals(m.expectedUp, getStereoBondCount( flippedChemical, Bond.Stereo.UP));
            } catch (IOException e) {
                System.err.println("Error processing molfile "+ e.getMessage());
                Assert.fail("Error processing molfile " + m.molfileName);
            }
        });
    }

    @Test
    public void testEquivalentTo() throws IOException {
        generateEquivalenceTestData().forEach(m->{
            CdkChemicalImpl impl = (CdkChemicalImpl)m.chemical1.getImpl();
            try {
                impl.generateCoordinates();
                m.chemical2.getImpl().generateCoordinates();
                boolean actual = impl.equivalentTo(m.chemical2);
                if( actual != m.expectedToBeEquivalent) {
                    System.out.printf("equivalence error on %s and %s \n", m.chemical1.toMol(), m.chemical2.toMol());
                }
                Assert.assertEquals(m.expectedToBeEquivalent, actual);
            } catch (IOException | MolwitchException e) {
                System.err.println("Error processing structure "+ e.getMessage());
                Assert.fail();
            }
        });
    }
    private int getStereoBondCount(Chemical chemical, Bond.Stereo requiredStereo){
        int count = 0;
        for(int iBond = 0;iBond<chemical.getBondCount(); iBond++) {
            if(chemical.getBond(iBond).getStereo() == requiredStereo) count++;
        }
        return count;
    }

    private List<TestMol> generateTestMolSet(){
        return Arrays.asList(
                //new TestMol("1-Phenylethanol-R", 1, 0),
                //new TestMol("1-Phenylethanol-S", 0, 1),
                new TestMol("1-Phenylethanol-racemic", 0, 0),
                //new TestMol("33W7SJ9TBX", 0, 2),
                new TestMol("VG7S7JRA56", 6, 3),
                new TestMol("DRR5D9W4K6", 2, 4),
                new TestMol("N6WK7SF4JA", 5, 3),
                new TestMol("G783UGT4GL", 1, 1),
                new TestMol("L6N99T5XF6", 0, 0)
        );
    }

    private List<TestMol> generateEpimericTestMolSet(){
        return Arrays.asList(
                new TestMol("1-Phenylethanol-racemic", 0, 0),
                new TestMol("R949HH57L8", 2, 4),
                new TestMol("ZU7D9H38JX", 0, 2),
                new TestMol("M6T5F9662S", 3,6)
        );
    }
    private List<MolEquivalenceTestData> generateEquivalenceTestData() throws IOException {
        return Arrays.asList(
                new MolEquivalenceTestData(Chemical.parseMol("CCCCCNC"),
                        Chemical.parseMol("CCCCCNC"), true),
                new MolEquivalenceTestData(Chemical.parseMol("CCCCCNC"), Chemical.parseMol("CCCCCNCC"), false),
                new MolEquivalenceTestData(Chemical.parseMol("CC(O)CCCNC"), Chemical.parseMol("CC(O)CCCNC"), true),
                new MolEquivalenceTestData(Chemical.parseMol("C[C@H](O)CCCNC"), Chemical.parseMol("C[C@@H](O)CCCNC"), false),
                new MolEquivalenceTestData(Chemical.parseMol("C[C@@H](O)CCCNC"), Chemical.parseMol("C[C@H](O)CCCNC"), false),
                new MolEquivalenceTestData(Chemical.parseMol("C[C@@H](O)CCOCNC"), Chemical.parseMol("C[C@@H](O)CCOCNC"), true)
        );
    }
}
