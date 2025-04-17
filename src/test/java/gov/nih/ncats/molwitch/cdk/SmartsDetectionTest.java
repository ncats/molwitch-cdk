package gov.nih.ncats.molwitch.cdk;

import org.junit.Assert;
import org.junit.Test;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class SmartsDetectionTest {

    private static class smartsTestData {
        public static smartsTestData of(String input, boolean expected) {
            return new smartsTestData(input, expected);
        }

        public smartsTestData(String testInput, boolean expectedToBeSmarts) {
            this.testInput = testInput;
            this.expectedToBeSmarts = expectedToBeSmarts;
        }

        private String testInput;
        private boolean expectedToBeSmarts;

        public String getTestInput() {
            return testInput;
        }

        public void setTestInput(String testInput) {
            this.testInput = testInput;
        }

        public boolean isExpectedToBeSmarts() {
            return expectedToBeSmarts;
        }

        public void setExpectedToBeSmarts(boolean expectedToBeSmarts) {
            this.expectedToBeSmarts = expectedToBeSmarts;
        }
    }

    @Test
    public void testSmartsEstimates() {

        for( smartsTestData data : generateTestData()) {
            boolean actual = CdkChemical2FactoryImpl.looksLikeSmarts(data.getTestInput());
            Assert.assertEquals(String.format("must detect that input %s is %b",
                    data.getTestInput(), data.isExpectedToBeSmarts()), data.isExpectedToBeSmarts(), actual);
        }
    }

    List<smartsTestData> generateTestData() {
        return Arrays.asList(smartsTestData.of("c1ccccc1", false),
                smartsTestData.of("c1ccccc1C#N", false),
                smartsTestData.of("[CX4]", true),
                smartsTestData.of("CCCOC[CX4]", true),
                smartsTestData.of("[CX3]=[OX1]", true),
                smartsTestData.of("NCCC[CX3]=[OX1]", true),
                smartsTestData.of("CC(C)(C)OC(=O)NC1CCN(CC1)c2c3cc(ccc3ncc2Cl)-c4cccc(C#N)c4O", false)
        );
    }
}
