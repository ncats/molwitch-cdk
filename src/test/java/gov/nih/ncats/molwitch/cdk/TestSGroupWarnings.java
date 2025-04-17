package gov.nih.ncats.molwitch.cdk;

import gov.nih.ncats.molwitch.Chemical;
import org.apache.commons.io.IOUtils;
import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;
import java.util.List;

public class TestSGroupWarnings {

    @Test
    public void cleanStyreneNoWarnings() throws IOException {
        String molfileText = IOUtils.toString(
                this.getClass().getResourceAsStream("/mols/polystyrene_clean.mol"),
                "UTF-8"
        );
        Chemical c1=Chemical.parse(molfileText);
        CdkChemicalImpl impl = (CdkChemicalImpl)c1.getImpl();
        List<String> warnings = impl.getSGroupWarnings();
        Assert.assertTrue(warnings.isEmpty());
    }

    @Test
    public void styreneWarnings() throws IOException {
        String molfileText = IOUtils.toString(
                this.getClass().getResourceAsStream("/mols/polyethylstyrene_really_messy.mol"),
                "UTF-8"
        );
        Chemical c1=Chemical.parse(molfileText);
        CdkChemicalImpl impl = (CdkChemicalImpl)c1.getImpl();
        List<String> warnings = impl.getSGroupWarnings();
        Assert.assertFalse(warnings.isEmpty());
    }

    @Test
    public void styreneWarnings2() throws IOException {
        String molfileText = IOUtils.toString(
                this.getClass().getResourceAsStream("/mols/polyethylstyrene_messy.mol"),
                "UTF-8"
        );
        Chemical c1=Chemical.parse(molfileText);
        CdkChemicalImpl impl = (CdkChemicalImpl)c1.getImpl();
        List<String> warnings = impl.getSGroupWarnings();
        Assert.assertFalse(warnings.isEmpty());
    }

    @Test
    public void styreneWarnings3() throws IOException {
        String molfileText = IOUtils.toString(
                this.getClass().getResourceAsStream("/mols/another.mol"),
                "UTF-8"
        );
        Chemical c1=Chemical.parse(molfileText);
        CdkChemicalImpl impl = (CdkChemicalImpl)c1.getImpl();
        List<String> warnings = impl.getSGroupWarnings();
        Assert.assertFalse(warnings.isEmpty());
    }

    @Test
    public void doublePolymerNoWarnings() throws IOException {
        String molfileText = IOUtils.toString(
                this.getClass().getResourceAsStream("/mols/double_polymer1.mol"),
                "UTF-8"
        );
        Chemical c1=Chemical.parse(molfileText);
        CdkChemicalImpl impl = (CdkChemicalImpl)c1.getImpl();
        List<String> warnings = impl.getSGroupWarnings();
        Assert.assertTrue(warnings.isEmpty());
    }

    @Test
    public void doublePolymerWarnings() throws IOException {
        String molfileText = IOUtils.toString(
                this.getClass().getResourceAsStream("/mols/double_polymer2.mol"),
                "UTF-8"
        );
        Chemical c1=Chemical.parse(molfileText);
        CdkChemicalImpl impl = (CdkChemicalImpl)c1.getImpl();
        List<String> warnings = impl.getSGroupWarnings();
        Assert.assertFalse(warnings.isEmpty());
    }

    @Test
    public void multipleGroupNoWarnings() throws IOException {
        String molfileText = IOUtils.toString(
                this.getClass().getResourceAsStream("/mols/calcium benzoate monohydrate.mol"),
                "UTF-8"
        );
        Chemical c1=Chemical.parse(molfileText);
        CdkChemicalImpl impl = (CdkChemicalImpl)c1.getImpl();
        List<String> warnings = impl.getSGroupWarnings();
        Assert.assertTrue(warnings.isEmpty());
    }

    @Test
    public void multipleGroupWarnings() throws IOException {
        String molfileText = IOUtils.toString(
                this.getClass().getResourceAsStream("/mols/calcium benzoate monohydrate-18 confusing.mol"),
                "UTF-8"
        );
        Chemical c1=Chemical.parse(molfileText);
        CdkChemicalImpl impl = (CdkChemicalImpl)c1.getImpl();
        List<String> warnings = impl.getSGroupWarnings();
        Assert.assertFalse(warnings.isEmpty());
    }

}
