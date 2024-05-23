import gov.nih.ncats.molwitch.Chemical;
import gov.nih.ncats.molwitch.cdk.CdkChemicalImpl;
import org.apache.commons.io.IOUtils;
import org.junit.Assert;
import org.junit.Test;

import java.io.IOException;
import java.util.List;

public class TestSGroupWarnings {

    @Test
    public void cleanStyreneNoWarnings() throws IOException {
        String molfileText = IOUtils.toString(
                this.getClass().getResourceAsStream("mols/polystyrene_clean.mol"),
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
                this.getClass().getResourceAsStream("mols/polyethylstyrene_really_messy.mol"),
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
                this.getClass().getResourceAsStream("mols/polyethylstyrene_messy.mol"),
                "UTF-8"
        );
        Chemical c1=Chemical.parse(molfileText);
        CdkChemicalImpl impl = (CdkChemicalImpl)c1.getImpl();
        List<String> warnings = impl.getSGroupWarnings();
        Assert.assertFalse(warnings.isEmpty());
    }
}

