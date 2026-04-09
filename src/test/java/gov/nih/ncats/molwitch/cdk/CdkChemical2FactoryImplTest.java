package gov.nih.ncats.molwitch.cdk;

import org.junit.Assert;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.Arguments;
import org.junit.jupiter.params.provider.MethodSource;

import java.util.stream.Stream;

public class CdkChemical2FactoryImplTest {

    @ParameterizedTest
    @MethodSource("inputData")
    void testMultipleSmarts(String smiles, boolean expectedToBeSmarts) {
        boolean result = CdkChemical2FactoryImpl.looksLikeSmarts(smiles);
        Assert.assertEquals(expectedToBeSmarts, result);
    }

    private static Stream<Arguments> inputData() {
        return Stream.of(
                Arguments.of("[C-]#[N+]C(F)(F)F", false),
                Arguments.of("ccccccC(=O)O", false),
                Arguments.of("*-N(-[#1])(-[#1])", true),
                Arguments.of("[CX3H1](=O)[#6]", true),
                Arguments.of("[CX3](=[OX1])[F,Cl,Br,I]", true),
                Arguments.of("[C]-S-[C]", false),
                Arguments.of("N(=NN1)C=N1", false),
                Arguments.of("n(nn1)cn1", false),
                Arguments.of("c1[nH]nnn1", false),
                Arguments.of("c1ccccc1C#N", false),
                Arguments.of("[CX4]", true),
                Arguments.of("CCCOC[CX4]", true),
                Arguments.of("[CX3]=[OX1]", true),
                Arguments.of("NCCC[CX3]=[OX1]", true),
                Arguments.of("CC(C)(C)OC(=O)NC1CCN(CC1)c2c3cc(ccc3ncc2Cl)-c4cccc(C#N)c4O", false));
        }
    }
