**W**hich **I**nternal **T**oolkit for **CH**emicals

MolWitch implementation using [CDK](https://cdk.github.io/)


# Available on Maven Central
Usually, one needs to add 2 dependencies:
This adds the API.
```
<dependency>
  <groupId>gov.nih.ncats</groupId>
  <artifactId>molwitch</artifactId>
  <version>0.6.1</version>
</dependency>
```

There also needs to be a molwitch implementation

To add CDK molwitch implementation (this version uses CDK 2.6):
```
<dependency>
        <groupId>gov.nih.ncats</groupId>
        <artifactId>molwitch-cdk</artifactId>
        <version>1.0.8</version>
</dependency>
```

## Current API Contract Compliance
Results from running the latest code on Molwitch-cdk using the [API Contract](https://github.com/ncats/molwitch-apitests)

| Feature | Compliance Level | Comments|
 | ------ | ---------- | ---------- |
| Extended Tetrahedral| PARTIALLY |  |
| Fingerprint| FULLY |  |
| fullInchi | FULLY ( 998 ) |  |
| fullInchi | NOT_COMPLIANT ( 2 ) |  |
| Valence Error | FULLY ( 1 ) |  |
| Valence Error | PARTIALLY ( 1 ) |  |
| Valence Error | NOT_COMPLIANT ( 1 ) |  Hypervalent Hydrogen Incorrect Valence |
| parse mol wierd parity| FULLY |  |
| inchiKey | FULLY ( 998 ) |  |
| inchiKey | NOT_COMPLIANT ( 2 ) |  |
| Remove Non Descript Hydrogens| FULLY |  |
| Inchi| FULLY |  |
| Default Fingerprinter| FULLY |  |
| Mol Parser| FULLY |  |
| Create Chemical| FULLY |  |
| MolSearcher| FULLY |  |
| Write Mol| FULLY |  |
| Problematic Smiles| FULLY |  |
| Chemical Source| FULLY |  |
| Clone Chemical| FULLY |  |
| Atom Alias| FULLY |  |
| mol parser unknown format| FULLY |  |
| Hetero Atom Tetrahedral| FULLY |  |
| Tetrahedral| FULLY |  |
| Atom Path Traversal| FULLY |  |
| Atom Coords| FULLY |  |
| Isotopes| FULLY |  |
| Cis/Trans| FULLY |  |
