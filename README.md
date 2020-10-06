**W**hich **I**nternal **T**oolkit for **CH**emicals

MolWitch implementation using [CDK](https://cdk.github.io/)


# Available on Maven Central
Usually, one needs to add 2 dependencies:
This adds the API.
```
<dependency>
  <groupId>gov.nih.ncats</groupId>
  <artifactId>molwitch</artifactId>
  <version>0.6.0</version>
</dependency>
```

There also needs to be a molwitch implementation

To add CDK molwitch implementation (NOTE: currently a snapshot until CDK 2.4 is released) :
```
<dependency>
        <groupId>gov.nih.ncats</groupId>
        <artifactId>molwitch-cdk</artifactId>
        <version>1.0.4-SNAPSHOT</version>
</dependency>
```

## Current API Contract Compliance
Results from running the latest code on Molwitch-jchem using the [API Contract](https://github.com/ncats/molwitch-apitests)

| Feature | Compliance Level | Comments|
 | ------ | ---------- | ---------- |
| Extended Tetrahedral| FULLY |  |
| Fingerprint| FULLY |  |
| fullInchi| FULLY |  |
| Valence Error | FULLY ( 1 ) |  |
| Valence Error | PARTIALLY ( 1 ) |  |
| Valence Error | NOT_COMPLIANT ( 1 ) |  Hypervalent Hydrogen Incorrect Valence |
| parse mol wierd parity| FULLY |  |
| inchiKey| FULLY |  |
| Remove Non Descript Hydrogens| FULLY |  |
| Inchi| FULLY |  |
| Default Fingerprinter| FULLY |  |
| Mol Parser| FULLY |  |
| Create Chemical| FULLY |  |
| MolSearcher| FULLY |  |
| Write Mol| FULLY |  |
| Problematic Smiles| FULLY |  |
| Chemical Source| PARTIALLY |  |
| Clone Chemical| FULLY |  |
| Atom Alias| FULLY |  |
| mol parser unknown format| FULLY |  |
| Tetrahedral| FULLY |  |
| Atom Path Traversal| FULLY |  |
| Atom Coords| FULLY |  |
| Cis/Trans| FULLY |  |