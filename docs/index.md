# Welcome to QECC's documentation

QECC is a tool for quantum error correcting codes developed as part of the [Munich Quantum Toolkit (MQT)](https://mqt.readthedocs.io) by the [Chair for Design Automation](https://www.cda.cit.tum.de/) at the [Technical University of Munich](https://www.tum.de).

We recommend you to start with the {doc}`installation instructions <Installation>`.
Then proceed to read the {doc}`reference documentation <api/mqt/qecc/index>`.
If you are interested in the theory behind QECC, have a look at the publications in the {doc}`publication list <Publications>`.

We appreciate any feedback and contributions to the project. If you want to contribute, you can find more information in
the {doc}`Contribution <Contributing>` guide. If you are having trouble with the installation or the usage of QECC,
please let us know at our {doc}`Support <Support>` page or by reaching out to us at
[quantum.cda@xcit.tum.de](mailto:quantum.cda@xcit.tum.de).

---

```{toctree}
:hidden:

self
```

```{toctree}
:maxdepth: 1
:caption: User Guide

Installation
LightsOutDecoder
StatePrep
CatStates
Encoders
AnalogInfo
Publications
```

```{toctree}
:caption: Developers
:glob:
:maxdepth: 1

Contributing
DevelopmentGuide
Support
```

```{toctree}
:maxdepth: 3
:caption: API Reference

api/mqt/qecc/index
```
