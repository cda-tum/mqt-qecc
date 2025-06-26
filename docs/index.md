# Welcome to QECC's documentation

QECC is a tool for quantum error correcting codes developed as part of the _{doc}`Munich Quantum Toolkit (MQT) <mqt:index>`_ [^1].

We recommend you to start with the {doc}`installation instructions <Installation>`.
Then proceed to read the {doc}`reference documentation <api/mqt/qecc/index>`.
If you are interested in the theory behind QECC, have a look at the publications in the {doc}`publication list <Publications>`.

We appreciate any feedback and contributions to the project. If you want to contribute, you can find more information in
the {doc}`Contribution <Contributing>` guide. If you are having trouble with the installation or the usage of QECC,
please let us know at our {doc}`Support <Support>` page or by reaching out to us at
[quantum.cda@xcit.tum.de](mailto:quantum.cda@xcit.tum.de).

[^1]:
    The _[Munich Quantum Toolkit (MQT)](https://mqt.readthedocs.io)_ is a collection of software tools for quantum computing developed by the [Chair for Design Automation](https://www.cda.cit.tum.de/) at the [Technical University of Munich](https://www.tum.de/) as well as the [Munich Quantum Software Company (MQSC)](https://munichquantum.software).
    Among others, it is part of the [Munich Quantum Software Stack (MQSS)](https://www.munich-quantum-valley.de/research/research-areas/mqss) ecosystem, which is being developed as part of the [Munich Quantum Valley (MQV)](https://www.munich-quantum-valley.de) initiative.

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
