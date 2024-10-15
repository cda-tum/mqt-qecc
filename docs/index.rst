Welcome to QECC's documentation
================================

QECC is a tool for quantum error correcting codes developed as part of the `Munich Quantum Toolkit (MQT) <https://mqt.readthedocs.io>`_ by the `Chair for Design Automation <https://www.cda.cit.tum.de/>`_ at the `Technical University of Munich <https://www.tum.de>`_. It builds upon `MQT Core <https://github.com/cda-tum/mqt-core>`_, which forms the backbone of the MQT.

We recommend you to start with the :doc:`installation instructions <Installation>`.
Then proceed to read the :doc:`reference documentation <library/Library>`.
If you are interested in the theory behind QECC, have a look at the publications in the :doc:`publication list <Publications>`.

We appreciate any feedback and contributions to the project. If you want to contribute, you can find more information in
the :doc:`Contribution <Contributing>` guide. If you are having trouble with the installation or the usage of QECC,
please let us know at our :doc:`Support <Support>` page or by reaching out to us at
`quantum.cda@xcit.tum.de <mailto:quantum.cda@xcit.tum.de>`_.

----

 .. toctree::
    :hidden:

    self

 .. toctree::
    :maxdepth: 2
    :caption: User Guide
    :glob:

    Installation
    LightsOutDecoder
    EccFramework
    StatePrep
    Publications

 .. toctree::
    :maxdepth: 2
    :caption: Developers
    :glob:

    Contributing
    DevelopmentGuide
    Support

 .. toctree::
    :maxdepth: 6
    :caption: API Reference
    :glob:

    library/Library
