Controlled Vocabulary Interface
-------------------------------

This module provides a limited interface to the :title-reference:`HUPO PSI-MS` controlled
vocabulary through the :mod:`psims` module. This module is used heavily internally to generate
static data lists.

.. automodule:: ms_deisotope.data_source.metadata.cv

    .. autoclass:: Term
        :members: is_a

    .. data:: cv_psims

        A lazily loaded mapping over the PSI-MS controlled vocabulary. It may be queried
        using names, aliases, or accession numbers, and the term graph may be traversed.
