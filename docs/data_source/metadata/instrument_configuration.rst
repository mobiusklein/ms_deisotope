Instrument Configuration
------------------------

A mass spectrometer may be configured to operate in multiple modes during
a single data acquisition run. A program may need to know which configurations
were available, and formats like :title-reference:`mzML` are able record this
information. To retrieve the configuration used for a particular scan, see
:attr:`ms_deisotope.data_source.common.Scan.instrument_configuration`

.. automodule:: ms_deisotope.data_source.metadata.instrument_components

    .. autoclass:: Component
        :members:
        
        .. automethod:: is_a

    .. autoclass:: ComponentGroup
        :members:

    .. autoclass:: InstrumentInformation
        :members: