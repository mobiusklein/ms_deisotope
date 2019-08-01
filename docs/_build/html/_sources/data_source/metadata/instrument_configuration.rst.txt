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


    .. data:: ionization_types

        .. exec::

            from ms_deisotope.data_source.metadata.instrument_components import ionization_types
            from rst_table import as_rest_table

            rows = [("Ionization Type", "ID", "Category")]
            for ionization_type in ionization_types:
                rows.append((ionization_type.name, ionization_type.id, ionization_type.specialization[0]))
            print(as_rest_table(rows))


    .. data:: detector_types

        .. exec::

            from ms_deisotope.data_source.metadata.instrument_components import detector_types
            from rst_table import as_rest_table

            rows = [("Detector Type", "ID", "Category")]
            for detector_type in detector_types:
                rows.append((detector_type.name, detector_type.id, detector_type.specialization[0]))
            print(as_rest_table(rows))


    .. data:: analyzer_types

        .. exec::

            from ms_deisotope.data_source.metadata.instrument_components import analyzer_types
            from rst_table import as_rest_table

            rows = [("Analyzer Type", "ID", "Category")]
            for analyzer_type in analyzer_types:
                rows.append((analyzer_type.name, analyzer_type.id, analyzer_type.specialization[0]))
            print(as_rest_table(rows))


    .. data:: inlet_types

        .. exec::

            from ms_deisotope.data_source.metadata.instrument_components import inlet_types
            from rst_table import as_rest_table

            rows = [("Inlet Type", "ID", "Category")]
            for inlet_type in inlet_types:
                rows.append((inlet_type.name, inlet_type.id, inlet_type.specialization[0]))
            print(as_rest_table(rows))
