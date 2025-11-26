import logging

class LoggerMixin:
    def __init__(self, logger_level=None):
        if logger_level is None:
            logger_level = logging.INFO
        self.logger = logging.getLogger(self.__class__.__name__)
        if not self.logger.hasHandlers():
            self.logger.setLevel(logger_level)


class DataBundle(LoggerMixin):
    """
    A container for holding all data needed by plot components.

    Attributes are added dynamically by the loaders
    depending on the data needed by the plot components.
    """
    def __init__(self, **kwargs):
        object.__setattr__(self, '_initialized_attributes', set())
        # Add additional data fields here (not strictly necessary)

        for key, value in kwargs.items():
            setattr(self, key, value)

        self.processed_loaders = []
        super().__init__()
    
    def __getattr__(self, name):
        loaders = ', '.join([loader.__name__ for loader in self.processed_loaders])
        initialized_attributes = ', '.join(self._initialized_attributes)
        raise AttributeError(f"Data loaded with {loaders} is missing attribute '{name}'. Available attributes: {initialized_attributes}")
    
    def __setattr__(self, name, value):
        if name != 'logger' and name in self._initialized_attributes:
            self.logger.warning(f"Warning: Attribute '{name}' is overridden.")
        self._initialized_attributes.add(name)
        super().__setattr__(name, value)

    def __repr__(self):
        display_attrs = ', '.join(
            [
                x for x in self._initialized_attributes 
                if x not in ('logger', 'processed_loaders', 'interval')
            ]
        )
        return f"DataBundle({display_attrs})"

    def copy(self):
        """
        Create a shallow copy of the data bundle.
        """
        data = DataBundle()
        for attr in self._initialized_attributes:
            setattr(data, attr, getattr(self, attr))
        for loader in self.processed_loaders:
            data.processed_loaders.append(loader)
        return data


class RequiredArgument:
    def __repr__(self):
        return 'Required loader arg'
