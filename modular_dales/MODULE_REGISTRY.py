MODULE_REGISTRY = {}
SINGLETON_REGISTRY = {}  # Keep track of singleton module names
SPECIAL_SERIALIZING_REGISTRY = (
    {}
)  # For modules that require special handling during serialization


def register_module(cls, singleton=False):
    """Class decorator to register simulation modules."""
    MODULE_REGISTRY[cls.__name__] = cls
    return cls


def register_singleton(cls):
    """Class decorator to register a singleton simulation module."""
    SINGLETON_REGISTRY[cls.__name__] = cls
    return cls


def register_special_serializing(cls):
    """Class decorator to register a module that requires special handling during serialization."""
    SPECIAL_SERIALIZING_REGISTRY[cls.__name__] = cls
    return cls
