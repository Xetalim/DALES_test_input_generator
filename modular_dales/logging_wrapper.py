# Source - https://stackoverflow.com/a
# Posted by HanSooloo, modified by community. See post 'Timeline' for change history
# Retrieved 2025-12-16, License - CC BY-SA 4.0
import logging
import pathlib
import yaml


# At the beginning of every .py file in the project
def logwrap(fn):
    from functools import wraps
    import inspect

    @wraps(fn)
    def wrapper(*args, **kwargs):
        logger = logging.getLogger(fn.__module__)
        module = None
        for arg in args:
            if hasattr(arg, "module_name"):
                module = arg.module_name
                break
        else:
            for arg in kwargs.values():
                if hasattr(arg, "module_name"):
                    module = arg.module_name
        if module:
            logger.debug("Entering %s (module: %s)", fn.__name__, module)
        else:
            logger.debug("Entering %s", fn.__name__)
        try:
            out = fn(*args, **kwargs)
        except Exception as e:
            logger.error("Function %s failed with exception %s: ", fn.__name__, e.args)
            raise e
        if module:
            logger.debug("Exiting %s (module: %s)", fn.__name__, module)
        else:
            logger.debug("Exiting %s", fn.__name__)
        # Return the return value
        return out

    return wrapper


def setup_logging(config_path="logging.yaml"):
    import logging.config

    path = pathlib.Path(config_path)
    if path.exists():
        with path.open() as f:
            user_cfg = yaml.safe_load(f)

        logging.config.dictConfig(user_cfg)
    else:
        logging.basicConfig()
