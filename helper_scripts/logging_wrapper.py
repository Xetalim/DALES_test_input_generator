# Source - https://stackoverflow.com/a
# Posted by HanSooloo, modified by community. See post 'Timeline' for change history
# Retrieved 2025-12-16, License - CC BY-SA 4.0
import logging


# At the beginning of every .py file in the project
def logwrap(fn):
    from functools import wraps
    import inspect

    @wraps(fn)
    def wrapper(*args, **kwargs):
        logger = logging.getLogger(fn.__module__)
        logger.debug("Entering %s", fn.__name__)
        try:
            out = fn(*args, **kwargs)
        except Exception as e:
            logger.error("Function %s failed with exception %s: ", fn.__name__, e.args)
            raise e

        logger.debug("Exiting %s", fn.__name__)
        # Return the return value
        return out

    return wrapper
