"""
This module contains utility functions for working with strings
"""

def grepl(s, l, index=False, ignorecase=False):
    """
        Return all elements of list, l which contains string, s. If index is True
        return indices instead of elements.
    """
    from collections import deque
    # TODO: add functionality for ignorecase
    out = deque()
    try:
        assert(isinstance(l, list))
    except AssertionError:
        raise ValueError('Input variable is not a list')
    for i, e in enumerate(l):
        if s in e:
            out.append(i if index else e)
    return out
# END
