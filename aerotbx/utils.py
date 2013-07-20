"""Aerospace Toolbox / utils.py

Utility functions used by the entire package.

"""

import scipy as sp
from types import ListType, IntType, LongType, FloatType, TupleType

def to_ndarray(item):
    """convert any item to a numpy ndarray"""
    
    return type(item), sp.array(item, sp.float64, ndmin=1)

def from_ndarray(itemtype, *items):
    """convert any numpy array back to its original item type"""

    if itemtype in [IntType, LongType, FloatType]:
        val = tuple(i.flat[0] for i in items)
    elif itemtype is ListType:
        val = tuple(i.tolist() for i in items)
    elif itemtype is TupleType:
        val = tuple(tuple(i.tolist()) for i in items)
    elif itemtype is sp.matrix:
        val = tuple(sp.matrix(i, copy=False) for i in items)
    else:
        val = tuple(i for i in items)
    return val if len(items) > 1 else val[0]