import scipy as sp
from types import ListType, IntType, LongType, FloatType, TupleType

def to_ndarray(item):
    return type(item), sp.array(item, sp.float64, ndmin=1)

def from_ndarray(itemtype, itemarray):
    if itemtype in [IntType, LongType, FloatType]:
        return itemarray.flat[0]
    elif itemtype is ListType:
        return itemarray.tolist()
    elif itemtype is TupleType:
        return tuple(itemarray.tolist())
    elif itemtype is sp.matrix:
        return sp.matrix(itemarray, copy=False)
    else:
        return itemarray
