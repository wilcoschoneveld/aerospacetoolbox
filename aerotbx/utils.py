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

def from_ndarray2(itemtype, *items):
    if itemtype in [IntType, LongType, FloatType]:
        return [i.flat[0] for i in items]
    elif itemtype is ListType:
        return [i.tolist() for i in items]
    elif itemtype is TupleType:
        return [tuple(i.tolist()) for i in items]
    elif itemtype is sp.matrix:
        return [sp.matrix(i, copy=False) for i in items]
    else:
        return [i for i in items]
