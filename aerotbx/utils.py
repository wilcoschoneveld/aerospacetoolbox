import scipy as sp
from types import ListType, IntType, LongType, FloatType, TupleType

#convert any items to numpy arrays
def to_ndarray(item):
    return type(item), sp.array(item, sp.float64, ndmin=1)

#convert any numpy array back to original item type
def from_ndarray(itemtype, *items):
    if itemtype in [IntType, LongType, FloatType]:
        v = tuple(i.flat[0] for i in items)
    elif itemtype is ListType:
        v = tuple(i.tolist() for i in items)
    elif itemtype is TupleType:
        v = tuple(tuple(i.tolist()) for i in items)
    elif itemtype is sp.matrix:
        v = tuple(sp.matrix(i, copy=False) for i in items)
    else:
        v = tuple(i for i in items)
    return v if len(items) > 1 else v[0]
