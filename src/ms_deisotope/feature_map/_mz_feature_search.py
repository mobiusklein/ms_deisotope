from typing import Iterator, List, Sequence, TypeVar, Union

T = TypeVar("T")

def binary_search_with_flag(array, mz, error_tolerance=1e-5):
    lo = 0
    n = hi = len(array)
    while hi != lo:
        mid = (hi + lo) // 2
        x = array[mid]
        err = (x.mz - mz) / mz
        if abs(err) <= error_tolerance:
            i = mid - 1
            # Begin Sweep forward
            while i > 0:
                x = array[i]
                err = (x.mz - mz) / mz
                if abs(err) <= error_tolerance:
                    i -= 1
                    continue
                else:
                    break
            low_end = i
            i = mid + 1

            # Begin Sweep backward
            while i < n:
                x = array[i]
                err = (x.mz - mz) / mz
                if abs(err) <= error_tolerance:
                    i += 1
                    continue
                else:
                    break
            high_end = i
            return list(range(low_end, high_end)), True
        elif (hi - lo) == 1:
            return [mid], False
        elif err > 0:
            hi = mid
        elif err < 0:
            lo = mid
    return [0], False


def binary_search(array, mz, error_tolerance=1e-5):
    """Binary search an ordered array of objects with :attr:`mz`
    using a PPM error tolerance of `error_tolerance`

    Parameters
    ----------
    array : list
        An list of objects, sorted over :attr:`mz` in increasing order
    mz : float
        The mz to search for
    error_tolerance : float, optional
        The PPM error tolerance to use when deciding whether a match has been found

    Returns
    -------
    int:
        The index in `array` of the best match
    """
    lo = 0
    n = hi = len(array)
    while hi != lo:
        mid = (hi + lo) // 2
        x = array[mid]
        err = (x.mz - mz) / mz
        if abs(err) <= error_tolerance:
            best_index = mid
            best_error = abs(err)
            i = mid - 1
            while i >= 0:
                x = array[i]
                err = abs((x.mz - mz) / mz)
                if err < best_error:
                    best_error = err
                    best_index = i
                i -= 1

            i = mid + 1
            while i < n:
                x = array[i]
                err = abs((x.mz - mz) / mz)
                if err < best_error:
                    best_error = err
                    best_index = i
                i += 1
            return best_index
        elif (hi - lo) == 1:
            return None
        elif err > 0:
            hi = mid
        elif err < 0:
            lo = mid
    return 0


def binary_search_exact(array, mz):
    lo = 0
    hi = len(array)
    while hi != lo:
        mid = (hi + lo) // 2
        x = array[mid]
        err = (x.mz - mz)
        if err == 0:
            return mid
        elif (hi - lo) == 1:
            return mid
        elif err > 0:
            hi = mid
        else:
            lo = mid
    return 0


def search_sweep(array, mz, error_tolerance=1e-5):
    lo = 0
    n = hi = len(array)
    while hi != lo:
        mid = (hi + lo) // 2
        x = array[mid]
        err = (x.mz - mz) / mz
        if abs(err) <= error_tolerance:
            i = mid - 1
            while i >= 0:
                x = array[i]
                err = abs((x.mz - mz) / mz)
                if err > error_tolerance:
                    break
                i -= 1
            start = i + 1
            i = mid + 1
            while i < n:
                x = array[i]
                err = abs((x.mz - mz) / mz)
                if err > error_tolerance:
                    break
                i += 1
            end = i
            return (start, end)
        elif (hi - lo) == 1:
            return None
        elif err > 0:
            hi = mid
        elif err < 0:
            lo = mid
    return 0


def binary_search_with_flag_neutral(array, neutral_mass, error_tolerance=1e-5):
    lo = 0
    n = hi = len(array)
    while hi != lo:
        mid = (hi + lo) // 2
        x = array[mid]
        err = (x.neutral_mass - neutral_mass) / neutral_mass
        if abs(err) <= error_tolerance:
            i = mid - 1
            # Begin Sweep forward
            while i > 0:
                x = array[i]
                err = (x.neutral_mass - neutral_mass) / neutral_mass
                if abs(err) <= error_tolerance:
                    i -= 1
                    continue
                else:
                    break
            low_end = i
            i = mid + 1

            # Begin Sweep backward
            while i < n:
                x = array[i]
                err = (x.neutral_mass - neutral_mass) / neutral_mass
                if abs(err) <= error_tolerance:
                    i += 1
                    continue
                else:
                    break
            high_end = i
            return list(range(low_end, high_end)), True
        elif (hi - lo) == 1:
            return [mid], False
        elif err > 0:
            hi = mid
        elif err < 0:
            lo = mid
    return [0], False


def binary_search_neutral(array, neutral_mass, error_tolerance=1e-5):
    """Binary search an ordered array of objects with :attr:`neutral_mass`
    using a PPM error tolerance of `error_toler

    Parameters
    ----------
    array : list
        An list of objects, sorted over :attr:`neutral_mass` in increasing order
    neutral_mass : float
        The neutral_mass to search for
    error_tolerance : float, optional
        The PPM error tolerance to use when deciding whether a match has been found

    Returns
    -------
    int:
        The index in `array` of the best match
    """
    lo = 0
    n = hi = len(array)
    while hi != lo:
        mid = (hi + lo) // 2
        x = array[mid]
        err = (x.neutral_mass - neutral_mass) / neutral_mass
        if abs(err) <= error_tolerance:
            best_index = mid
            best_error = abs(err)
            i = mid - 1
            while i >= 0:
                x = array[i]
                err = abs((x.neutral_mass - neutral_mass) / neutral_mass)
                if err < best_error:
                    best_error = err
                    best_index = i
                i -= 1

            i = mid + 1
            while i < n:
                x = array[i]
                err = abs((x.neutral_mass - neutral_mass) / neutral_mass)
                if err < best_error:
                    best_error = err
                    best_index = i
                i += 1
            return best_index
        elif (hi - lo) == 1:
            return None
        elif err > 0:
            hi = mid
        elif err < 0:
            lo = mid
    return 0


def search_sweep_neutral(array, neutral_mass, error_tolerance=1e-5):
    lo = 0
    n = hi = len(array)
    while hi != lo:
        mid = (hi + lo) // 2
        x = array[mid]
        err = (x.neutral_mass - neutral_mass) / neutral_mass
        if abs(err) <= error_tolerance:
            i = mid - 1
            while i >= 0:
                x = array[i]
                err = abs((x.neutral_mass - neutral_mass) / neutral_mass)
                if err > error_tolerance:
                    break
                i -= 1
            start = i + 1
            i = mid + 1
            while i < n:
                x = array[i]
                err = abs((x.neutral_mass - neutral_mass) / neutral_mass)
                if err > error_tolerance:
                    break
                i += 1
            end = i
            return (start, end)
        elif (hi - lo) == 1:
            return None
        elif err > 0:
            hi = mid
        elif err < 0:
            lo = mid
    return 0, 0


def binary_search_exact_neutral(array, neutral_mass):
    lo = 0
    hi = len(array)
    while hi != lo:
        mid = (hi + lo) // 2
        x = array[mid]
        err = (x.neutral_mass - neutral_mass)
        if err == 0:
            return mid
        elif (hi - lo) == 1:
            return mid
        elif err > 0:
            hi = mid
        else:
            lo = mid
    return 0


class _FeatureIndex(Sequence[T]):
    def __iter__(self) -> Iterator[T]:
        return iter(self.features)

    def __getitem__(self, i) -> Union[T, List[T]]:
        return self.features[i]

    def __len__(self):
        return len(self.features)

    def __nonzero__(self):
        return bool(self.features)

    def __bool__(self):
        return bool(self.features)

    def __repr__(self):
        return "{self.__class__.__name__}(<{size} features>)".format(self=self, size=len(self))


class MZIndex(_FeatureIndex[T]):
    def __init__(self, features):
        self.features = sorted(features, key=lambda x: x.mz)

    def find_all(self, mz: float, error_tolerance: float=2e-5) -> List[T]:
        bounds = search_sweep(self.features, mz, error_tolerance)
        if bounds is not None:
            lo, hi = bounds
            return self[lo:hi]
        else:
            return []


class NeutralMassIndex(_FeatureIndex[T]):
    def __init__(self, features):
        self.features = sorted(features, key=lambda x: x.neutral_mass)

    def find_all(self, mass: float, error_tolerance: float=2e-5) -> List[T]:
        bounds = search_sweep_neutral(self.features, mass, error_tolerance)
        if bounds is not None:
            lo, hi = bounds
            return self[lo:hi]
        else:
            return []


try:
    from ms_deisotope._c.feature_map.feature_map import binary_search_with_flag
except ImportError:
    pass
