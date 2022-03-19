import numpy as np

def scale(arr, mn, mx, arr_min = None, arr_max = None):
    if len(arr)>1:
        arr = np.array(arr)

        if arr_min is None:
            arr_min = arr.min()
        else:
            arr = np.where(arr<arr_min, arr_min, arr)

        if arr_max is None:
            arr_max = arr.max()
        else:
            arr = np.where(arr>arr_max, arr_max, arr)

        return (mx-mn)*((arr-arr_min)/(arr_max-arr_min))+mn
    elif len(arr)==1:
        return np.array([mx])
    else:
        return arr