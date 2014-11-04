
"""Yes I know it's not Pythonic to type check, but for pylinac, I don't see many alternatives."""

# The following is adapted from: https://wiki.python.org/moin/PythonDecoratorLibrary#Type_Enforcement_.28accepts.2Freturns.29
def type_accept(*types):
    """Function decorator enforcing input types.

    :param types: expected types (e.g. str, int).
    """

    def type_decor(func):  # func == function to decorate
        def new_func(self, *args, **kwargs):  # the decorated function's arguments'
            if len(args) != 0 or len(kwargs) != 0:
                # For each argument passed, check it against the type assertion
                for arg_list in [args, kwargs.values()]:
                    for idx, arg in enumerate(arg_list):
                        if type(arg) != types[idx]:  # for single type comparisons
                            if len(types[idx]) > 1 and (type(arg) in types[idx]):  # for multiple type comparisons
                                pass
                            else:
                                raise TypeError("{} was expected to be of type {} but was of type {}".format(arg, str(types[idx]).split("'")[1],
                                                                                                             str(type(arg)).split("'")[1]))
            return func(self, *args, **kwargs)
        return new_func
    return type_decor

def value_accept(*val_accept):
    """Function decorator enforcing input types.

    :param val_accept: expected values in format (lower range, upper range). E.g. (1, 10) means accept a value between
        1 and 10.
    """
    def decorator(func):  # func == function to decorate
        def new_func(self, *args, **kwargs):  # the decorated function's arguments'
            # For each argument passed, check it against the type assertion
            if len(args) != 0 or len(kwargs) != 0:
                for arg_list in [args, list(kwargs.values())]:
                    for idx, val in enumerate(arg_list):
                        if arg_list[idx] == None:
                            pass
                        elif type(arg_list[idx]) in (float, int):
                            if arg_list[idx] < val_accept[idx][0] and arg_list[idx] > val_accept[idx][1]:  #
                        #  for
                        # upper/lower
                        # range comparisons
                                raise ValueError("{:f} needs to be between {:f} and {:f}".format(arg_list[idx], val[0], val[1]))
                        elif type(arg_list[idx]) == str:
                            if arg_list[idx] not in val:
                                raise KeyError("{} not one of these options: {}".format(arg_list[idx], val))
            return func(self, *args, **kwargs)
        return new_func
    return decorator