import numpy as np


def checkinput(Mi, zi, ci, verbose=None):
    """ Check and convert any input scalar or array to numpy array """

    # How many halo redshifts provided?
    if verbose:
        print('How many redshift output provided?')
    zi = convert_to_np_array(zi)

    if verbose:
        print('{} redshift outputs'.format(len(zi)))
        print('How many halo masses provided?')

    Mi = convert_to_np_array(Mi)
    Mi_shape = Mi.shape
    if verbose:
        print(Mi_shape[0])
    if Mi_shape[1] == len(zi):
        if verbose:
            print('Redshift & MH array dimensions are ok')
    if Mi_shape[1] != len(zi):
        if verbose:
            print('Redshift & MH array dimensions must match')
        raise Exception('Redshift & MH array dimensions must match. Got instead {} and {}'.format(Mi_shape[1], len(zi)))

    # How many halo concentrations provided?
    if verbose:
        print('How many halo concentrations provided?')
    ci = convert_to_np_array(ci)
    if verbose:
        print(len(ci))
    if Mi_shape[0] == len(ci):
        if verbose:
            print('Concentration & MH(z=0) array dimensions are ok')
    else:
        if verbose:
            print('Concentration & MH(z=0) array dimensions must match')
        raise Exception('Concentration & MH(z=0) array dimensions must match. \
            Got instead {} and {}'.format(Mi_shape[0], len(ci)))
    return Mi, zi, ci


def convert_to_np_array(x):
    """ Converts the given input to a numpy array. It takes into account that the input can already be a numpy array,
    an iterable or a number.

    :param x: A numpy array, iterable or number
    :return: numpy array
    """

    if not isinstance(x, np.ndarray):
        try:
            x = np.array(x)
        except TypeError:
            x = np.array([x])
        except Exception as e:
            raise e

    return x


def bestfit(x, gamma):
    """ Best-fitting function, with gamma unkown """
    alpha_1 = x[0]
    alpha_2 = x[1]
    x = x[2:len(x)]
    alpha = alpha_1-gamma*alpha_2
    f = alpha*x/3.+gamma*(np.exp(x/3.)-1.0)
    return f
