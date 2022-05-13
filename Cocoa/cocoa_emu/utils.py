def get_params_from_sample(sample, labels):
    """
    Format arrays into cocoa params
    """
    assert len(sample)==len(labels), "Length of the labels not equal to the length of samples"
    params = {}
    for i, label in enumerate(labels):
        param_i = sample[i]
        params[label] = param_i
    return params

def get_params_list(samples, labels):
    params_list = []
    for i in range(len(samples)):
        params = get_params_from_sample(samples[i], labels)
        params_list.append(params)
    return params_list

def get_params_from_lhs_sample(unit_sample, lhs_prior):
    """
    Format unit LHS arrays into cocoa params
    """
    assert len(unit_sample)==len(lhs_prior), "Length of the labels not equal to the length of samples"
    params = {}
    for i, label in enumerate(lhs_prior):
        lhs_min = lhs_prior[label]['min']
        lhs_max = lhs_prior[label]['max']
        param_i = lhs_min + (lhs_max - lhs_min) * unit_sample[i]
        params[label] = param_i
    return params

def get_lhs_params_list(samples, lhs_prior):
    params_list = []
    for i in range(len(samples)):
        params = get_params_from_lhs_sample(samples[i], lhs_prior)
        params_list.append(params)
    return params_list