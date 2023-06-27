from cobaya.likelihoods.base_classes import planck_clik
from cobaya.log import LoggedError, get_logger

class TEEE(planck_clik.PlanckClik): 
    pass

_planck_get_data_path = planck_clik.get_data_path

def get_data_path(name):
    log = get_logger(name)
    if "spt" not in name.lower():
        return _planck_get_data_path(name)
    log.info("override default get_data_path from %s"%(_planck_get_data_path.__module__))
    return "spt_data"

planck_clik.get_data_path = get_data_path
