# A selection of functions to compute RMS

import numpy as np
from . import utilities


def L2_on_vector(resid_vector, sigma_vector):
    """
    Take the RMS on a simple vector. Chi2 is computed with associated sigmas. Pure math function.
    """
    rms_mm = np.multiply(np.sqrt(np.mean(np.square(resid_vector))), 1000);
    added_sum = 0;
    for i in range(len(resid_vector)):
        chi2 = np.square(resid_vector[i]) / np.square(sigma_vector[i]);
        added_sum = added_sum + chi2
    reported_chi2 = np.sqrt(added_sum / len(sigma_vector));
    return rms_mm, reported_chi2;


def L1_on_vector(resid_vector, sigma_vector):
    """
    Take the L1 norm on a simple vector. Chi2 is computed with associated sigmas. Pure math function.
    """
    rms_mm = np.multiply(np.mean(np.abs(resid_vector)), 1000);
    added_sum = 0;
    for i in range(len(resid_vector)):
        normalized = np.abs(resid_vector[i]) / sigma_vector[i];
        added_sum = added_sum + normalized
    reported_normalized = added_sum / len(sigma_vector);
    return rms_mm, reported_normalized;


def filter_disp_pts_to_vector(disp_points, components=("E", "N", "U")):
    """
    Filter disp_points by type and then expand them into a vector of obs and a matching vector of sigmas.
    """
    all_misfits, all_sigmas = [], [];
    for item in disp_points:
        if "E" in components:
            all_misfits.append(item.dE_obs);
            all_sigmas.append(item.Se_obs);
        if "N" in components:
            all_misfits.append(item.dN_obs);
            all_sigmas.append(item.Sn_obs);
        if "U" in components:
            all_misfits.append(item.dU_obs);
            all_sigmas.append(item.Su_obs);
    return all_misfits, all_sigmas;


def obs_vs_model_L2_aggregate(obs_disp_points, model_disp_points):
    """
    L2 norm on observed vs modeled disp_point objects. Unpacks into vectors, then computes L2 norm.
    """
    resid = utilities.subtract_disp_points(obs_disp_points, model_disp_points);   # make residual points
    res_cgps = filter_disp_pts_to_vector(utilities.filter_to_meas_type(resid, "continuous"), ("E", "N", "U"));
    res_sgps = filter_disp_pts_to_vector(utilities.filter_to_meas_type(resid, "survey"), ("E", "N"));
    res_lev = filter_disp_pts_to_vector(utilities.filter_to_meas_type(resid, "leveling"), ("U",));
    res_tg = filter_disp_pts_to_vector(utilities.filter_to_meas_type(resid, "tide_gage"), ("U",));
    resid_data = res_cgps[0] + res_sgps[0] + res_lev[0] + res_tg[0];
    resid_sigma = res_cgps[1] + res_sgps[1] + res_lev[1] + res_tg[1];
    rms_mm, reported_chi2 = L2_on_vector(resid_data, resid_sigma);
    return rms_mm, reported_chi2;


def obs_vs_model_L2_horiz(obs_disp_points, model_disp_points):
    """
    L2 norm on observed vs modeled disp_point objects. Unpacks into vectors, then computes L2 norm.
    """
    resid = utilities.subtract_disp_points(obs_disp_points, model_disp_points);   # make residual points
    res_cgps = filter_disp_pts_to_vector(utilities.filter_to_meas_type(resid, "continuous"), ("E", "N"));
    res_sgps = filter_disp_pts_to_vector(utilities.filter_to_meas_type(resid, "survey"), ("E", "N"));
    resid_data = res_cgps[0] + res_sgps[0];
    resid_sigma = res_cgps[1] + res_sgps[1];
    rms_mm, reported_chi2 = L2_on_vector(resid_data, resid_sigma);
    return rms_mm, reported_chi2;


def obs_vs_model_L2_vertical(obs_disp_points, model_disp_points):
    """
    L2 norm on observed vs modeled disp_point objects. Unpacks into vectors, then computes L2 norm.
    """
    resid = utilities.subtract_disp_points(obs_disp_points, model_disp_points);   # make residual points
    res_cgps = filter_disp_pts_to_vector(utilities.filter_to_meas_type(resid, "continuous"), ("U",));
    res_lev = filter_disp_pts_to_vector(utilities.filter_to_meas_type(resid, "leveling"), ("U",));
    res_tg = filter_disp_pts_to_vector(utilities.filter_to_meas_type(resid, "tide_gage"), ("U",));
    resid_data = res_cgps[0] + res_lev[0] + res_tg[0];
    resid_sigma = res_cgps[1] + res_lev[1] + res_tg[1];
    rms_mm, reported_chi2 = L2_on_vector(resid_data, resid_sigma);
    return rms_mm, reported_chi2;


def obs_vs_model_L2_misfit(obs_disp_points, model_disp_points):
    """
    Implement one definition of model misfit: L2 norm
    """
    resid_pts = utilities.subtract_disp_points(obs_disp_points, model_disp_points);  # make residual points
    all_resids_m, all_sigmas = filter_disp_pts_to_vector(resid_pts, ("E", "N", "U"));
    horiz_resids_m, horiz_sigmas = filter_disp_pts_to_vector(resid_pts, ("E", "N"));
    all_L2_norm = np.sqrt(np.sum(np.square(all_resids_m)));  # L2 norm
    all_misfits_norm = np.divide(np.square(all_resids_m), np.square(all_sigmas));  # chi2 as a vector
    avg_misfit_norm = np.nanmean(all_misfits_norm);
    horiz_L2_norm = np.sqrt(np.sum(np.square(horiz_resids_m)));  # L2 norm
    horiz_misfits_norm = np.divide(np.square(horiz_resids_m), np.square(horiz_sigmas));  # chi2 as a vector
    avg_horiz_norm = np.nanmean(horiz_misfits_norm);
    return all_L2_norm, avg_misfit_norm, horiz_L2_norm, avg_horiz_norm;


def obs_vs_model_L1_misfit(obs_disp_points, model_disp_points):
    """
    Implement one definition of model misfit: L1 norm
    """
    resid_pts = utilities.subtract_disp_points(obs_disp_points, model_disp_points);  # make residual points
    all_misfits_m, all_sigmas = filter_disp_pts_to_vector(resid_pts, ("E", "N", "U"));
    horiz_misfits_m, horiz_sigmas = filter_disp_pts_to_vector(resid_pts, ("E", "N"));
    avg_misfit_norm = np.nanmean(np.divide(np.abs(all_misfits_m), all_sigmas))
    avg_horiz_norm = np.nanmean(np.divide(np.abs(horiz_misfits_m), horiz_sigmas));
    return np.abs(np.nanmean(all_misfits_m)), avg_misfit_norm, np.nanmean(np.abs(horiz_misfits_m)), avg_horiz_norm;


def compute_rms(disp_points):
    """
    Compute the bare RMS value of all displacements within a list of disp_points
    """
    values, _sigmas = filter_disp_pts_to_vector(disp_points);
    return np.sqrt(np.mean(np.square(values)));
