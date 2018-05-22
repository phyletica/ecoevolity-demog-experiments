#! /usr/bin/env python

import sys
import os
import glob

from yaml import load, dump
try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

import project_util

def get_yaml_config(path):
    with open(path, 'r') as stream:
        config = load(stream, Loader = Loader)
    return config

def convert_for_variable_sites_only(config):
    if "global_comparison_settings" in config:
        config["global_comparison_settings"]["constant_sites_removed"] = True
    for c in config["comparisons"]:
        comparison = c["comparison"]
        comparison["constant_sites_removed"] = True
    return config

def estimating_population_sizes(config):
    try:
        global_estimate = int(config["global_comparison_settings"]["parameters"]["population_size"]["estimate"])
    except:
        global_estimate = 1
    comp_estimate = []
    for c in config["comparisons"]:
        comparison = c["comparison"]
        try:
            comp_estimate.append(comparison["parameters"]["population_size"]["estimate"])
        except:
            comp_estimate.append(global_estimate)
    return sum(comp_estimate) > 0

def update_time_size_rate_scalers(config):
    if not estimating_population_sizes(config):
        return config
    config["operator_settings"]["operators"]["TimeSizeRateScaler"]["weight"] = 3.0
    return config


def main_cli(argv = sys.argv):
    config_path_pattern = os.path.join(project_util.VAL_DIR, "*", "*", "simcoevolity-sim-*.yml") 

    for config_path in sorted(glob.glob(config_path_pattern)):
        dir_path = os.path.dirname(config_path)
        config_name = os.path.basename(config_path)
        new_config_name = "var-only-" + config_name
        new_config_path = os.path.join(dir_path, new_config_name)
        if os.path.exists(new_config_path):
            sys.stderr.write(
                    "Skipping config that already exists: {0}\n".format(
                            new_config_path))
            continue

        config = get_yaml_config(config_path)
        config = convert_for_variable_sites_only(config)
        config = update_time_size_rate_scalers(config)

        assert(not os.path.exists(new_config_path))
        with open(new_config_path, 'w') as stream:
            dump(config, stream, Dumper = Dumper, indent = 4)
    

if __name__ == "__main__":
    main_cli()
