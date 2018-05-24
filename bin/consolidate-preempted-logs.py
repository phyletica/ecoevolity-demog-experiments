#! /usr/bin/env python

import sys
import os
import re
import glob
import argparse

import project_util


batch_number_pattern = re.compile(r'batch(?P<batch_number>\d+)')
sim_number_pattern = re.compile(r'-sim-(?P<sim_number>\d+)-')
run_number_pattern = re.compile(r'-run-(?P<sim_number>\d+)\.log')

def line_count(path):
    count = 0
    with open(path) as stream:
        for line in stream:
            count += 1
    return count

def get_run_number(log_path):
    run_number_matches = run_number_pattern.findall(log_path)
    assert len(run_number_matches) == 1
    run_number_str = run_number_matches[0]
    return int(run_number_str)

def consolidate_preempted_logs(
        target_run_number = 1,
        number_of_samples = 1501,
        batch_dir_name = None):
    number_of_lines = number_of_samples + 1
    val_sim_dirs = glob.glob(os.path.join(project_util.VAL_DIR, '0*'))
    for val_sim_dir in sorted(val_sim_dirs):
        sim_name = os.path.basename(val_sim_dir)
        batch_dirs = glob.glob(os.path.join(val_sim_dir, "batch*"))
        for batch_dir in sorted(batch_dirs):
            if batch_dir_name and (os.path.basename(batch_dir) != batch_dir_name):
                sys.stderr.write("Skipping {0}\n".format(batch_dir))
                continue
            batch_number_matches = batch_number_pattern.findall(batch_dir)
            assert len(batch_number_matches) == 1
            batch_number_str = batch_number_matches[0]
            batch_number = int(batch_number_str)

            sh_paths = glob.glob(os.path.join(batch_dir,
                    "*simcoevolity-sim-*-config-run-{0}-qsub.sh".format(
                            target_run_number)))
            if not sh_paths:
                sys.stderr.write("WARNING: No qsub files found for\n"
                        "    Simulation: {0}\n"
                        "    Batch: {1}\n"
                        "    Target run: {2}\n    Skipping!!\n".format(
                            sim_name,
                            batch_number,
                            target_run_number))
                continue
            for sh_path in sorted(sh_paths):
                posterior_path = sh_path.replace(
                        "-run-{0}-qsub.sh".format(target_run_number),
                        "-state-run-{0}.log".format(target_run_number))
                if not os.path.exists(posterior_path):
                    sys.stderr.write("WARNING: Missing log: {0}\n".format(posterior_path))
                    sys.stdout.write("{0}\n".format(sh_path))
                    continue
                sim_number_matches = sim_number_pattern.findall(posterior_path)
                assert len(sim_number_matches) == 1
                sim_number_str = sim_number_matches[0]
                sim_number = int(sim_number_str)
                posterior_file = os.path.basename(posterior_path)
                prefix = posterior_file.split("-sim-")[0]
                gp = os.path.join(batch_dir,
                        "{0}-sim-{1}-config-state-run-{2}.log*".format(
                                prefix,
                                sim_number_str,
                                target_run_number))
                target_state_log_paths = glob.glob(gp)
                assert (len(target_state_log_paths) == 1), (
                        "Multiple matches to {0!r}".format(gp))
                target_state_log_path = target_state_log_paths[0]
                gp = os.path.join(batch_dir,
                        "{0}-sim-{1}-config-operator-run-{2}.log*".format(
                                prefix,
                                sim_number_str,
                                target_run_number))
                target_op_log_paths = glob.glob(gp)
                assert (len(target_op_log_paths) == 1), (
                        "Multiple matches to {0!r}".format(gp))
                target_op_log_path = target_op_log_paths[0]
                state_log_path_pattern = os.path.join(batch_dir,
                        "{0}-sim-{1}-config-state-run-*.log*".format(
                                prefix,
                                sim_number_str))
                state_log_paths = glob.glob(state_log_path_pattern)
                op_log_path_pattern = os.path.join(batch_dir,
                        "{0}-sim-{1}-config-operator-run-*.log*".format(
                                prefix,
                                sim_number_str))
                op_log_paths = glob.glob(op_log_path_pattern)
                assert (len(state_log_paths) == len(op_log_paths)), (
                        "{0} matches for {1!r} and {2} for {3!r}".format(
                            len(state_log_paths),
                            state_log_path_pattern,
                            len(op_log_paths),
                            op_log_path_pattern))
                assert (target_state_log_path in state_log_paths), (
                        "Target {0!r} not in matches".format(
                                target_state_log_path))
                assert (target_op_log_path in op_log_paths), (
                        "Target {0!r} not in matches".format(
                                target_op_log_path))
                run_numbers = sorted(get_run_number(p) for p in state_log_paths)
                assert (run_numbers == sorted(get_run_number(p) for p in op_log_paths))
                extra_run_numbers = [rn for rn in run_numbers if rn > target_run_number]
                if len(extra_run_numbers) < 1:
                    if line_count(target_state_log_path) != number_of_lines:
                        sys.stderr.write(
                                "WARNING: Target log is incomplete, but there are no extra runs\n"
                                "    Simulation: {0}\n"
                                "    Batch: {1}\n"
                                "    Rep: {2}\n"
                                "    Target run: {3}\n    Skipping!!\n".format(
                                        sim_name,
                                        batch_number,
                                        sim_number,
                                        target_run_number))
                        sys.stdout.write("{0}\n".format(sh_path))
                    continue
                else:
                    if line_count(target_state_log_path) >= number_of_lines:
                        sys.stderr.write(
                                "WARNING: Target log is complete, but there are extra runs\n"
                                "    Simulation: {0}\n"
                                "    Batch: {1}\n"
                                "    Rep: {2}\n"
                                "    Target run: {3}\n    Skipping!!\n".format(
                                        sim_name,
                                        batch_number,
                                        sim_number,
                                        target_run_number))
                        sys.stdout.write("{0}\n".format(sh_path))
                        continue
                    completed_run_number = extra_run_numbers.pop(-1)
                    completed_state_log_pattern = os.path.join(batch_dir,
                            "{0}-sim-{1}-config-state-run-{2}.log*".format(
                                    prefix,
                                    sim_number_str,
                                    completed_run_number))
                    completed_state_log_paths = glob.glob(completed_state_log_pattern)
                    assert (len(completed_state_log_paths) == 1), (
                            "Multiple matches to complete state log {0!r}".format(
                                    completed_state_log_pattern))
                    completed_state_log_path = completed_state_log_paths[0]
                    completed_op_log_pattern = os.path.join(batch_dir,
                            "{0}-sim-{1}-config-operator-run-{2}.log*".format(
                                    prefix,
                                    sim_number_str,
                                    completed_run_number))
                    completed_op_log_paths = glob.glob(completed_op_log_pattern)
                    assert (len(completed_op_log_paths) == 1), (
                            "Multiple matches to complete op log {0!r}".format(
                                    completed_state_log_pattern))
                    completed_op_log_path = completed_op_log_paths[0]
                    if line_count(completed_state_log_path) != number_of_lines:
                        sys.stderr.write(
                                "WARNING: could not find completed log for\n"
                                "    Simulation: {0}\n"
                                "    Batch: {1}\n"
                                "    Rep: {2}\n"
                                "    Target run: {3}\n    Skipping!!\n".format(
                                        sim_name,
                                        batch_number,
                                        sim_number,
                                        target_run_number))
                        sys.stdout.write("{0}\n".format(sh_path))
                        continue
                    os.rename(completed_state_log_path, target_state_log_path)
                    os.rename(completed_op_log_path, target_op_log_path)
                    for n in extra_run_numbers:
                        sp = os.path.join(batch_dir,
                                "{0}-sim-{1}-config-state-run-{2}.log*".format(
                                        prefix,
                                        sim_number_str,
                                        n))
                        state_purge_paths = glob.glob(sp)
                        assert (len(state_purge_paths) == 1), (
                                "Multiple matches to incomplete state log {0!r}".format(
                                        sp))
                        state_purge_path = state_purge_paths[0]
                        op = os.path.join(batch_dir,
                                "{0}-sim-{1}-config-operator-run-{2}.log*".format(
                                        prefix,
                                        sim_number_str,
                                        n))
                        op_purge_paths = glob.glob(op)
                        assert (len(op_purge_paths) == 1), (
                                "Multiple matches to incomplete op log {0!r}".format(
                                        op))
                        op_purge_path = op_purge_paths[0]
                        os.remove(state_purge_path)
                        os.remove(op_purge_path)
                    

def main_cli(argv = sys.argv):
    parser = argparse.ArgumentParser()

    parser.add_argument('-r', '--run-number',
            action = 'store',
            type = int,
            default = 1,
            help = 'Target run number for consolidation.')
    parser.add_argument('-n', '--number-of-samples',
            action = 'store',
            type = int,
            default = 1501,
            help = ('Number of MCMC samples that should be found in the '
                    'completed log file of each analysis.'))
    parser.add_argument('-b', '--batch-dir',
            action = 'store',
            type = str,
            default = None,
            help = ('Batch directory name.'))

    if argv == sys.argv:
        args = parser.parse_args()
    else:
        args = parser.parse_args(argv)

    consolidate_preempted_logs(
            target_run_number = args.run_number,
            number_of_samples = args.number_of_samples,
            batch_dir_name = args.batch_dir)


if __name__ == "__main__":
    main_cli()
