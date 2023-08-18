#!/usr/bin/env python

import argparse
import logging
import os.path
import pathlib
import sys
import time as pytime
import numpy as np
import pandas as pd
import yaml

from ase.io.vasp import read_vasp, write_vasp_xdatcar
from jsonschema import validate


def acns_logger():
    log_formatter = logging.Formatter('ACNS %(levelname)s [%(asctime)s] | %(message)s')
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    return root_logger, log_formatter


def acns_welcome(name, autors, emails):
    logging.info("--------------------------------------------------")
    logging.info("      _/_/      _/_/_/  _/      _/    _/_/_/   ")
    logging.info("   _/    _/  _/        _/_/    _/  _/          ")
    logging.info("  _/_/_/_/  _/        _/  _/  _/    _/_/       ")
    logging.info(" _/    _/  _/        _/    _/_/        _/      ")
    logging.info("_/    _/    _/_/_/  _/      _/  _/_/_/         ")
    logging.info("")
    logging.info("---- Australian Centre for Neutron Scattering ----")
    logging.info("")
    logging.info(" ==============================================")
    logging.info(" " + name)
    logging.info(" Author(s): " + ",".join(autors))
    logging.info(" " + ",".join(emails))
    logging.info(" ==============================================")
    logging.info("")


schema = """
type: object
properties:
  structure:
    type: object
    properties:
        input_poscar:
           type: string
        output_xdatcar:
           type: string
  evolution:
    type: object
    properties:
      final_time: 
        type: number
        minimum: 0
      delta_time:
        type: number
        minimum: 0
  groups:
    type: array
    items:
      type: object
      properties:
          name: 
            type: string
          indices: 
            type: array
            items:
              type: integer
          frequency:
            type: array
            items:
              type: number
            minItems: 3
            maxItems: 3
          amplitude: 
            type: array
            items:
              type: number
            minItems: 3
            maxItems: 3
  supercell:
    type: object
    properties:
      dimensions: 
        type: array
        items:
          type: integer
        minItems: 3
        maxItems: 3
      parameters:
        type: array
        items:
          type: object
          properties:
            supercell_index: 
              type: array
              items:
                type: integer
              minItems: 3
              maxItems: 3
            groups:
              type: array
              items:
                type: object
                properties:
                  name:
                    type: string
                  phase: 
                    type: array
                    items:
                      type: number
                    minItems: 3
                    maxItems: 3
                  scale: 
                    type: array
                    items:
                      type: number
                    minItems: 3
                    maxItems: 3
"""


class Group:
    def __init__(self):
        self.amplitude = np.zeros(3)
        self.frequency = np.zeros(3)
        self.phase = np.zeros(3)
        self.scale = np.ones(3)
        self.indices = []

    def __str__(self):
        return f"amplitude:{self.amplitude}, frequency:{self.frequency}, phase:{self.phase}, scale:{self.scale}, indices:{self.indices} "

    def __repr__(self):
        return self.__str__()


if __name__ == "__main__":

    # store start time for benchmarking
    start_time = pd.to_datetime(pytime.time(), unit="s")

    # setup logger
    root_logger, log_formatter = acns_logger()

    # setup arguments
    epilog_text = "Australian Centre for Neutron Scattering - Scientific Computing"
    parser = argparse.ArgumentParser(description='region indexer for unknown chromosomes', epilog=epilog_text,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-s', '--silent', action='store_true', help='Starts in silent mode, no message will be output.')
    parser.add_argument('-d', '--debug', action='store_true', help='Shows debug info')
    parser.add_argument('-c', '--config', type=pathlib.Path, required=True, help='Path to configuration file')

    # parse arguments and set logger
    args = parser.parse_args()
    consoleHandler = logging.StreamHandler(sys.stdout)

    if args.debug:
        root_logger.setLevel(logging.DEBUG)

    if args.silent:
        consoleHandler.setLevel(logging.ERROR)

    consoleHandler.setFormatter(log_formatter)
    root_logger.addHandler(consoleHandler)

    acns_welcome("TrajMaker", ["Pablo Galaviz"], ["galavizp@ansto.gov.au"])

    logging.info("Configuration file: %s" % args.config)

    if os.path.exists(args.config) and not os.path.isdir(args.config):

        if args.config.name.endswith(".yaml"):
            with open(args.config.name, 'r') as file_yaml:
                config_yaml = yaml.safe_load(file_yaml)
        schema_yaml = yaml.safe_load(schema)
        try:
            validate(config_yaml, schema_yaml)
        except:
            logging.error("Error in configuration file: %s.", )
            exit(-1)

        structure_options = config_yaml["structure"]

        ase_atoms = read_vasp(structure_options["input_poscar"])
        groups_options = config_yaml["groups"]
        total_atoms = len(ase_atoms.get_positions())

        supercell = config_yaml["supercell"]
        dimensions = supercell["dimensions"].copy()
        for d in dimensions:
            if d <= 0 or type(d) != int:
                logging.error("Supercell dimensions should be and integer >= 1")
                exit(-1)
        supercell_dict = dict()
        for item in supercell["parameters"]:
            key = "-".join("%d" % x for x in item["supercell_index"])
            group_dict = dict()
            for group in item["groups"]:
                name = group.pop("name")
                group_dict[name] = group
            supercell_dict[key] = group_dict

        index_mapping = np.arange(np.prod(dimensions) * total_atoms)
        dimensions += [-1]
        index_mapping = index_mapping.reshape(dimensions)

        group_list = list()
        for group in groups_options:
            indices = group["indices"]
            if np.min(indices) < 1 or np.max(indices) > total_atoms:
                logging.error("Atom index out of bounds for group %s", group)
                logging.error("Index should be between 1 and %d", total_atoms)
                exit(-1)
            else:
                amplitude = np.array(group["amplitude"])
                frequency = np.array(group["frequency"])
                for a in np.arange(dimensions[0]):
                    for b in np.arange(dimensions[1]):
                        for c in np.arange(dimensions[2]):
                            cell_key = "%d-%d-%d" % (a + 1, b + 1, c + 1)
                            new_group = Group()
                            new_group.indices = [index_mapping[a, b, c, i - 1] for i in indices]
                            new_group.amplitude = np.array(amplitude)
                            new_group.frequency = np.array(frequency)
                            if cell_key in supercell_dict.keys():
                                sc = supercell_dict[cell_key]
                                group_name = group["name"]
                                if group_name in sc.keys():
                                    new_group.phase = np.array(sc[group_name]["phase"])
                                    new_group.scale = np.array(sc[group_name]["scale"])
                            group_list.append(new_group)

        evolution_options = config_yaml["evolution"]
        final_time = evolution_options["final_time"]
        delta_time = evolution_options["delta_time"]
        cell_list = list()
        for time in np.arange(start=0, stop=final_time, step=delta_time):
            logging.debug("calculating time %f", time)
            new_ase_atoms = ase_atoms.copy() * supercell["dimensions"]
            positions = new_ase_atoms.get_positions()
            for group in group_list:
                positions[group.indices] += group.scale * group.amplitude * np.sin(
                    2 * np.pi * (group.frequency * time - group.phase))
            new_ase_atoms.set_positions(positions)
            cell_list.append(new_ase_atoms)
        logging.info("Saving %d frames to %s", len(cell_list), structure_options["output_xdatcar"])
        write_vasp_xdatcar(structure_options["output_xdatcar"], cell_list, label='Trajectory')
    else:
        logging.error("Configuration file %s does not exists" % args.config)

    logging.info("Total computation time: %s", str(pd.to_datetime(pytime.time(), unit="s") - start_time))
